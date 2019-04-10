from multiprocessing.pool import ThreadPool
import os
from typing import List, Iterable

from ..helpers.primerpair import PrimerPairSet
from ..helpers.primer import Gene
from ..helpers import run_and_feed


class Blaster:
    def __init__(self, config):
        self.config = config

    def run_task(self, input: str):
        print(f"BLASTing primer {input.split('_')[0]}")
        p = run_and_feed(
            self.config.blast_path,
            reward=1,
            gapopen=5,  # Steffen had -5 for both here?
            gapextend=5,
            outfmt=6,  # tabular
            db=os.path.join(self.config.blast_dbpath, self.config.blast_dbname),
            _input_str=input,
            _long_arg_prefix="-",
        )
        return p.stdout.strip()

    def blast_primer_set(self, sequence_sets: Iterable[Gene], linkers):
        blast_strings: List[str] = []
        for gene in sequence_sets:
            for amplicon in gene.amplicons:
                blast_string = ""
                for pair in amplicon.primer_set:
                    blast_string += (
                        f">{pair.name}_fwd\n{linkers[0]}{pair.fwd.sequence}\n\n"
                    )
                    blast_string += (
                        f">{pair.name}_rev\n{linkers[1]}{pair.rev.sequence}\n\n"
                    )
                blast_strings.append(blast_string)

        pool = ThreadPool(self.config.max_threads)
        results = []
        for s in blast_strings:
            results.append(pool.apply_async(self.run_task, (s,)))
        pool.close()
        pool.join()

        strings = []
        for r in results:
            strings.append(r.get())

        index = 0
        for gene in sequence_sets:
            for amplicon in gene.amplicons:
                blast_hits = []
                blast_out_string = strings[index]
                if blast_out_string == "":
                    for pair in amplicon.primer_set:
                        pair.fwd.blast_hits = 0
                        pair.rev.blast_hits = 0
                else:
                    for line in blast_out_string.split("\n"):
                        act_result = line.strip().split("\t")
                        if float(act_result[10]) < 0.1:
                            blast_hits.append(act_result[0])

                    for pair in amplicon.primer_set:
                        pair.fwd.blast_hits = blast_hits.count(pair.name + "_fwd")
                        pair.rev.blast_hits = blast_hits.count(pair.name + "_rev")

                amplicon.primer_set = PrimerPairSet(
                    amplicon.primer_set.name,
                    [
                        pair
                        for pair in amplicon.primer_set
                        if pair.fwd.blast_hits <= self.config.blast_max_hits
                        and pair.rev.blast_hits <= self.config.blast_max_hits
                    ],
                )
                print(gene.name, len(amplicon.primer_set))

                if len(amplicon.primer_set) == 0:
                    print(
                        "WARNING: At least one amplicon was removed from"
                        f"gene {gene.name} due to too many BLAST hits.\n"
                    )

                index += 1

            gene.amplicons = [
                amplicon for amplicon in gene.amplicons if len(amplicon.primer_set) != 0
            ]
            if len(gene.amplicons) == 0:
                print("Remove", gene.name)
                raise Exception(
                    "No primers left for "
                    + gene.name
                    + ". Consider less restrictive BLAST settings."
                )
