from __future__ import print_function
from functools import partial
import multiprocessing
from collections import deque
from multiprocessing.pool import ThreadPool
import shlex
import os
import subprocess
from werkzeug.utils import secure_filename
from helpers.primerpair import PrimerPair


class Blaster:

    def __init__(self, config):
        self.config = config


    def run_task(self, input):
        cmd = self.config.blast_path + " -p blastn -r 1 -G -5 -E -5 -m 8 -d " + os.path.join(self.config.blast_dbpath, self.config.blast_dbname)
        print(cmd)
        print("BLASTing primer " + str.split(str(input), "_")[0])
        p = subprocess.Popen(shlex.split(cmd), stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        r_string = p.communicate(str(input))[0].strip()
        return r_string

    def blast_primer_set(self, sequence_sets, linkers):

        blast_strings = []
        for gene in sequence_sets:
            for amplicon in gene.amplicons:
                blast_string = ""
                for pair in amplicon.primer_set.set:
                    blast_string += ">" + pair.name + "_fwd\n" + linkers[0] + pair.fwd.sequence + "\n\n"
                    blast_string += ">" + pair.name + "_rev\n" + linkers[1] + pair.rev.sequence + "\n\n"
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
                blat_out_string = strings[index]
                if blat_out_string == "":
                    pair.fwd.blast_hits = 0
                    pair.rev.blast_hits = 0
                else:
                    blat_out_string = blat_out_string.split("\n")

                    for line in blat_out_string:
                        act_result = line.strip().split("\t")
                        if float(act_result[10]) < 0.1:
                            blast_hits.append(act_result[0])

                    for i, pair in enumerate(amplicon.primer_set.set):
                        pair.fwd.blast_hits = blast_hits.count(pair.name + "_fwd")
                        pair.rev.blast_hits = blast_hits.count(pair.name + "_rev")

                amplicon.primer_set.set = [pair for pair in amplicon.primer_set.set if pair.fwd.blast_hits <= self.config.blast_max_hits and pair.rev.blast_hits <= self.config.blast_max_hits]
                print(gene.name + " " + str(len(amplicon.primer_set.set)))

                if len(amplicon.primer_set.set) == 0:
                    print("WARNING: At least one amplicon was removed from gene " + gene.name + " due to too many BLAST hits.\n")

                index += 1

            gene.amplicons = [amplicon for amplicon in gene.amplicons if len(amplicon.primer_set.set) != 0]
            if len(gene.amplicons) == 0:
                print("Remove " + gene.name)
                raise Exception("No primers left for " + gene.name + ". Consider less restrictive BLAST settings.")

