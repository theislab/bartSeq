import re
import sys
import shlex
import subprocess

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

from ..helpers.p3_parser import P3Parser
from ..helpers.primerpair import PrimerPair, PrimerPairSet
from ..helpers.primer import Primer, Amplicon, Gene, ExcludedRegion, TargetRegion


class PrimerPredictor:
    def __init__(self, config, input_handle, predefined_handle):
        self.config = config
        self.input_handle = input_handle
        self.predefined_handle = predefined_handle

    def parse_predefined_pairs(self, predefined_sets):
        for record in SeqIO.parse(self.predefined_handle, "fasta"):

            cur_id = record.id.split("_")[0]
            seq = str(record.seq)

            if seq.find("&") == -1:
                print("Please specify fwd and rev primer sequences by separating them with '&' for the predefined "
                      f"primer {record.id}.")
                sys.exit(1)

            seqs = seq.split("&")
            if len(seqs) != 2:
                print(("Exactly two primer sequences (fwd&rev) have to provided for the predefined "
                      f"primer {record.id}."))
                sys.exit(1)

            pair_ind = 0
            if cur_id in predefined_sets:
                act_set = predefined_sets[cur_id]
                for pair in act_set.set:
                    ind = int(pair.name.split("_")[1])
                    if ind > pair_ind:
                        pair_ind = ind
                pair_ind += 1
            else:
                act_set = predefined_sets[cur_id] = PrimerPairSet(cur_id)

            act_set.set.append(PrimerPair(
                Primer(seqs[0], 0),
                Primer(seqs[1], 0, reverse=True),
                cur_id + "_" + str(pair_ind),
                True,
            ))

    def predict_primer_set(self):
        predefined_sets = dict()
        if self.predefined_handle is not None:
            self.parse_predefined_pairs(predefined_sets)

        out_genes = []
        for record in SeqIO.parse(self.input_handle, "fasta"):
            gene = Gene(record.id)
            sequence = str(record.seq)
            for sel_sequence in re.split(r'//', sequence):
                s = re.sub(r"[\[\]<>]", "", sel_sequence)
                amplicon = Amplicon(s)

                if record.id in predefined_sets:
                    amplicon.primer_set = predefined_sets[record.id]
                    gene.append(amplicon)
                    del predefined_sets[record.id]
                    continue

                input_string = ""
                input_string += "SEQUENCE_ID=" + record.id + "\n"
                input_string += "SEQUENCE_TEMPLATE=" + s + "\n"

                if sel_sequence.find("<") >= 0 and sel_sequence.find(">") >= 0:
                    input_string += "SEQUENCE_EXCLUDED_REGION="
                    spl_sequence = re.split(r"[<>]", sel_sequence.replace("[", "").replace("]", ""))
                    for i in range(0, len(spl_sequence) - 1, 2):
                        start = 0
                        for j in range(0, i + 1):
                            start += len(spl_sequence[j])
                        input_string += str(start + 1) + "," + str(len(spl_sequence[i + 1])) + " "
                        amplicon.add_feature(
                            ExcludedRegion(FeatureLocation(start + 1, start + len(spl_sequence[i + 1]))))
                    input_string += "\n"

                sel_sequence = sel_sequence.replace("<", "")
                sel_sequence = sel_sequence.replace(">", "")

                if sel_sequence.find("[") >= 0 and sel_sequence.find("]") >= 0:
                    input_string += "SEQUENCE_TARGET="
                    spl_sequence = re.split(r"[\[\]]", sel_sequence)
                    for i in range(0, len(spl_sequence) - 1, 2):
                        start = 0
                        for j in range(0, i + 1):
                            start += len(spl_sequence[j])
                        input_string += str(start + 1) + "," + str(len(spl_sequence[i + 1])) + " "
                        amplicon.add_feature(TargetRegion(FeatureLocation(start + 1, start + len(spl_sequence[i + 1]))))
                    input_string += "\n"

                input_string += "P3_FILE_FLAG=0\n"
                input_string += "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + self.config.p3_thermo_path + "\n="

                cmd = self.config.p3_path + " -p3_settings_file=" + self.config.p3_config_path
                args = shlex.split(cmd)

                print(input_string)

                p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

                p3_output = p.communicate(input_string)[0].strip()
                print(p3_output)

                m = re.search(r'(?<=PRIMER_ERROR=)\w+', p3_output)
                if m is not None:
                    raise Exception(f"Error for sequence (Probably no primer found in region): {record.id}: {m.group(0)}\n Start NEW Primerprediction.")

                primer_set = PrimerPairSet(record.id)
                P3Parser.parse_p3_information(primer_set, p3_output)

                if len(primer_set) == 0:
                    print(f"WARNING: No primer found for {record.id} sequence {i+1}.")
                    continue

                amplicon.primer_set = primer_set
                gene.append(amplicon)

            if len(gene) == 0:
                raise Exception(f"No primer found for {gene.name}. Consider less restrictive Primer3 settings.")
            out_genes.append(gene)
        for key in predefined_sets:
            print("WARNING: No input sequence could be found for the predefined primer", key)

        return out_genes
