from __future__ import print_function
import re
import shlex
from Bio import SeqIO
import subprocess
from helpers.p3_parser import P3Parser
from helpers.primer import *

class P3SeqResult:

    def __init__(self, spacing, interval):
        self.spacing = spacing
        self.interval = interval
        self.warning = ""
        self.error = ""

class P3Seq:

    def __init__(self, config, input_handle):
        self.config = config
        self.input_handle = input_handle


    def run(self, spacing_range, interval_range):
        output = dict()
        for record in SeqIO.parse(self.input_handle, "fasta"):
            out_seq = Gene(record.id)
            for i, spacing in enumerate(spacing_range):
                for j, interval in enumerate(interval_range):
                    sequence = str(record.seq)

                    ampl = Amplicon(sequence)
                    ampl.spacing = spacing
                    ampl.interval = interval

                    input_string = ""
                    input_string += "SEQUENCE_ID=" + record.id + "_" + spacing + "_" + interval + "\n"
                    input_string += "PRIMER_TASK=pick_sequencing_primers\n"
                    input_string += "PRIMER_SEQUENCING_SPACING=" + spacing + "\n"
                    input_string += "PRIMER_SEQUENCING_INTERVAL=" + interval + "\n"

                    if sequence.find("<") >= 0 and sequence.find(">") >= 0:
                        input_string += "SEQUENCE_EXCLUDED_REGION="
                        spl_sequence = re.split("\<|\>",sequence.replace("[","").replace("]",""))
                        for i in xrange(0,len(spl_sequence)-1,2):
                            start = 0
                            for j in xrange(0,i+1):
                                start += len(spl_sequence[j])
                            input_string += str(start+1) + "," + str(len(spl_sequence[i+1])) + " "
                            ampl.add_feature(ExcludedRegion(FeatureLocation(start+1, start + len(spl_sequence[i+1]))))
                        input_string += "\n"


                    sequence = sequence.replace("<", "")
                    sequence = sequence.replace(">", "")

                    if sequence.find("[") >= 0 and sequence.find("]") >= 0:
                        input_string += "SEQUENCE_TARGET="
                        spl_sequence = re.split("\[|\]", sequence)
                        for i in xrange(0,len(spl_sequence)-1, 2):
                            start = 0
                            for j in xrange(0,i+1):
                                start += len(spl_sequence[j])
                            input_string += str(start+1) + "," + str(len(spl_sequence[i+1])) + " "
                            ampl.add_feature(TargetRegion(FeatureLocation(start+1, start + len(spl_sequence[i+1]))))
                        input_string += "\n"


                    sequence = sequence.replace("]", "")
                    sequence = sequence.replace("[", "")

                    ampl.sequence = sequence

                    input_string += "SEQUENCE_TEMPLATE=" + sequence + "\n"
                    input_string += "P3_FILE_FLAG=0\n"
                    input_string += "PRIMER_EXPLAIN_FLAG=1\n"
                    input_string += "PRIMER_NUM_RETURN=20\n"
                    input_string += "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + self.config.p3_thermo_path + "\n="

                    cmd = self.config.p3_path + " -p3_settings_file=" + self.config.p3_config_path
                    args = shlex.split(cmd)
                    p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                    p3_output = p.communicate(input_string)[0].strip()

                    primer_set_fwd = PrimerSet(record.id + "_" + spacing + "_" + interval + "_left")
                    primer_set_rev = PrimerSet(record.id + "_" + spacing + "_" + interval + "_right")

                    m = re.search('(?<=PRIMER_ERROR=).+', p3_output)
                    ampl.error = ""
                    ampl.warning = ""
                    if m is not None:
                        ampl.error = m.group(0)
                    elif p3_output != "":
                        mw = re.search('(?<=PRIMER_WARNING=).+', p3_output)
                        if mw is not None:
                            ampl.warning += mw.group(0)
                        P3Parser.parse_p3seq_information(primer_set_fwd, primer_set_rev, p3_output, ampl, self.config.p3_config_path)

                    ampl.primer_set_fwd = primer_set_fwd
                    ampl.primer_set_rev = primer_set_rev
                    
                    out_seq.amplicons.append(ampl)

            output[record.id] = out_seq
        return output

