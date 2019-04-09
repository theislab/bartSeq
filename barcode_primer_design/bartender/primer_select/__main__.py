import argparse
import os
from contextlib import suppress

from .ps_configuration import PsConfiguration
from . import predict_primerset, optimize, output


parser = argparse.ArgumentParser(description="Run the PrimerSelect pipeline.")

parser.add_argument(
    "input",
    help="Input file containing the sequences in FASTA format. "
    "The FASTA headers indicate the sequence ID and have to be unique.",
    type=str,
)
parser.add_argument(
    "-predefined",
    dest="predefined",
    help="Input file containing sequences in FASTA format for predefined primer pairs. "
    "The primer pair sequences have to provided as 'fwdseq&revseq'. "
    "The FASTA header of a given primer has to be specified according to the "
    "corresponding input sequence ID. "
    "If you want to specify more than one primer pair per input "
    "sequence, please add '_0', '_1' ... to the sequence ID.",
    type=str,
    default="",
)

args = parser.parse_args()

if not os.path.isfile("../config.cfg"):
    PsConfiguration.write_standard_config("../config.cfg")

with open("../config.cfg", "r") as config_handle:
    config = PsConfiguration.read_config(config_handle)

linkers = ["ATGCGCATTC", "AGCGTAACCT"]

predef_handle = open(args.predefined, "r") if args.predefined else suppress()
with predef_handle, open(args.input, "r") as input_handle:
    primer_sets = predict_primerset(config, input_handle, predef_handle, linkers)

print("Blast done")
opt_result = optimize(config, primer_sets, linkers)
output_ = output(opt_result, primer_sets)
print(output_)
