from .blaster import Blaster
from .optimizer import Optimizer
from .primer_predictor import PrimerPredictor
from .rnacofolder import Cofolder


class PrimerSelect:
    @staticmethod
    def output(arrangements, sequence_set):
        output_string = ""

        for i, run in enumerate(arrangements[0:5]):
            output_string += "Rank " + str(i + 1) + ": Sum MFE= -" + str(run[0]) + "\n"

            v = run[1]
            w = run[2]
            for j, seq in enumerate(sequence_set):
                amplicon = w[j]
                pset = seq.amplicons[amplicon].primer_set
                pair = pset.set[v[j]]
                output_string += (
                    pset.name
                    + " (sequence "
                    + str(amplicon + 1)
                    + ")"
                    + "\tfwd: "
                    + pair.fwd.sequence.upper()
                    + "\trev: "
                    + pair.rev.sequence.upper()
                    + "\tBLAST hits: "
                    + str(pair.fwd.blast_hits)
                    + " / "
                    + str(pair.rev.blast_hits)
                    + "\n"
                )
            output_string += "------------------\n"
        return output_string

    @staticmethod
    def optimize(config, sequence_set, linkers):
        cofolder = Cofolder(config)
        cofolder.cofold(sequence_set, linkers)
        optimizer = Optimizer(config, sequence_set)
        opt_result = optimizer.optimize()
        return opt_result

    @staticmethod
    def predict_primerset(config, input_handle, predefined_handle, linkers):
        primer_predictor = PrimerPredictor(config, input_handle, predefined_handle)

        primer_sets = primer_predictor.predict_primer_set()
        blaster = Blaster(config)
        blaster.blast_primer_set(primer_sets, linkers)
        return primer_sets
