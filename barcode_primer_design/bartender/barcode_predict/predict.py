from __future__ import print_function
import itertools
from multiprocessing import freeze_support

from barcode_predict.filters import *


def main():
    seq_length = 10
    sequences = []
    for it in itertools.product("ACTG", repeat=seq_length):
        sequences.append("".join(it))
    print("Number of initial sequences: ", len(sequences))
    sequences = SequenceFilters.gc_content(sequences, 5)
    print("Number of GC filtered sequences: ", len(sequences))
    sequences = SequenceFilters.repeats(sequences, 2, 3)
    print("Number of repeat filtered sequences: ", len(sequences))
    SequenceFilters.similarity(sequences)

    # linker_sequences = []
    # for it in itertools.product("ACTG", repeat=10):
    #     linker_sequences.append("".join(it))
    # linker_sequences = SequenceFilters.gc_content(linker_sequences, 6)
    # linker_sequences = SequenceFilters.repeats(linker_sequences, 2, 3)
    # linker_sequences = random.sample(linker_sequences, 1000)
    # SequenceFilters.similarity(linker_sequences)


if __name__ == "__main__":
    freeze_support()
    main()