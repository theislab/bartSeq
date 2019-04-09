import itertools

from .filters import gc_content, repeats, similarity


def main():
    seq_length = 10
    sequences = []
    for it in itertools.product("ACTG", repeat=seq_length):
        sequences.append("".join(it))
    print("Number of initial sequences:", len(sequences))
    sequences = gc_content(sequences, 5)
    print("Number of GC filtered sequences:", len(sequences))
    sequences = repeats(sequences, 2, 3)
    print("Number of repeat filtered sequences:", len(sequences))
    similarity(sequences)
