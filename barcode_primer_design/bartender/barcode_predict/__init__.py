import itertools
from logging import getLogger

from .filters import gc_content, repeats, similarity


log = getLogger(__name__)


def main():
    seq_length = 10
    sequences = []
    for it in itertools.product("ACTG", repeat=seq_length):
        sequences.append("".join(it))
    log.info("Number of initial sequences: %s", len(sequences))
    sequences = gc_content(sequences, 5)
    log.info("Number of GC filtered sequences: %s", len(sequences))
    sequences = repeats(sequences, 2, 3)
    log.info("Number of repeat filtered sequences: %s", len(sequences))
    similarity(sequences)
