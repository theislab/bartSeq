from .primer import Primer


class PrimerPairSet:
    def __init__(self, name):
        self.name = name
        self.set = []

    def __len__(self):
        return len(self.set)

    def append(self, primer_pair):
        self.set.append(primer_pair)


class PrimerPair:
    def __init__(self, fwd: Primer, rev: Primer, name: str, predefined: bool = False):
        self.fwd = fwd
        self.rev = rev
        self.name = name
        self.predefined = predefined

    @property
    def fwd_seq(self) -> str:
        return self.fwd.sequence

    @property
    def rev_seq(self) -> str:
        return self.rev.sequence
