from __future__ import print_function
import ConfigParser
import StringIO
import re
from helpers.primer import Primer


class PrimerPairSet:
    def __init__(self, name):
        self.name = name
        self.set = []

    def __len__(self):
        return len(self.set)

    def append(self, primer_pair):
        self.set.append(primer_pair)


class PrimerPair:
    def __init__(self, fwd, rev, name, predefined=False):
        self.fwd = fwd
        self.rev = rev
        self.name = name
        self.predefined = predefined

    def fwd_seq(self):
        return self.fwd.sequence

    def rev_seq(self):
        return self.rev.sequence