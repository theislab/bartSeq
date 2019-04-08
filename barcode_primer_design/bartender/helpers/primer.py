from Bio.SeqFeature import SeqFeature, FeatureLocation

class Gene:
    def __init__(self, name):
        self.name = name
        self.amplicons = []

    def append(self, amplicon):
        self.amplicons.append(amplicon)

    def remove(self, amplicon):
        self.amplicons.remove(amplicon)

    def __len__(self):
        return len(self.amplicons)

class Amplicon:
    def __init__(self, sequence):
        self.sequence = sequence
        self.features = []

    def add_feature(self, feature):
        self.features.append(feature)

    def __len__(self):
        return len(self.sequence)


class ExcludedRegion(SeqFeature):
    def __init__(self, location):
        SeqFeature.__init__(self, location, None)

class TargetRegion(SeqFeature):
    def __init__(self, location):
        SeqFeature.__init__(self, location, None)


class Primer(SeqFeature):
    """
    start is a position in python style, i.e. as an index.
    """
    def __init__(self, sequence, start, length=None, reverse=False):
        self.sequence = sequence
        if length is None:
            length = len(sequence)
        if not reverse:
            # Python-style indexing, think range(start, start+length) or seq[start:start+length]
            location = FeatureLocation(start, start + length)
            strand = +1
        else:
            # The same applies here: We are handling boundaries
            location = FeatureLocation(start - length, start)
            strand = -1
        SeqFeature.__init__(self, location, strand=strand)

    def __len__(self):
        return len(self.sequence)

class PrimerSet:

    def __init__(self, name):
        self.name = name
        self.set = []

    def __len__(self):
        return len(self.set)

    def append(self, primer):
        self.set.append(primer)
