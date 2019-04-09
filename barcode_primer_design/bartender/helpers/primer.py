from typing import Optional, List

from Bio.SeqFeature import SeqFeature, FeatureLocation


class Amplicon:
    def __init__(self, sequence: str):
        self.sequence = sequence
        self.features: List[SeqFeature] = []
        self.primer_set_fwd: Optional[PrimerSet] = None
        self.primer_set_rev: Optional[PrimerSet] = None

    def add_feature(self, feature: SeqFeature):
        self.features.append(feature)

    def __len__(self):
        return len(self.sequence)


class Gene:
    def __init__(self, name: str):
        self.name = name
        self.amplicons: List[Amplicon] = []

    def append(self, amplicon: Amplicon):
        self.amplicons.append(amplicon)

    def remove(self, amplicon: Amplicon):
        self.amplicons.remove(amplicon)

    def __len__(self) -> int:
        return len(self.amplicons)


class ExcludedRegion(SeqFeature):
    def __init__(self, location: FeatureLocation):
        super().__init__(self, location, None)


class TargetRegion(SeqFeature):
    def __init__(self, location: FeatureLocation):
        super().__init__(self, location, None)


class Primer(SeqFeature):
    """
    start is a position in python style, i.e. as an index.
    """

    def __init__(
        self,
        sequence: str,
        start: int,
        length: Optional[int] = None,
        reverse: bool = False,
    ):
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
        super().__init__(self, location, strand=strand)

    def __len__(self) -> int:
        return len(self.sequence)


class PrimerSet:
    def __init__(self, name: str):
        self.name = name
        self.set: List[Primer] = []

    def __len__(self) -> int:
        return len(self.set)

    def append(self, primer: Primer):
        self.set.append(primer)
