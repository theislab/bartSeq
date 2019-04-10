from collections.abc import MutableSequence
from typing import Optional, List, Iterator, Iterable

from Bio.SeqFeature import SeqFeature, FeatureLocation


class Amplicon:
    def __init__(self, sequence: str):
        from ..helpers.primerpair import PrimerPairSet

        self.sequence = sequence
        self.features: List[SeqFeature] = []
        self.primer_set: Optional[PrimerPairSet] = None
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
        super().__init__(location, None)


class TargetRegion(SeqFeature):
    def __init__(self, location: FeatureLocation):
        super().__init__(location, None)


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
        super().__init__(location, strand=strand)

    def __len__(self) -> int:
        return len(self.sequence)


class PrimerSet(MutableSequence, Iterable[Primer]):
    def __init__(self, name: str, items: Iterable[Primer] = ()):
        self.name = name
        self._set: List[Primer] = list(items)

    def insert(self, i: int, o: Primer) -> None:
        self._set.insert(i, o)

    def __getitem__(self, i: int) -> Primer:
        return self._set[i]

    def __setitem__(self, i: int, o: Primer) -> None:
        self._set[i] = o

    def __delitem__(self, i: int) -> None:
        del self._set[i]

    def __iter__(self) -> Iterator[Primer]:
        return iter(self._set)

    def __len__(self) -> int:
        return len(self._set)
