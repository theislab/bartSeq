from collections.abc import MutableSequence
from typing import List, Iterator, Iterable

from .primer import Primer


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


class PrimerPairSet(MutableSequence, Iterable[PrimerPair]):
    def __init__(self, name: str, items: Iterable[PrimerPair] = ()):
        self.name = name
        self._set: List[PrimerPair] = list(items)

    def insert(self, i: int, o: PrimerPair) -> None:
        self._set.insert(i, o)

    def __getitem__(self, i: int) -> PrimerPair:
        return self._set[i]

    def __setitem__(self, i: int, o: PrimerPair) -> None:
        self._set[i] = o

    def __delitem__(self, i: int) -> None:
        del self._set[i]

    def __iter__(self) -> Iterator[PrimerPair]:
        return iter(self._set)

    def __len__(self) -> int:
        return len(self._set)
