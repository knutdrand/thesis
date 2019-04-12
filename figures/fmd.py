import numpy as np


class FMIndex:
    def __init__(self, offsets, occ, idxs):
        self._offsets = offsets+1
        self._occ = occ
        self._idxs = idxs

    def backward(self, s, e, c):
        print(s, e, c)
        idx = lookup[c]
        new_s = self._occ[s, idx]
        new_e = self._occ[e, idx]
        return new_s + self._offsets[idx], new_e + self._offsets[idx]

    def find_seq(self, seq):
        seq = seq[::-1]
        idx = lookup[seq[0]]
        s, e = self._offsets[idx:idx+2]
        for c in seq[1:]:
            s, e = self.backward(s, e, c)
        return self._idxs[s:e]
