import numpy as np
import subprocess
from collections import namedtuple


lookup = {"A": 0, "C": 1, "G": 2, "T": 3}
rev_c = {"A": "T", "C": "G", "G": "C", "T": "A", "$": "$"}

BiInt = namedtuple("BiInt", ["k", "l", "s"])


class FMIndex:
    def __init__(self, offsets, occ, idxs):
        self._offsets = np.asanyarray(offsets+1, dtype="int")
        self._occ = np.asanyarray(occ, dtype="int")
        self._idxs = idxs

    def backward(self, s, e, c):
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

    def backward_ext(self, bint, char):
        rev_char = rev_c[char]
        idx = lookup[char]
        ks = self._offsets[:-1] + self._occ[bint.k, :]
        ss = self._occ[bint.k + bint.s, :]-self._occ[bint.k, :]
        rev_idx = lookup[rev_char]
        new_l = bint.l+bint.s - np.sum(ss[:4-rev_idx])
        return BiInt(ks[idx], new_l, ss[idx])

    def forward_ext(self, bint, char):
        rev_bint = BiInt(bint.l, bint.k, bint.s)
        new_rev_bint = self.backward_ext(rev_bint, rev_c[char])
        return BiInt(new_rev_bint.l, new_rev_bint.k, new_rev_bint.s)

    def get_bint(self, c):
        idx = lookup[c]
        rev_idx = lookup[rev_c[c]]
        k = self._offsets[idx]
        s = self._offsets[idx+1]-k
        l = self._offsets[rev_idx]
        return BiInt(k, l, s)

    def find_smem(self, P, i):
        bint = self.get_bint(P[i])
        longest_bints = []
        for j in range(1, i+1):
            new_bint = self.backward_ext(bint, P[i-j])
            if new_bint.s != bint.s:
                longest_bints.append(bint)
            bint = new_bint
            if bint.s == 0:
                break
        else:
            longest_bints.append(bint)
        smems = []
        prev = longest_bints[::-1]
        for j in range(i+1, len(P)):
            curr = []
            longer_found = False
            for bint in prev:
                new_bint = self.forward_ext(bint, P[j])
                if new_bint.s == 0:
                    if not longer_found:
                        smems.append(bint)
                        longer_found = True
                else:
                    curr.append(new_bint)
                    longer_found = True
            prev = curr
            if not prev:
                break
        else:
            smems.append(prev[0])
        print(smems)
        return smems


def sort_suffices(seq):
    suffices = [seq[i:]+seq[:i] for i in range(len(seq))]
    t = sorted((suffix, i) for i, suffix in enumerate(suffices))
    suffices, idxs = zip(*t)
    return idxs, suffices


def get_occurance_matrix(bwt):
    occ = np.zeros((len(bwt)+1, 4))
    for i, c in enumerate(bwt):
        if c == "$":
            continue
        occ[i+1, lookup[c]] += 1
    occ = np.cumsum(occ, axis=0)
    return occ


def get_char_counts(sorted_seq):
    sorted_seq = np.array(sorted_seq)
    t = np.flatnonzero(sorted_seq[1:] != sorted_seq[:-1])
    return t


def format_table(array):
    return "\n".join(" && ".join(str(c) for c in row) for row in array)


def table1(idxs, suffices, occ, counts):
    occ = np.asanyarray(occ, dtype="int")
    header = r"""\documentclass{article}
\usepackage{colortbl}
\begin{document}
\begin{table}
\begin{tabular}{|r|>{\columncolor[gray]{0.8}}r|%s|>{\columncolor[gray]{0.8}}r| r |rrrr|}
\hline
SA & &  %s BWT & OCC & A & C & G & T \\
\hline""" % ("r"*(len(suffices[0])-2), "&   "*(len(suffices[0])-2))
    rows = [str(i)+ " & " + " & ".join(str(c) for c in row) for i, row in zip(idxs, suffices)]
    occ_rows = ["&&" + "&".join(str(n) for n in row)+"\\\\" for row in occ]
    rows = [row + occ_row for row, occ_row in zip(rows, occ_rows)]
    rows = [row+"\n \\hline" if i in counts else row for i, row in enumerate(rows)]
    body = "\n".join(rows)
    body = body.replace("$", r"\$")
    footer = r"\end{tabular}\end{table}\end{document}"
    return header + "\n" + body + "\n" + footer

if __name__ == "__main__":
    # seq = "ATGTCATTGGA"
    seq = "ATTGAC"
    seq = seq+"$" + "".join(rev_c[c] for c in seq[::-1]) + "$"
    idxs, suffices = sort_suffices(seq)
    # print(format_table(suffices))

    bwt = [s[-1] for s in suffices]
    occ = get_occurance_matrix(bwt)
    counts = get_char_counts([s[0] for s in suffices] + ["$"])
    table = table1(idxs, suffices, occ, counts)
    # open("fm.tex", "w").write(table)
    # subprocess.call(["pdflatex", "fm.tex"])
    # print(counts)
    fm = FMIndex(counts, occ, np.array(idxs))
    for i in range(6):
        fm.find_smem("ATTCAC", i)
    # 
    # print fm.find_seq("ACGT")
    # print fm.find_seq("CG")
    # for s in suffices:
    #     print s
