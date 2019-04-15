from Bio.pairwise2 import align
import numpy as np
import subprocess
from pylatex import tikz_matrix, body


class AlignmentGraph:
    header = r"""\documentclass[convert={density=300,outext=.png}]{standalone}
\usepackage{tikz}
\usetikzlibrary{matrix}
\begin{document}
\begin{tikzpicture}[shorten >=1pt,->]
  \def\f{1.2};
  \def\fy{1}

\tikzstyle{vertex}=[circle,fill=black!25,minimum size=17*\f pt,inner sep=0pt]"""
    footer = """\end{tikzpicture}
\end{document}
"""

    def __init__(self, row_a, row_b, x=0, y=0, node_id=0):
        self._row_a = row_a
        self._row_b = row_b
        self._x = x
        self._y = y
        self.lines = []
        self._node_idx = node_id
        self._matrix = np.array([[""]*len(row_a)]*2)

    def node(self, c, x, y):
        line = r"""\node[vertex][] (P%s) at (%s*\f,%s) {%s};""" % (self._node_idx+1, self._x+x, -self._y + y, c)
        self._node_idx += 1
        # self.lines.append(line)
        self._matrix[y, x] = c
        return "%s-%s-%s" % (self.name, y+1, x+1) #self._node_idx

    def edge(self, i, j):
        line = r"\draw (%s) -- (%s);" % (i, j)
        self.lines.append(line)

    def del_edge(self, i, j):
        line = r"\path [->] (%s) edge[bend right=60] node {} (%s);" % (i, j)
        self.lines.append(line)

    def generate_tikz(self, name, args=""):
        self.name = name
        prev_a_id = None
        prev_b_id = None
        del_a = False
        del_b = False
        for x, (c_a, c_b) in enumerate(zip(self._row_a, self._row_b)):
            y = 0
            if c_a != "-":
                a_id = self.node(c_a, x, y)
                y += 1
                if prev_a_id is not None:
                    if del_a:
                        self.del_edge(prev_a_id, a_id)
                    else:
                        self.edge(prev_a_id, a_id)
                prev_a_id = a_id
                del_a = False
            else:
                del_a = True
            if c_b != "-":
                if c_a != c_b:
                    b_id = self.node(c_b, x, y)
                else:
                    b_id = a_id
                if prev_b_id is not None:
                    if del_b:
                        self.del_edge(prev_b_id, b_id)
                    else:
                        self.edge(prev_b_id, b_id)
                prev_b_id = b_id
                del_b = False
            else:
                del_b = True
        return tikz_matrix(name, body(self._matrix), args) + "\n" + "\n".join(self.lines)


    def to_file(self, filename):
        f = open(filename, "w")
        f.write(self.header + "\n")
        # f.write(self._row_a + "\n")
        f.write(tikz_matrix("seq", body(self._matrix))+"\n")
        for line in self.lines:
            f.write(line+"\n")
        f.write(self.footer)
        f.close()
        subprocess.call(["pdflatex", filename])

if __name__ == "__main__":
    seq_a = "ACGGG"
    seq_b = "ATGG"
    # seq_b = "ACCGGTAAGGT"
    # seq_a = "ACGGTGTACAT"
    # seq_b = "ACCGGTAAGGT"

    row_a, row_b, _, _, _ =  align.globalms(seq_a, seq_b, 0, -1, -1, -1)[0]
    print(row_a)
    print(row_b)
    graph = AlignmentGraph(row_a, row_b)
    graph.generate_tikz()
    graph.to_file("testing.tex")

#for line in graph.lines:
#
