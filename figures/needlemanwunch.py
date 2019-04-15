import numpy as np
from pylatex import Figure, body


GAP = -1
HEADER = r"""\documentclass{article}
\usepackage{tikz}
\usetikzlibrary{matrix}
\begin{document}
\begin{tikzpicture}
  \matrix (m) [matrix of math nodes, align=right, row sep=0.1em,column sep=0.1em, minimum width=1em, minimum height=1em, nodes={align=right}]
"""
FOOTER = r"""
\end{tikzpicture}
\end{document}
"""

matrix_header = r"\matrix (m) [matrix of math nodes, align=right, row sep=0.1em,column sep=0.1em, minimum width=1em, minimum height=1em, nodes={align=right}]"


def tikz_matrix(name, body, args=""):
    return r"\matrix (%s) [matrix of math nodes, align=right, row sep=0.1em,column sep=0.1em, minimum width=1em, minimum height=1em, nodes={align=right}, %s]{%s};" % (name, args, body)


def score_matrix(seq_a, seq_b, is_global=True):
    match = np.array([[c_a == c_b for c_b in seq_b] for c_a in seq_a])
    scores = np.where(match, 0, -1)
    D = np.zeros((len(seq_a)+1, len(seq_b)+1), dtype="int")
    D[0, :] = GAP*np.arange(len(seq_b)+1)
    if is_global:
        D[:, 0] = GAP*np.arange(len(seq_a)+1)
    for i in range(len(seq_a)):
        for j in range(len(seq_b)):
            D[i+1, j+1] = max(D[i, j+1]+GAP,
                              D[i+1, j]+GAP,
                              D[i, j]+scores[i, j])
    return D


def score_matrix_g(graph_a, graph_b, is_global=True):
    seq_a = graph_a._labels
    seq_b = graph_b._labels
    match = np.array([[c_a == c_b for c_b in seq_b] for c_a in seq_a])
    scores = np.where(match, 0, -1)
    D = np.zeros((len(seq_a)+1, len(seq_b)+1), dtype="int")
    D[0, 1:] = GAP*graph_b._dists
    if is_global:
        D[1:, 0] = GAP*graph_a.dists
    for i in range(1, len(seq_a)+1):
        for j in range(1, len(seq_b)+1):
            d_a = max(D[graph_a._rev_edges[i], j])+GAP
            d_b = max(D[i, graph_b._rev_edges[j]])+GAP
            d_m = max(D[k, l] + scores[k, l]
                      for k in graph_a._rev_edges[i]
                      for l in graph_b._rev_edges[j])
            return max(d_a, d_b, d_m)
    return D


def backtrack(D, is_global=True):
    path = []
    i, j = D.shape
    i -= 1
    j -= 1
    if not is_global:
        i = np.argmax(D[:, -1])
    idx = (i, j)
    while True:
        path.append(idx)
        i, j = idx
        if j > 0 and D[i, j-1]+GAP == D[idx]:
            idx = (i, j-1)
            continue
        elif i > 0 and D[i-1, j]+GAP == D[idx]:
            idx = (i-1, j)
            continue
        if j == 0 and (not is_global or i == 0):
            break
        idx = (i-1, j-1)
    return path


def format_matrix(D, path, seq_a, seq_b):
    path = set(path)
    lines = []
    color = r"\cellcolor{blue!25}"
    for i, row in enumerate(D):
        line = []
        for j, d in enumerate(row):
            if (i, j) in path:
                line.append(color + str(d))
            else:
                line.append(str(d))
        lines.append("& " + " && ".join(line) + r"\\")
        if i != D.shape[0]-1:
            lines.append(seq_a[i]+r"\\")
    first_row = "&& " + " && ".join(seq_b) + " &"
    print(first_row)
    header = "c"*(D.shape[1]+D.shape[1])
    print(header)
    return "\n".join(lines)


def format_matrix_tikz(D, path, seq_a, seq_b):
    path_set = set(path)
    lines = []
    path_color = "red!25"
    c_text = r"\textcolor{red}{%s}"
    b_text = r"\textcolor{black!25}{%s}"
    color = r"\cellcolor{blue!25}"
    elements = [[c_text % c if (i, j) in path_set else b_text % c for j, c in enumerate(row)] for i, row in enumerate(D)]
    m = tikz_matrix("m", body(elements))
    seq_template = r"\path[-stealth] (m-%s-1) edge [draw=none] node [left=0.5em] {%s} (m-%s-1);"
    seq_template_b = r"\path[-stealth] (m-%s-%s) edge [draw=none] node [above=0.5em] {%s} (m-%s-%s);"
    seq_a_text = "\n".join(seq_template %
                           (i+1, c, i+2) for i, c in enumerate(seq_a))
    seq_b_text = "\n".join(seq_template_b % (1, i+1, c, 1, i+2) for i, c in enumerate(seq_b))

    edge_template = r"\path[-stealth] (m-%s-%s) edge[color=red] (m-%s-%s);"
    edges = "\n".join([edge_template % (i+1, j+1, k+1, l+1) for (i, j), (k, l) in zip(path[:-1], path[1:])])
    return "\n".join((m, seq_a_text, seq_b_text, edges))


def get_alignment(seq_a, seq_b, path):
    path = path[::-1]
    a, b = ([], [])
    prev_i, prev_j = path[0]
    for i, j in path[1:]:
        if i > prev_i:
            a.append(seq_a[prev_i])
        else:
            a.append("-")
        if j > prev_j:
            b.append(seq_b[prev_j])
        else:
            b.append("-")
        prev_i, prev_j = (i, j)
    return "".join(a), "".join(b)


def main(seq_a, seq_b, is_global, figure):
    D = score_matrix(seq_a, seq_b, is_global=False)
    path = backtrack(D, is_global)
    f = format_matrix_tikz(D, path, seq_a, seq_b)
    A, B = get_alignment(seq_a, seq_b, path)
    block = tikz_matrix("block", body([A, B]), "below=of m")
    start = tikz_matrix("seqs", body([seq_a, seq_b]), "above=of m")
    figure.write(f + "\n" + block + "\n" + start)
    # open("needlemanwunch.tex", "w").write(HEADER + f + "\n" + block + FOOTER)
    " TODO: Add hlines and vlines to denote which letters are in the block "

if __name__ == "__main__":
    seq_a = "ACGTTAC"
    seq_b = "ACTTGA"
    figure = Figure("needleman_a", "needleman_figures")
    # main(seq_a, seq_b, True, figure)
    ref = "GTGG" + seq_a + "GCGGT"
    figure2 = Figure("needleman_b", "needleman_figures")
    main(ref, seq_b, False, figure2)
