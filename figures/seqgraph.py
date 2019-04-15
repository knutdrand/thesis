import numpy as np
from pylatex import body, tikz_matrix, Figure
GAP = -1


def score_matrix_g(graph_a, graph_b, is_global=True):
    seq_a = graph_a._labels
    seq_b = graph_b._labels
    match = np.array([[c_a == c_b for c_b in seq_b] for c_a in seq_a])
    scores = np.where(match, 0, -1)
    D = np.zeros((len(seq_a), len(seq_b)), dtype="int")
    D[0, 1:] = GAP*graph_b._dists[:-1]
    if is_global:
        D[1:, 0] = GAP*graph_a.dists[:-1]
    for i in range(1, len(seq_a)):
        for j in range(1, len(seq_b)):
            d_a = max(D[graph_a._rev_edges[i], j])+GAP
            d_b = max(D[i, graph_b._rev_edges[j]])+GAP
            d_m = max(D[k, l] + scores[i, j]
                      for k in graph_a._rev_edges[i]
                      for l in graph_b._rev_edges[j])
            D[i, j] = max(d_a, d_b, d_m)
    return D


def backtrack_g(D, graph, graph2, is_global=True):
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
        if not graph2._rev_edges[j] and not graph._rev_edges[i]:
            break

        found = False
        for l in graph2._rev_edges[j]:
            if D[i, l]+GAP == D[idx]:
                idx = (i, l)
                found = True
                break
        if found:
            continue
        for k in graph._rev_edges[i]:
            if D[i-1, j]+GAP == D[idx]:
                idx = (k, j)
                found = True
                break
        m = 0 if graph._labels[i] == graph2._labels[j] else -1
        for k in graph._rev_edges[i]:
            if found:
                break
            for l in graph2._rev_edges[j]:
                if D[k, l] + m == D[idx]:
                    idx = (k, l)
                    found = True
                    break
    return path


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
    a.append(seq_a[i])
    b.append(seq_b[j])
    return "".join(a), "".join(b)


class SeqGraph:
    def __init__(self, labels, edges):
        self._labels = labels
        print(self._labels)
        self._edges = edges
        print(self._edges)
        self._rev_edges = self.get_preds()
        self.dists = self.get_dists()
        self._dists = self.dists

    def get_dists(self):
        dists = np.zeros(len(self._labels))
        for i in range(1, len(self._labels)):
            prevs = self._rev_edges[i]
            dists[i] = min(dists[prevs])+1
        return dists

    def get_preds(self):
        return {node: [i for i, js in self._edges.iteritems() if node in js]
                for node in range(len(self._labels))}

    def get_edge(self, i, j, horizontal=True):
        d = "left" if horizontal else "right"
        coord = "m-1-%s" if horizontal else "m-%s-1"
        if j == i+1:
            tmplt = r"\path [->] (%s) edge node {} (%s);" % (coord, coord)
            # tmplt = r"\draw (m-1-%s) -- (m-1-%s);"
        else:
            tmplt = r"\path [->] (%s) edge[bend %s=60] node {} (%s);" % (coord, d, coord)
        return tmplt % (i+2, j+2)

    def to_tikz(self, horizontal=True):
        lines = []
        labels = [self._labels] if horizontal else [[c] for c in self._labels]
        lines.append(tikz_matrix("m", body(labels)))
        lines.extend(self.get_edge(i, j, horizontal) for i, js in self._edges.iteritems() for j in js)
        return "\n".join(lines)

    def get_tikz_edges(self, horizontal=True):
        lines = [self.get_edge(i, j, horizontal) for i, js in self._edges.iteritems() for j in js]
        return "\n".join(lines)

if __name__ == "__main__":
    labels = "#ACTGGG"
    edges = dict(enumerate([[1], [2, 3], [4], [4], [5, 6], [6]]))
    seq_graph2 = SeqGraph(labels, edges)
    labels2 = "#ATAG"
    edges2 = dict(zip(range(4), ([i] for i in range(1, 5))))
    edges2[1].append(3)
    seq_graph = SeqGraph(labels2, edges2)
    m = score_matrix_g(seq_graph, seq_graph2)
    path = backtrack_g(m, seq_graph, seq_graph2)
    print path
    figure = Figure("seqgraph_a", "seqalign_figs")
    c_text = r"\textcolor{red}{%s}"
    b_text = r"\textcolor{black!25}{%s}"
    m = [[c_text % c if (i, j) in path else b_text % c for j, c in enumerate(row)]
         for i, row in enumerate(m)]
    elems = [list("."+labels)] + [[c] + row for c, row in zip(labels2, m)]
    edges2 = seq_graph2.get_tikz_edges()
    edges = seq_graph.get_tikz_edges(False)
    path_edge_template = r"\path[-stealth] (m-%s-%s) edge[color=red] (m-%s-%s);"
    path_edges = "\n".join([path_edge_template % (i+2, j+2, k+2, l+2) for (i, j), (k, l) in zip(path[:-1], path[1:])])

    A, B = get_alignment(seq_graph._labels, seq_graph2._labels, path)
    alignment_matrix = tikz_matrix("am", body([A, B]), "below=of m")
    content = "\n".join([tikz_matrix("m", body(elems)), edges, edges2, path_edges, alignment_matrix])
    content = content.replace("#", r"\#")
    figure.write(content)
    # figure.write(seq_graph.to_tikz(False))
    # print(edges)
