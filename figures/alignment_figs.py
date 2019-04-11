import alignment
from Bio.pairwise2 import align


def get_graph(seq1, seq2, x=0, y=0, node_id=0):
    row_a, row_b, _, _, _ =  align.globalms(seq1, seq2, 0, -1, -1, -0.5)[0]
    print(row_a)
    print(row_b)
    return alignment.AlignmentGraph(row_a, row_b, x, y, node_id)


orig_seq = "ACGTACGTACGT"
seq_a = "ACGGACGTCCGT"
seq_aa = "ACGGATCGTCCCT"
seq_ab = "CGGATCTTCCCT"
seq_b =  "ACGTACCGTAGGT"
seq_ba = "ATGTACCGTACGT"
seq_bb =  "ACGTACCGTTTAGGT"

lines = []
x_dist = 15
y_dist = 1
node_id = 0
for i, seq in enumerate([seq_aa, seq_ab, seq_ba, seq_bb]):
    x = (i // 2)*x_dist
    y = (i % 2)*y_dist
    print(x, y)
    graph = alignment.AlignmentGraph(seq, seq, x=x, y=y, node_id=node_id)
    node_id = graph._node_idx
    lines.extend(graph.lines)

graph_a = get_graph(seq_aa, seq_ab, 0, 4, node_id)
node_id = graph_a._node_idx
lines.extend(graph_a.lines)
graph_b = get_graph(seq_ba, seq_bb, x_dist, 4, node_id)
lines.extend(graph_b.lines)
node_id = graph_b._node_idx
graph = get_graph(seq_a, seq_b, x_dist//2, 8, node_id)
lines.extend(graph.lines)
graph.lines = lines
graph.to_file("graph_msa.tex")
for line in lines:
    print(line)


# graph_b.to_file("b_graph.tex")
# get_graph(seq_aa, seq_ab).to_file("a_graph.tex")

