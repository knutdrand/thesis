from Bio.pairwise2 import align
import alignment


def get_graph(seq1, seq2, x=0, y=0, node_id=0):
    row_a, row_b, _, _, _ =  align.globalms(seq1, seq2, 0, -1, -1, -0.5)[0]
    print(row_a)
    print(row_b)
    return alignment.AlignmentGraph(row_a, row_b, x, y, node_id)


ref_seq =    "ACGTTACGGTACCGTA"
ref_seq2 =   "ACCGTTGCGGGTACCGTA"
get_graph(ref_seq, ref_seq2).to_file("ref_graph.tex")


query_seq = "TTGCGGAC"
smems2 = ["TTGCGG", "AC"]

smems = [
    "TT",
    "TT",
    "G",
    "CGG",
    "CGG",
    "CGG",
    "CGG",
    "AC",
    "CC"]

smems = [
    ("TT", 0, 3),
    ("G", 2, [2, 7, 8, 13]),
    ("CGG", 3, 6),
    ("AC", 6, [0, 5, 10, 15])]

ref_row, query_row = align.globalms(ref_seq, query_seq, 0, -1, -1, -1, penalize_end_gaps=(True, False))[0][:2]
print(ref_row)
print(query_row)
