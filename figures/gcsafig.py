from alignment import AlignmentGraph
from seqgraph import SeqGraph
from pylatex import Figure

if __name__ == "__main__":
    graph1 = AlignmentGraph("GCACG", "GCTCG").generate_tikz("seq")
    edges = {0: [1], 1: [2, 3], 2: [4], 3: [4], 4: [5]}
    seqgraph1 = SeqGraph("GCATCG", edges)
    figure = Figure("gcsafiga", "gcsafig")
    figure.write(graph1)
