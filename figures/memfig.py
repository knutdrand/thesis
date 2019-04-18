from pylatex import Figure, tikz_matrix, body

rev_c = {"A": "T", "C": "G", "G": "C", "T": "A", "$": "$"}
colors = ["red", "green", "blue", "yellow"]


def create_figure(seq, interval_dict, query):
    seq_fig = tikz_matrix("seq", body(
        [seq, "".join(rev_c[c] for c in seq)]))
    query_fig = tikz_matrix("q", body([query]), "below=of seq")
    rec_1 = r"\node[color=%s, draw, dashed, rounded corners, fit=(seq-1-%s) (seq-1-%s)] {};"
    rec_2 = r"\node[color=%s, draw, dashed, rounded corners, fit=(seq-2-%s) (seq-2-%s)] {};"
    rec_q = r"\node[color=%s, draw, dashed, rounded corners, fit=(q-1-%s) (q-1-%s)] {};"
    lines = []
    N = len(seq)
    i = 0
    for (s, e), (intervals, reverse_intervals) in interval_dict.iteritems():
        lines.append(rec_q % (colors[i], s+1, e))
        for (start, end) in intervals:
            lines.append(rec_1 % (colors[i], start+1, end))

        for start, end in reverse_intervals:
            lines.append(rec_2 % (colors[i], N-end+1, N-start))
        i += 1
    content = seq_fig + "\n" + query_fig + "\n" + "\n".join(lines)
    content = content.replace("$", r"\$")
    return content


if __name__ == "__main__":
    fig = Figure("memfig", "memfigs")
    fig.write(create_figure("ATTCT", {1: ([(0, 3), (2, 5)], [(1, 3)])}))
