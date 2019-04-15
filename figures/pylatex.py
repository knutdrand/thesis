import subprocess


def body(rows):
    return "\\\\ \n".join(" & ".join(row) for row in rows)+r"\\"


def tikz_matrix(name, body, args=""):
    return r"\matrix (%s) [matrix of math nodes, align=right, row sep=0.1em,column sep=0.1em, minimum width=1em, minimum height=1em, nodes={align=right}, %s]{%s};" % (name, args, body)

class Figure:
    def __init__(self, name, master):
        self.file_name = name + ".tex"
        self.master = master + ".tex"

    def write(self, content):
        open(self.file_name, "w").write(content)
        subprocess.call(["pdflatex", self.master])
