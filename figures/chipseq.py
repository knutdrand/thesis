def reverse_compliment(seq):
    lookup = {"A": "T", "T": "A",
              "G": "C", "C": "G"}
    return "".join(lookup[c] for c in seq[::-1])

seqs = ["TCGAGCTGTA", "CTACTCGAGC", "TTGATAGATA", "TGTAGTTTGA"]
for seq in seqs:
    print reverse_compliment(seq)


TACAGCTCGA
GCTCGAGTAG
TATCTATCAA
TCAAACTACA

