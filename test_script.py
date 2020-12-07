import sys
import os

seq_map = {}
with open("./fasta/6TDB_aln.fasta", "r") as f:
        seqText = ""
        key = ""

        for sq in f:
            if '>' in sq:
                if seqText != "":
                    seq_map[key] = seqText
                    seqText = "" 
                if '6TDB' not in sq:
                    seq_map[sq.split("|")[1]] = ''
                    key = sq.split("|")[1]
                else:
                    seq_map['6TDB'] = ''
                    key = '6TDB'

            else:
                seqText += sq.replace('\n',"")
        seq_map[key] = seqText
print(seq_map)
