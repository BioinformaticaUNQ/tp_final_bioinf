from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline 
import logomaker as lm
import matplotlib.pyplot as plt



# file = open("./fasta/6TDB.fasta", "w")

# blast_records = NCBIXML.parse(open("6TDB.xml"))
# original = open("6TDB.fasta","r")
# for line in original:
#     file.writelines(line)
# for blast_record in blast_records:
#         for alignment in blast_record.alignments:
#             file.writelines(">" + alignment.hit_id)
#             file.writelines("\n")
#             for hsp in alignment.hsps:
#                 sbjct_no_gaps = hsp.sbjct.replace("-","")
#                 file.writelines(sbjct_no_gaps)
#                 file.writelines("\n")
outputPath = "./fasta/6TDB_Alligned.fasta"
clustalomega_cline = ClustalOmegaCommandline(infile = "./fasta/6TDB.fasta", outfile = outputPath,force = True)
clustalomega_cline()
raw_seqs =[]
with open(outputPath, "r") as f:
    seqText = ""
    for sq in f:
        if '>' in sq:
            if seqText != "":
                raw_seqs.append(seqText)
                seqText = "" 
            raw_seqs.append(sq)
        else:
            seqText += sq
seqs = [seq.strip() for seq in raw_seqs if ('#' not in seq) and ('>') not in seq]

counts_mat = lm.alignment_to_matrix(seqs)
counts_mat.head()
lm.Logo(counts_mat)
plt.savefig("test.png")
plt.show()

# print the executable command
#print(clustalomega_cline)
#             file.writelines(hit.id)
#             file.writelines(",")
# file.close()
# ids = dict([(x, True) for x in open('ids.txt','r').read().split(',')])
# print(ids)
# for record in SeqIO.parse(open('./fasta/6TDB.fasta','r'),'fasta'):
#     print(record)
#     if record.id in ids:
#         print(">" + record.id + "\n" + str(record.seq))
#             #for hsp in hit:
#                 #print(hsp)
#                 #file.writelines(">" +hit.description)
#                 #file.writelines("\n")
#                 #file.writelines(hsp.hit_string)
#                 #file.writelines("\n")