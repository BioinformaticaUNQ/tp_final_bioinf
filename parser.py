from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline 

file = open("./fasta/6TDB.fasta", "w")

blast_records = NCBIXML.parse(open("6TDB.xml"))
for blast_record in blast_records:
        for alignment in blast_record.alignments:
            file.writelines(">" + alignment.hit_id)
            file.writelines("\n")
            for hsp in alignment.hsps:
                sbjct_no_gaps = hsp.sbjct.replace("-","")
                file.writelines(sbjct_no_gaps)
                file.writelines("\n")

clustalomega_cline = ClustalOmegaCommandline(infile = "./fasta/6TDB.fasta", outfile = "./fasta/6TDB_Alligned.fasta", outfmt = 'phylip', verbose = True, auto = False)
clustalomega_cline()
# print the executable command
print(clustalomega_cline)
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