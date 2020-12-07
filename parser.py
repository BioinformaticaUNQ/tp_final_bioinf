from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline 
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file


# p = PDBParser()
# structure = p.get_structure("6TDB", "./6TDB.pdb")
# model = structure[0]
# dssp = DSSP(model, "./6TDB.pdb")
# print(len(list(dssp.keys())))
# # DSSP data is accessed by a tuple (chain_id, res_id)
# for a_key in list(dssp.keys()):
 
# (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
# NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
# NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)
    # if dssp[a_key][2] != '-':
    #     print(dssp[a_key])
    #print(dssp[a_key][2])

# file = open("./fasta/6TDB.fasta", "w")
dssp_tuple = dssp_dict_from_pdb_file("6TDB.pdb")
dssp_dict = dssp_tuple[0]
#EL PRIMER VALOR DE LA TUPLA ES UN DICCIONARIO (TUPLA KEY, DATA DE LA ESTRUCTURA)
#EL SEGUNDO VALOR ES LA LISTA DE KEYS QUE SON DEL FORMATO ("CADENA",('',NRO DE RESIDUO,''))
#LAS CADENAS ESTAN SEPARADAS , POR EJ LA CADENA A SON MUCHAS KEYS TODAS EMPEZANDO CON A PERO CON DISTINTO NUMERO DE RESIDUO
AChain = []
for key in dssp_tuple[1]:
    if(key[0] == 'A'):
        #OBTENGO TODAS LAS KEYS DE LA CADENA A
        AChain.append(key)
for chainPart in AChain:
    #OBTENGO LAS ESTRUCTURAS DE LA CADENA A
    print(dssp_dict[chainPart][0])
# print(dssp_dict[('A', (' ', 272, ' '))])
# print(dssp_dict[('A', (' ', 273, ' '))])
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
# outputPath = "./fasta/6TDB_Alligned.fasta"
# clustalomega_cline = ClustalOmegaCommandline(infile = "./fasta/6TDB.fasta", outfile = outputPath,force = True)
# clustalomega_cline()
# raw_seqs =[]
# with open(outputPath, "r") as f:
#     seqText = ""
#     for sq in f:
#         if '>' in sq:
#             if seqText != "":
#                 raw_seqs.append(seqText)
#                 seqText = "" 
#             raw_seqs.append(sq)
#         else:
#             seqText += sq
# seqs = [seq.strip() for seq in raw_seqs if ('#' not in seq) and ('>') not in seq]

# counts_mat = lm.alignment_to_matrix(seqs)
# counts_mat.head()
# lm.Logo(counts_mat)
# plt.savefig("test.png")
# plt.show()

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