from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIWWW, NCBIXML

import os

def blastp_query(pdb_id,evalue,coverage,data,sequence):
    blast = pdb_id + ".xml"
    cline = NcbiblastpCommandline(query=pdb_id + '.fasta', db="./db/pdbaa",
                              evalue=evalue, out=blast, outfmt=5,qcov_hsp_perc=coverage)
    cline()
    blast_records = NCBIXML.parse(open(blast))
    if not os.path.exists("./fasta"):
        os.mkdir("./fasta")
    
    all_seq_fasta = "./fasta/"+pdb_id+".fasta"
    file = open(all_seq_fasta, "w")
    file.writelines(data + '\n')
    file.writelines(sequence + '\n')
    identity_perc = 40
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            identity = alignment.hsps[0].identities
            align_length = alignment.hsps[0].align_length
            percentage =  identity / align_length * 100
            if (percentage >= identity_perc):
                file.writelines(">" + alignment.hit_id)
                file.writelines("\n")
                for hsp in alignment.hsps:
                    sbjct_no_gaps = hsp.sbjct.replace("-","")
                    file.writelines(sbjct_no_gaps)
                    file.writelines("\n")
    
    file.close()
    return all_seq_fasta