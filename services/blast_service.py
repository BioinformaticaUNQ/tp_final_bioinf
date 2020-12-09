from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIWWW, NCBIXML

import os

def blastp_query(pdb_id,evalue,coverage,data,sequence,input_path):
    blast = pdb_id + ".xml"
    cline = NcbiblastpCommandline(query=input_path + "/" + pdb_id + '.fasta', db="./db/pdbaa",
                              evalue=evalue, out=blast, outfmt=5,qcov_hsp_perc=coverage)
    cline()
    blast_records = NCBIXML.parse(open(blast))
    if not os.path.exists("./fasta"):
        os.mkdir("./fasta")
    
    all_seq_fasta = input_path + "/" + pdb_id + "_blastp.fasta"

    pdbs_to_process = []
    identity_perc = 40
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            identity = alignment.hsps[0].identities
            align_length = alignment.hsps[0].align_length
            percentage =  identity / align_length * 100
            if (percentage >= identity_perc):
                pdbs_to_process.append(alignment.hit_id.split('|')[1])
                # file.writelines(">" + alignment.hit_id)
                # file.writelines("\n")
                
                # for hsp in alignment.hsps:
                #     sbjct_no_gaps = hsp.sbjct.replace("-","")
                #     file.writelines(sbjct_no_gaps)
                #     file.writelines("\n")
    get_complete_seqs = "blastdbcmd -db ./db/pdbaa -entry {} -out {}".format(
            ",".join(pdbs_to_process),
            all_seq_fasta
        )
    
    os.system(get_complete_seqs)
    file = open(all_seq_fasta, "a")

    with open(input_path + "/" + pdb_id + ".fasta", "r") as f:    
        for line in f:
            file.writelines(line)
    file.close()
    return all_seq_fasta