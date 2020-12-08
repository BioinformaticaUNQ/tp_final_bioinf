from Bio.Align.Applications import ClustalOmegaCommandline
import logging 

def run_clustal(pdb_id,fasta_seq,outputPath):
    outputPath = outputPath + "/"+pdb_id+"_aln.fasta"
    clustalomega_cline = ClustalOmegaCommandline(infile = fasta_seq, outfile = outputPath,force = True)
    clustalomega_cline()
    logging.info("Se ejecuto clustal omega con los siguientes parametros")
    logging.info("infile = " + fasta_seq + ", outfile = " + outputPath + ",force = True")
    return outputPath