from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from tkinter import messagebox
import logging 

def run_clustal(pdb_id,fasta_seq,outputPath):
    outputPath = outputPath + "/"+pdb_id+"_aln.fasta"
    seqs = list(SeqIO.parse(fasta_seq,"fasta"))
    if(len(seqs) <= 1):
        messagebox.showerror("INFO", "No se encontraron secuencias para alinear")
        return
    clustalomega_cline = ClustalOmegaCommandline(infile = fasta_seq, outfile = outputPath,force = True)
    clustalomega_cline()
    logging.info("Se ejecuto clustal omega con los siguientes parametros")
    logging.info("infile = " + fasta_seq + ", outfile = " + outputPath + ",force = True")
    return outputPath