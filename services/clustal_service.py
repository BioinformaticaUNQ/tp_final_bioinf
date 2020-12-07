from Bio.Align.Applications import ClustalOmegaCommandline 

def run_clustal(pdb_id,fasta_seq,outputPath):
    print(outputPath)
    outputPath = outputPath + "/"+pdb_id+"_aln.fasta"
    clustalomega_cline = ClustalOmegaCommandline(infile = fasta_seq, outfile = outputPath,force = True)
    clustalomega_cline()
    return outputPath