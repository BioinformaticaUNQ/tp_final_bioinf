from tkinter import *
from Bio.Blast import NCBIWWW, NCBIXML
import urllib.request

import os
os.system('clear')

root = Tk()
root.title("TP FINAL")
root.geometry("800x600")

def fasta():
    print('Beginning file download with urllib2...')
    pdb_id = str(myTextbox.get())
    url = "https://files.rcsb.org/download/"+ pdb_id +".pdb"
    urllib.request.urlretrieve(url, pdb_id + ".pdb")
    url2 = "https://www.rcsb.org/fasta/entry/"+pdb_id+"/download"
    urllib.request.urlretrieve(url2, pdb_id + ".fasta")
    fasta_string = open(pdb_id + ".fasta").read()
    print(fasta_string)
    result_handle = NCBIWWW.qblast("blastp", "pdb", fasta_string.split('\n')[1],
        short_query=True,
        hitlist_size=50,
        expect=1000,
        word_size=7,
        nucl_reward=1,
        matrix_name='BLOSUM62',
        nucl_penalty=-3,
        threshold=0.05,
        gapcosts="11 1")
    blast = pdb_id + ".xml"

    save_clk = open(blast, "w")
    save_clk.write(result_handle.read())    
    save_clk.close()
    
    blast_records = NCBIXML.parse(open(blast))
    print(blast_records)
    results = []
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            aligned_sequence = ''
            aligned_matches = ''
            for hsp in alignment.hsps:
                aligned_sequence = aligned_sequence + hsp.sbjct
                aligned_matches = aligned_matches + hsp.match

            print(alignment.title.split('>')[0])
            print(aligmened_sequence.replace('-', ''))
            print(aligmened_matches)
            print(aligmened_sequence)
            


myLabel = Label(root, text="Ingresar codigo PDB")
myLabel.pack()

myTextbox = Entry(root, width=50)
myTextbox.pack()

myButton = Button(root, text="Procesar", command = fasta)
myButton.pack()

root.mainloop()