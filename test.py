from tkinter import *
from tkinter import messagebox
from Bio.Blast import NCBIWWW, NCBIXML
import urllib.request
import pymol
import os

root = Tk()
root.title("TP FINAL")
root.geometry("800x600")

def fasta():
    print('Beginning file download with urllib2...')
    pdb_id = str(myTextbox.get())
    url = "https://files.rcsb.org/download/"+ pdb_id +".pdb"
    try:
        urllib.request.urlretrieve(url, pdb_id + ".pdb")
    except Exception as e:
        print(str(e))
        if(str(e) == "HTTP Error 404: Not Found"):
            messagebox.showerror("Error", "Código pdb inválido")
        else:
            messagebox.showerror("Error", "Hubo un error al obtener la proteína")
        return
    url2 = "https://www.rcsb.org/fasta/entry/"+pdb_id+"/download"
    urllib.request.urlretrieve(url2, pdb_id + ".fasta")
    fasta_string = open(pdb_id + ".fasta").read()
    print(fasta_string)
    result_handle = NCBIWWW.qblast("blastp", "pdb", fasta_string.split('\n')[1],
        word_size=7)
    blast = pdb_id + ".xml"
    save_clk = open(blast, "w")
    save_clk.write(result_handle.read())    
    save_clk.close()
    
    blast_records = NCBIXML.parse(open(blast))
    print(blast_records)
    owd = os.getcwd()
    if not os.path.exists("./fasta"):
        os.mkdir("./fasta")
    all_seq_fasta = "./fasta/"+pdb_id+".fasta"
    file = open(all_seq_fasta, "w")
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            aligned_sequence = ''
            aligned_matches = ''
            for hsp in alignment.hsps:
                aligned_sequence = aligned_sequence + hsp.sbjct
                aligned_matches = aligned_matches + hsp.match
            file.writelines(alignment.title)
            file.writelines("\n")
            file.writelines(aligned_sequence)
            file.writelines("\n")
            print(alignment.title.split('>')[0])
            print(aligned_sequence)
            print(aligned_matches)
    file.close()
            


myLabel = Label(root, text="Ingresar codigo PDB")
myLabel.pack()

myTextbox = Entry(root, width=50)
myTextbox.pack()

myButton = Button(root, text="Procesar", command = fasta)
myButton.pack()

root.mainloop()