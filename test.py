from tkinter import *
from tkinter import messagebox
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO
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
    result_handle = NCBIWWW.qblast("blastp", "pdb", fasta_string)
    blast = pdb_id + ".xml"
    with open(blast, "w") as out_handle:
        out_handle.write(result_handle.read())

    result_handle.close()
    blast_records = NCBIXML.parse(open(blast))
    print(blast_records)
    owd = os.getcwd()
    if not os.path.exists("./fasta"):
        os.mkdir("./fasta")
    all_seq_fasta = "./fasta/"+pdb_id+".fasta"
    file = open(all_seq_fasta, "w")
    # for qresult in SearchIO.parse(pdb_id + ".xml", "blast-xml"):
    #     for hit in qresult.hits:
    #         print(hit)
    #         for hsp in hit:
    #             print(hsp)
    #             file.writelines(">" +hit.description)
    #             file.writelines("\n")
    #             #file.writelines(hsp.hit_string)
    #             file.writelines("\n")

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            file.writelines(">" + alignment.hit_id)
            file.writelines("\n")
            for hsp in alignment.hsps:
                sbjct_no_gaps = hsp.sbjct.replace("-","")
                file.writelines(sbjct_no_gaps)
                file.writelines("\n")
    file.close()
            


myLabel = Label(root, text="Ingresar codigo PDB")
myLabel.pack()

myTextbox = Entry(root, width=50)
myTextbox.pack()

myButton = Button(root, text="Procesar", command = fasta)
myButton.pack()

root.mainloop()