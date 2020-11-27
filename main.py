from tkinter import *
from tkinter import messagebox
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO
import urllib.request
import pymol
from Bio.Align.Applications import ClustalOmegaCommandline 
import logomaker as lm
import matplotlib.pyplot as plt
import os
from PIL import Image,ImageTk

root = Tk()
root.title("TP FINAL")
root.geometry("1920x1080")

def blast_query(fasta_string,pdb_id):
    result_handle = NCBIWWW.qblast("blastp", "pdb", fasta_string.split('\n')[1],alignments=10,hitlist_size=10)
    blast = pdb_id + ".xml"
    with open(blast, "w") as out_handle:
        out_handle.write(result_handle.read())

    result_handle.close()
    return blast

def generate_3structure(pdb_id):
    pymol.cmd.load(pdb_id + ".pdb", pdb_id)
    pymol.cmd.disable("all")
    pymol.cmd.enable(pdb_id)
    print(pymol.cmd.get_names())
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)
    pymol.cmd.color('red', 'ss h')
    pymol.cmd.color('yellow', 'ss s')
    pymol.cmd.png("%s.png"%(pdb_id))
    #pymol.cmd.quit()
    loadPymol = Image.open(pdb_id + ".png")
    renderPymol = ImageTk.PhotoImage(loadPymol)
    imgPymol = Label(root,image=renderPymol)
    imgPymol.image = renderPymol
    imgPymol.place(x=0,y=400)

def generate_alignment_view(outputPath,pdb_id):
    raw_seqs =[]
    with open(outputPath, "r") as f:
        seqText = ""
        for sq in f:
            if '>' in sq:
                if seqText != "":
                    raw_seqs.append(seqText)
                    seqText = "" 
                raw_seqs.append(sq)
            else:
                seqText += sq
    seqs = [seq.strip() for seq in raw_seqs if ('#' not in seq) and ('>') not in seq]

    counts_mat = lm.alignment_to_matrix(seqs)
    counts_mat.head()
    crp_logo = lm.Logo(counts_mat, font_name = 'Arial Rounded MT Bold')

    # style using Axes methods
    crp_logo.ax.xaxis.set_ticks_position('none')
    crp_logo.ax.xaxis.set_tick_params(pad=-1)
    plt.savefig(pdb_id + "_aln.png")
    load = Image.open(pdb_id + "_aln.png")
    render = ImageTk.PhotoImage(load)
    img = Label(root,image=render)
    img.image = render
    img.place(x=0,y=80)

def fasta():
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
    
    blast = blast_query(fasta_string,pdb_id)

    blast_records = NCBIXML.parse(open(blast))

    owd = os.getcwd()
    if not os.path.exists("./fasta"):
        os.mkdir("./fasta")
    all_seq_fasta = "./fasta/"+pdb_id+".fasta"
    file = open(all_seq_fasta, "w")
    file.writelines(fasta_string.split('\n')[0] + '\n')
    file.writelines(fasta_string.split('\n')[1] + '\n')
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            file.writelines(">" + alignment.hit_id)
            file.writelines("\n")
            for hsp in alignment.hsps:
                sbjct_no_gaps = hsp.sbjct.replace("-","")
                file.writelines(sbjct_no_gaps)
                file.writelines("\n")
    file.close()
    outputPath = "./fasta/"+pdb_id+"_aln.fasta"
    clustalomega_cline = ClustalOmegaCommandline(infile = all_seq_fasta, outfile = outputPath,force = True)
    clustalomega_cline()
    
    generate_alignment_view(outputPath,pdb_id)

    generate_3structure(pdb_id)
    
            


myLabel = Label(root, text="Ingresar codigo PDB")
myLabel.pack()

myTextbox = Entry(root, width=50)
myTextbox.pack()

myButton = Button(root, text="Procesar", command = fasta)
myButton.pack()

root.mainloop()