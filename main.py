from tkinter import *
from tkinter import messagebox, ttk
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SearchIO
import urllib.request
import pymol
from Bio.Align.Applications import ClustalOmegaCommandline 
import logomaker as lm
import matplotlib.pyplot as plt
import os
from PIL import Image,ImageTk
import tarfile
import numpy as np
from pandas import DataFrame

root = Tk()
root.title("TP FINAL")
root.geometry("1920x1080")

#SCROLLBAR
mainFrame = Frame(root)
mainFrame.pack(fill = BOTH, expand = 1)
myCanvas = Canvas(mainFrame)
myCanvas.pack(side = LEFT, fill = BOTH, expand = 1)
myScrollBar = ttk.Scrollbar(mainFrame, orient = VERTICAL, command = myCanvas.yview)
myScrollBar.pack(side = RIGHT, fill = Y)
myCanvas.configure(yscrollcommand = myScrollBar.set)
myCanvas.bind('<Configure>', lambda e : myCanvas.configure(scrollregion = myCanvas.bbox("all")))
secondFrame = Frame(myCanvas)
myCanvas.create_window((0, 0), window = secondFrame, anchor = "nw")

#DESCARGA BDD
if not os.path.exists("./db"):
        os.mkdir("./db")
if not os.path.isfile("./db/pdbaa.pdb"):
    url = "https://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz"
    try:
        if not os.path.isfile("./db/db.tar.gz"):
            print("Descargando base de datos...")
            urllib.request.urlretrieve(url, "./db/db.tar.gz")
        tar = tarfile.open("./db/db.tar.gz")
        tar.extractall(path="./db")
        tar.close()
    except Exception as e:
        messagebox.showerror("Error", "Hubo un error al obtener la base de datos pdb")
        exit()

def getPDB():
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

    getFASTA(pdb_id)

def getFASTA(pdb_id):
    url2 = "https://www.rcsb.org/fasta/entry/"+pdb_id+"/download"
    urllib.request.urlretrieve(url2, pdb_id + ".fasta")
    getSequencesFromPDB(pdb_id)

def getSequencesFromPDB(pdb_id):
    fasta_string = open(pdb_id + ".fasta").read()
    sequences = []
    data = []
    rna = []
    index = 0

    for line in fasta_string.splitlines():
        if (oddNumber(index)):
            if not (isRNA(line)):
                sequences.append(line)
            else :
                rna.append(line)
        else:
            data.append(line)
        index += 1
    
    numberOfSequences = len(sequences)
    dataAndSequencesMap = dict(zip(data,sequences))

    putNumberOfSequencesLabel(pdb_id, numberOfSequences, dataAndSequencesMap, rna)

def isRNA(sequence):
    return len("".join(dict.fromkeys(sequence))) == 4

def oddNumber(number):
    return number % 2 != 0

def putNumberOfSequencesLabel(pdb_id, numberOfSequences, dataAndSequencesMap, rna):
    if (numberOfSequences > 1):
        numberOfSequenceLabel = Label(secondFrame, text = "El código " + pdb_id + " tiene " + str(numberOfSequences) + " secuencias")
        numberOfSequenceLabel.pack(pady = 5)
    else:
        numberOfSequenceLabel = Label(secondFrame, text = "El código " + pdb_id + " tiene 1 sola secuencia")
        numberOfSequenceLabel.pack(pady = 5)

    putSequencesLabelAndButtonsToChoose(pdb_id, dataAndSequencesMap, rna)

def putSequencesLabelAndButtonsToChoose(pdb_id, dataAndSequencesMap, rna):
    selectSequenceLabel = Label(secondFrame, text = "Seleccione la secuencia a procesar")
    selectSequenceLabel.pack()

    for k, v in dataAndSequencesMap.items():
        Button(secondFrame, text = v, command = lambda k = k, v = v : blast_query(pdb_id, k, v)).pack(pady = 5)

    if len(rna) > 1:
        rnaLabel = Label(secondFrame, text = "Se omitieron " + str(len(rna)) + " secuencias de RNA")
        rnaLabel.pack(pady = (30, 0))
    elif len(rna) == 1:
        rnaLabel = Label(secondFrame, text = "Se omitió 1 secuencia de RNA")
        rnaLabel.pack(pady = (30, 0))

def blast_query(pdb_id, data, sequence):
    print("Eligio la secuencia " + sequence)    
    blast = pdb_id + ".xml"
    cline = NcbiblastpCommandline(query=pdb_id + '.fasta', db="./db/pdbaa",
                              evalue=0.001, out=blast, outfmt=5,qcov_hsp_perc=80)
    cline()
    
    runClustal(pdb_id, blast, sequence, data)
    
def runClustal(pdb_id, blast, sequence, data):
    blast_records = NCBIXML.parse(open(blast))
    if not os.path.exists("./fasta"):
        os.mkdir("./fasta")
    
    all_seq_fasta = "./fasta/"+pdb_id+".fasta"
    file = open(all_seq_fasta, "w")
    file.writelines(data + '\n')
    file.writelines(sequence + '\n')

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
    counts_mat_list = np.array_split(counts_mat, 4)
    
    for df in counts_mat_list:
        crp_logo = lm.Logo(df, font_name = 'Arial Rounded MT Bold')

        # style using Axes methods
        crp_logo.ax.xaxis.set_ticks_position('none')
        crp_logo.ax.xaxis.set_tick_params(pad=-1)
        plt.savefig(pdb_id + "_aln.png")
        load = Image.open(pdb_id + "_aln.png")
        render = ImageTk.PhotoImage(load)
        img = Label(secondFrame,image=render)
        img.image = render
        img.pack(pady = (30, 0))

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
    imgPymol = Label(secondFrame,image=renderPymol)
    imgPymol.image = renderPymol
    imgPymol.pack(pady = (15, 0))

myLabel = Label(secondFrame, text="Ingresar código PDB")
myLabel.pack()

myTextbox = Entry(secondFrame, width=50)
myTextbox.pack()

myButton = Button(secondFrame, text="Procesar", command = getPDB)
myButton.pack(pady = (0, 30))

root.mainloop()