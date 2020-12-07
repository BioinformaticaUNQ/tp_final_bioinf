from tkinter import *
from tkinter import messagebox, ttk
import urllib.request
import logomaker as lm
import matplotlib.pyplot as plt
import os
from services import pymol_service,blast_service,clustal_service,logomaker_service,dssp_service
from PIL import Image,ImageTk
import tarfile
import numpy as np
from datetime import datetime
from pandas import DataFrame
import logging
from threading import Thread
import math

root = Tk()
root.title("TP FINAL")
root.geometry("800x600")
pdbs_to_process =[]

evalue = DoubleVar(value=0.004)
coverage = IntVar(value=80)
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

progress_bar = ttk.Progressbar(secondFrame,orient=HORIZONTAL,maximum=100)
search_label = Label(secondFrame,text="Por favor, espere...")

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
    input_path = "./ejecucion-" + datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
    os.mkdir(input_path)
    logging.basicConfig(filename=input_path + "/info.log",level=logging.INFO)
    pdb_id = str(pdbTextbox.get())
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

    getFASTA(pdb_id,input_path)

def getFASTA(pdb_id,input_path):
    url2 = "https://www.rcsb.org/fasta/entry/"+pdb_id+"/download"
    urllib.request.urlretrieve(url2, pdb_id + ".fasta")
    getSequencesFromPDB(pdb_id,input_path)

def getSequencesFromPDB(pdb_id,input_path):
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

    putNumberOfSequencesLabel(pdb_id, numberOfSequences, dataAndSequencesMap, rna,input_path)

def isRNA(sequence):
    return len("".join(dict.fromkeys(sequence))) == 4

def oddNumber(number):
    return number % 2 != 0

def putNumberOfSequencesLabel(pdb_id, numberOfSequences, dataAndSequencesMap, rna,input_path):
    if (numberOfSequences > 1):
        numberOfSequenceLabel = Label(secondFrame, text = "El código " + pdb_id + " tiene " + str(numberOfSequences) + " secuencias")
        numberOfSequenceLabel.pack(pady = 5)
    else:
        numberOfSequenceLabel = Label(secondFrame, text = "El código " + pdb_id + " tiene 1 sola secuencia")
        numberOfSequenceLabel.pack(pady = 5)

    putSequencesLabelAndButtonsToChoose(pdb_id, dataAndSequencesMap, rna,input_path)

def putSequencesLabelAndButtonsToChoose(pdb_id, dataAndSequencesMap, rna,input_path):
    selectSequenceLabel = Label(secondFrame, text = "Seleccione la cadena a procesar")
    selectSequenceLabel.pack()

    index = 0
    for k, v in dataAndSequencesMap.items():
        index += 1
        chain_label = Label(secondFrame, text="Cadena " + str(index))
        seq_label = Label(secondFrame,text = v)
        chain_label.pack(pady = (30, 0))
        seq_label.pack(pady = (30, 0))
        Button(secondFrame, text = "Ejecutar cadena "+str(index), command = lambda k = k, v = v : start_query(pdb_id, k, v,input_path)).pack(pady = 5)

    if len(rna) > 1:
        rnaLabel = Label(secondFrame, text = "Se omitieron " + str(len(rna)) + " secuencias de RNA")
        rnaLabel.pack(pady = (30, 0))
    elif len(rna) == 1:
        rnaLabel = Label(secondFrame, text = "Se omitió 1 secuencia de RNA")
        rnaLabel.pack(pady = (30, 0))


def start_query(pdb_id, data, sequence,input_path):
    t= Thread(target=blast_query(pdb_id, data, sequence,input_path))
    t.start()

def blast_query(pdb_id, data, sequence,input_path):
    search_label.pack(pady=(0,30))
    progress_bar.pack(pady=(0,30))
    progress_bar.start()
    root.update_idletasks()
    logging.info("Eligio la secuencia " + sequence)    
    logging.info("Se ejecuta la busqueda de proteinas homologas con los siguientes parametros:")
    logging.info("Porcentaje de identidad: >40%")
    logging.info("Porcentaje de coverage: >" + str(coverage.get()))
    logging.info("eValue: " + str(evalue.get()))
    logging.info("El resto de los valores son estandares de blastp") #Poner link de doc
    fasta_seq = blast_service.blastp_query(pdb_id,evalue.get(),coverage.get(),data,sequence)
    progress_bar.step(25)
    root.update_idletasks()
    align_and_generate_structures(pdb_id, fasta_seq, sequence, data,input_path)
    
def align_and_generate_structures(pdb_id, fasta_seq, sequence, data,input_path):
    output_path = clustal_service.run_clustal(pdb_id,fasta_seq,input_path)
    progress_bar.step(50)
    root.update_idletasks()
    generate_alignment_view(output_path,pdb_id,input_path)
    generate_3structure(pdb_id,input_path)

def generate_alignment_view(outputPath,pdb_id,input_path):
    pdbs = []
    raw_seqs =[]
    with open(outputPath, "r") as f:
        seqText = ""
        for sq in f:
            if '>' in sq:
                if 'pdb' in sq:
                    pdbs.append(sq.split("|")[1])
                if seqText != "":
                    raw_seqs.append(seqText)
                    seqText = "" 
                raw_seqs.append(sq)
            else:
                seqText += sq
    seqs = [seq.strip() for seq in raw_seqs if ('#' not in seq) and ('>') not in seq]
    pdbs.append(pdb_id)
    pdbs_to_process.extend(pdbs[-10:])
    second_structure_fasta = dssp_service.generate_2structures(pdbs_to_process,input_path,pdb_id)
    # primary_map = get_primary_map(input_path + '/' + pdb_id + '_aln.fasta', pdb_id)
    # secondary_fasta = dssp_service.generate_secondary_fasta(primary_map, input_path)
    raw_seqs2 =[]
    with open(second_structure_fasta, "r") as f:
        seqText = ""
        for sq in f:
            if '>' in sq:
                if seqText != "":
                    raw_seqs2.append(seqText)
                    seqText = "" 
                raw_seqs2.append(sq)
            else:
                seqText += sq
    seqs2 = [seq.strip() for seq in raw_seqs2 if ('#' not in seq) and ('>') not in seq]
    
    
    counts_mat = lm.alignment_to_matrix(seqs)
    divider = len(seqs[0]) / 35
    counts_mat_list = np.array_split(counts_mat, math.ceil(divider))
    
    counts_mat2 = lm.alignment_to_matrix(seqs2)
    divider2= len(seqs2[0]) / 40
    counts_mat_list2 = np.array_split(counts_mat2, math.ceil(divider2))

    alignment_label = Label(secondFrame,text="Alineamiento de estructura primaria")
    alignment_label.pack(pady=(0,30))
    for df in counts_mat_list:
        crp_logo = lm.Logo(df,color_scheme='skylign_protein')

        # style using Axes methods
        crp_logo.ax.xaxis.set_ticks_position('none')
        crp_logo.ax.xaxis.set_tick_params(pad=-1)
        plt.savefig(input_path + "/" + pdb_id + "_aln.png")
        load = Image.open(input_path + "/" + pdb_id + "_aln.png")
        render = ImageTk.PhotoImage(load)
        img = Label(secondFrame,image=render)
        img.image = render
        img.pack(pady = (30, 0))

    alignment_label2 = Label(secondFrame,text="Alineamiento de estructura secundaria")
    alignment_label2.pack(pady=(0,30))
    for df in counts_mat_list2:
        crp_logo = lm.Logo(df)
        # style using Axes methods
        crp_logo.ax.xaxis.set_ticks_position('none')
        crp_logo.ax.xaxis.set_tick_params(pad=-1)
        plt.savefig(pdb_id + "_aln_secondary.png")
        load = Image.open(pdb_id + "_aln_secondary.png")
        render = ImageTk.PhotoImage(load)
        img = Label(secondFrame,image=render)
        img.image = render
        img.pack(pady = (30, 0))
    progress_bar.step(75)
    root.update_idletasks()

def get_primary_map(filePath, pdb_id):
    seq_map = {}
    with open(filePath, "r") as f:
            seqText = ""
            key = ""

            for sq in f:
                if '>' in sq:
                    if seqText != "":
                        seq_map[key] = seqText
                        seqText = "" 
                    if pdb_id not in sq:
                        seq_map[sq.split("|")[1]] = ''
                        key = sq.split("|")[1]
                    else:
                        seq_map[pdb_id] = ''
                        key = pdb_id

                else:
                    seqText += sq.replace('\n',"")
            seq_map[key] = seqText
    return seq_map

def generate_3structure(pdb_id,input_path):
    progress_bar.step(90)
    root.update_idletasks()
    pymol_service.generate_3structure_image(pdb_id,pdbs_to_process,input_path)
    pymol_label = Label(secondFrame,text="Alineamiento de estructuras terciarias. Para verlo en detalle abrir el archivo .pse en la carpeta de la ejecucion actual.")
    pymol_label.pack(pady=(0,30))
    loadPymol = Image.open(input_path + "/" + pdb_id + ".png")
    renderPymol = ImageTk.PhotoImage(loadPymol)
    imgPymol = Label(secondFrame,image=renderPymol)
    imgPymol.image = renderPymol
    imgPymol.pack(pady = (15, 0))
    files = os.listdir("./")
    for file in files:
        if(file.endswith(".pdb") or file.endswith(".fasta") or file.endswith(".cif") or file.endswith(".xml")):
            os.remove(os.path.join("./",file))  
    files = os.listdir(input_path)
    for file in files:
        if(file.endswith(".pdb") or file.endswith(".png")):
            os.remove(os.path.join(input_path,file))  
    search_label.config(text="Busqueda finalizada")
    progress_bar.stop()

pdbLabel = Label(secondFrame, text="Ingresar código PDB")
pdbLabel.pack()

pdbTextbox = Entry(secondFrame, width=50)
pdbTextbox.pack()

evalueLabel = Label(secondFrame, text="eValue")
evalueLabel.pack()

evalueTextbox = Entry(secondFrame, width=25, textvariable=evalue)
evalueTextbox.pack()

coverageLabel = Label(secondFrame, text="coverage")
coverageLabel.pack()

coverageTextbox = Entry(secondFrame, width=25,textvariable=coverage)
coverageTextbox.pack()

myButton = Button(secondFrame, text="Procesar", command = getPDB)
myButton.pack(pady = (0, 30))

root.mainloop()