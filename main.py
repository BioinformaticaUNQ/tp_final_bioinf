from tkinter import *
from tkinter import messagebox, ttk
import urllib.request
import os
from services import pymol_service,blast_service,clustal_service,logomaker_service,dssp_service
from PIL import Image,ImageTk
import tarfile
import numpy as np
from datetime import datetime
from pandas import DataFrame
import logging
import math

#Configuracion de UI con scrollbars

root = Tk()
root.title("Visualizador de regiones conservadas de estructuras homólogas a distintos niveles")
root.geometry("800x600")

def escape_full_screen():
    root.attributes("-fullscreen", False)
root.bind("<Escape>",escape_full_screen())
pdbs_to_process =[]

evalue = DoubleVar(value=0.004)
coverage = IntVar(value=80)

mainFrame = Frame(root)
mainFrame.pack(fill = BOTH, expand = 1)
myCanvas = Canvas(mainFrame)
myCanvas.pack(side = LEFT, fill = BOTH, expand = 1)
myScrollBar = ttk.Scrollbar(mainFrame, orient = VERTICAL, command = myCanvas.yview)
myScrollBarHorizontal = ttk.Scrollbar(mainFrame, orient = HORIZONTAL, command = myCanvas.xview)
myScrollBarHorizontal.pack(side=BOTTOM,fill= X)
myScrollBar.pack(side = RIGHT, fill = Y)
myCanvas.configure(yscrollcommand = myScrollBar.set)
myCanvas.bind('<Configure>', lambda e : myCanvas.configure(scrollregion = myCanvas.bbox("all"),yscrollcommand=myScrollBar.set, xscrollcommand=myScrollBarHorizontal.set))
secondFrame = Frame(myCanvas)
myCanvas.create_window((0, 0), window = secondFrame, anchor = "nw")

progress_bar = ttk.Progressbar(secondFrame,orient=HORIZONTAL,maximum=100)
search_label = Label(secondFrame,text="Por favor, espere...")
btns = []
lbls = []

#Descarga base de datos PDB en la primera ejecucion
if not os.path.exists("./db"):
        os.mkdir("./db")
if not os.path.isfile("./db/pdbaa.pdb"):
    url = "https://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz"
    try:
        if not os.path.isfile("./db/db.tar.gz"):
            print("Descargando base de datos...")
            logging.info("Descargando base de datos...")
            urllib.request.urlretrieve(url, "./db/db.tar.gz")
        tar = tarfile.open("./db/db.tar.gz")
        tar.extractall(path="./db")
        tar.close()
    except Exception as e:
        logging.info("Hubo un error al obtener la base de datos pdb")
        messagebox.showerror("Error", "Hubo un error al obtener la base de datos pdb")
        exit()

#Descarga pdb del pdbid ingresado. Valida que sea correcto
def getPDB():
    input_path = "./ejecucion-" + datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
    os.mkdir(input_path)
    logging.basicConfig(filename=input_path + "/info.log",level=logging.INFO)
    pdb_id = str(pdbTextbox.get())
    url = "https://files.rcsb.org/download/"+ pdb_id +".pdb"
    try:
        urllib.request.urlretrieve(url, pdb_id + ".pdb")
        logging.info("Obtenido pdb de: " + url)
    except Exception as e:
        logging.info(str(e))
        print(str(e))
        if(str(e) == "HTTP Error 404: Not Found"):
            logging.info("Código pdb inválido")
            messagebox.showerror("Error", "Código pdb inválido")
        else:
            logging.info("Hubo un error al obtener la proteína")
            messagebox.showerror("Error", "Hubo un error al obtener la proteína")
        return

    getFASTA(pdb_id,input_path)

#Descarga fasta del pdb ingresado
def getFASTA(pdb_id,input_path):
    url2 = "https://www.rcsb.org/fasta/entry/"+pdb_id+"/download"
    urllib.request.urlretrieve(url2, input_path + "/" + pdb_id + ".fasta")
    logging.info("Obtenido fasta de: " + url2)
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
    lbls.append(numberOfSequenceLabel)
    putSequencesLabelAndButtonsToChoose(pdb_id, dataAndSequencesMap, rna,input_path)

def putSequencesLabelAndButtonsToChoose(pdb_id, dataAndSequencesMap, rna,input_path):
    selectSequenceLabel = Label(secondFrame, text = "Seleccione la cadena a procesar")
    lbls.append(selectSequenceLabel)
    selectSequenceLabel.pack()

    index = 0

    for button in btns:
        button.destroy()

    for label in lbls:
        label.destroy()

    for k, v in dataAndSequencesMap.items():
        index += 1
        chain_label = Label(secondFrame, text="Cadena " + str(index))
        seq_label = Label(secondFrame,text = v)
        lbls.append(chain_label)
        lbls.append(seq_label)
        chain_label.pack(pady = (30, 0))
        seq_label.pack(pady = (30, 0))
        btn = Button(secondFrame, text = "Ejecutar cadena "+str(index), command = lambda k = k, v = v : blast_query(pdb_id, k, v,input_path))
        btn.pack(pady = 5)
        btns.append(btn)

    if len(rna) > 1:
        rnaLabel = Label(secondFrame, text = "Se omitieron " + str(len(rna)) + " secuencias de RNA")
        lbls.append(rnaLabel)
        rnaLabel.pack(pady = (30, 0))
    elif len(rna) == 1:
        rnaLabel = Label(secondFrame, text = "Se omitió 1 secuencia de RNA")
        lbls.append(rnaLabel)
        rnaLabel.pack(pady = (30, 0))

#Busqueda de proteinas homologas con blastp
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
    logging.info("El resto de los valores son estandares de blastp")
    logging.info("https://biopython.readthedocs.io/en/latest/chapter_blast.html")
    if(coverage.get() < 0 or coverage.get() > 100):
        messagebox.showerror("Error", "El pocentaje de coverage debe estar entre 0 y 100")
        return
    if(evalue.get() < 0 or evalue.get() > 0.5):
        messagebox.showerror("Error", "El evalue esperado debe estar entre 0 y 0.5")
        return
    fasta_seq = blast_service.blastp_query(pdb_id,evalue.get(),coverage.get(),data,sequence,input_path)
    define_progress(25)
    align_and_generate_structures(pdb_id, fasta_seq, sequence, data,input_path)
    
#Alineamiento de estructura primaria con proteinas homologas usando clustal
def align_and_generate_structures(pdb_id, fasta_seq, sequence, data,input_path):
    logging.info("Se ejecuto un alineamiento multiple de la cadena problema junto con las proteinas homologas obtenidas previamente")
    logging.info("Esto se ejecuto con Clustal Omega, utilizando como parametro un archivo fasta con todas las cadenas y con el resto de valores por default")
    logging.info("http://www.clustal.org/omega/")
    output_path = clustal_service.run_clustal(pdb_id,fasta_seq,input_path)
    define_progress(50)
    generate_alignment_view(output_path,pdb_id,input_path)
    generate_3structure(pdb_id,input_path)

#Genera el alineamiento de las estructuras secundarias.
#Genera graficos con logomaker para mostrar los grados de conservacion de la estructura secundaria y primaria
def generate_alignment_view(outputPath,pdb_id,input_path):
    pdbs = []
    raw_seqs =[]
    with open(outputPath, "r") as f:
        seqText = ""
        for sq in f:
            if '>' in sq:
                if pdb_id not in sq:
                    for pdb in sq.split(">"):
                        if(pdb.split(" ")[0].split("_")[0] != ""):
                            pdbs.append(pdb.split(" ")[0].split("_")[0])
                if seqText != "":
                    raw_seqs.append(seqText)
                    seqText = "" 
                raw_seqs.append(sq)
            else:
                seqText += sq
    seqs = [seq.strip() for seq in raw_seqs if ('#' not in seq) and ('>') not in seq]
    pdbs_set = list(set(pdbs))
    pdbs_set.append(pdb_id)
    pdbs_to_process.extend(pdbs_set[-10:])
    second_structure_fasta = dssp_service.generate_2structures(pdbs_to_process,input_path,pdb_id)
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
    logomaker_service.primary_structure_conservation(seqs,pdb_id,secondFrame,logging,lbls)
    logomaker_service.secondary_structure_conservation(seqs2,pdb_id,secondFrame,logging,lbls)
    define_progress(75)


#Genera alineamiento de las estructuras secundarias con pymol
def generate_3structure(pdb_id,input_path):
    define_progress(90)

    logging.info("Se corre Pymol con las proteínas: " + str(pdbs_to_process))
    logging.info("Se utilizan los comandos fetch , alignto y el formato cartoon para el gráfico")
    logging.info("https://pymol.org/2/")
    pymol_service.generate_3structure_image(pdb_id,pdbs_to_process,input_path)

    logging.info("Se generó un espacio de trabajo Pymol en la carpeta: " + input_path)
    logging.info("Puede utilizarlo para ver en detalle las estructuras alineadas")

    pymol_label = Label(secondFrame,text="Alineamiento de estructuras terciarias. Para verlo en detalle abrir el archivo .pse en la carpeta de la ejecucion actual.")
    lbls.append(pymol_label)
    pymol_label.pack(pady=(0,30))
    loadPymol = Image.open(input_path + "/" + pdb_id + ".png")
    renderPymol = ImageTk.PhotoImage(loadPymol)
    imgPymol = Label(secondFrame,image=renderPymol)
    lbls.append(imgPymol)
    imgPymol.image = renderPymol
    imgPymol.pack(pady = (15, 0)) 

    #Finalizacion de la busqueda
    delete_unused_files(input_path)
    search_label.config(text="Búsqueda finalizada")
    lbls.append(search_label)
    lbls.append(progress_bar)
    progress_bar.stop()
    progress_bar.step(100)

#Actualiza la progressbar al valor indicado
def define_progress(value):
    progress_bar.step(value)
    root.update_idletasks()

#Se eliminan los archivos que ya fueron utilizados y no tienen utilidad para el usuario
def delete_unused_files(input_path):
    root.attributes("-fullscreen", True)
    root.update_idletasks()
    files = os.listdir("./")
    for file in files:
        if(file.endswith(".pdb") or file.endswith(".fasta") or file.endswith(".cif") or file.endswith(".xml")):
            os.remove(os.path.join("./",file))  
    files = os.listdir(input_path)
    for file in files:
        if(file.endswith(".pdb")):
            os.remove(os.path.join(input_path,file)) 


#UI principal

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