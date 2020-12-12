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

def escape_full_screen(event):
    root.attributes("-fullscreen", False)
root.bind("<Escape>",escape_full_screen)
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
            urllib.request.urlretrieve(url, "./db/db.tar.gz")
        tar = tarfile.open("./db/db.tar.gz")
        tar.extractall(path="./db")
        tar.close()
    except Exception as e:
        messagebox.showerror("Error", "Hubo un error al obtener la base de datos pdb")
        exit()

#Descarga pdb del pdbid ingresado. Valida que sea correcto
def getPDB():
    input_path = "./ejecucion-" + datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
    os.mkdir(input_path)
    logger = logging.getLogger(input_path)
    logger.setLevel(logging.INFO)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(input_path + "/info.log")
    fh.setLevel(logging.INFO)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    # logging.basicConfig(filename=input_path + "/info.log",filemode='w',level=logging.INFO)
    pdb_id = str(pdbTextbox.get())
    url = "https://files.rcsb.org/download/"+ pdb_id +".pdb"
    try:
        urllib.request.urlretrieve(url, pdb_id + ".pdb")
        logger.info("Obtenido pdb de: " + url)
    except Exception as e:
        logger.info(str(e))
        print(str(e))
        if(str(e) == "HTTP Error 404: Not Found"):
            logger.info("Código pdb inválido")
            messagebox.showerror("Error", "Código pdb inválido")
        else:
            logger.info("Hubo un error al obtener la proteína")
            messagebox.showerror("Error", "Hubo un error al obtener la proteína")
        return

    resetGUI()
    getFASTA(pdb_id,input_path,logger)

def resetGUI():
    for button in btns:
        button.destroy()

    for label in lbls:
        label.destroy()

#Descarga fasta del pdb ingresado
def getFASTA(pdb_id,input_path,logger):
    url2 = "https://www.rcsb.org/fasta/entry/"+pdb_id+"/download"
    urllib.request.urlretrieve(url2, input_path + "/" + pdb_id + ".fasta")
    logger.info("Obtenido fasta de: " + url2)
    getSequencesFromPDB(pdb_id,input_path,logger)

def getSequencesFromPDB(pdb_id,input_path,logger):
    fasta_string = open(input_path + "/" + pdb_id + ".fasta").read()
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

    putNumberOfSequencesLabel(pdb_id, numberOfSequences, dataAndSequencesMap, rna,input_path,logger)

def isRNA(sequence):
    return len("".join(dict.fromkeys(sequence))) == 4

def oddNumber(number):
    return number % 2 != 0

def putNumberOfSequencesLabel(pdb_id, numberOfSequences, dataAndSequencesMap, rna,input_path,logger):
    if (numberOfSequences > 1):
        numberOfSequenceLabel = Label(secondFrame, text = "El código " + pdb_id + " tiene " + str(numberOfSequences) + " secuencias")
        numberOfSequenceLabel.pack(pady = (0, 15))
    else:
        numberOfSequenceLabel = Label(secondFrame, text = "El código " + pdb_id + " tiene 1 sola secuencia")
        numberOfSequenceLabel.pack(pady = (0, 15))
    lbls.append(numberOfSequenceLabel)
    putSequencesLabelAndButtonsToChoose(pdb_id, dataAndSequencesMap, rna,input_path,logger)

def putSequencesLabelAndButtonsToChoose(pdb_id, dataAndSequencesMap, rna,input_path,logger):
    selectSequenceLabel = Label(secondFrame, text = "Cadenas")
    selectSequenceLabel.config(font=("Courier", 20))
    lbls.append(selectSequenceLabel)
    selectSequenceLabel.pack()

    index = 0

    for k, v in dataAndSequencesMap.items():
        index += 1
        chain_label = Label(secondFrame, text="Cadena " + str(index))
        seq_label = Label(secondFrame,text = v)
        lbls.append(chain_label)
        lbls.append(seq_label)
        chain_label.pack(pady = (5, 0))
        seq_label.pack(pady = (10, 0))
        btn = Button(secondFrame, text = "Ejecutar cadena "+str(index), command = lambda k = k, v = v : blast_query(pdb_id, k, v,input_path,logger))
        btn.pack(pady = (15,15))
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
def blast_query(pdb_id, data, sequence,input_path,logger):
    search_label.config(text = "Por favor, espere...")
    search_label.pack(pady=(0,30))
    progress_bar.pack(pady=(0,30))
    progress_bar.start()
    root.update_idletasks()

    logger.info("Eligio la secuencia " + sequence)    
    logger.info("Se ejecuta la busqueda de proteinas homologas con los siguientes parametros:")
    logger.info("Porcentaje de identidad: >40%")
    logger.info("Porcentaje de coverage: >" + str(coverage.get()))
    logger.info("eValue: " + str(evalue.get()))
    logger.info("El resto de los valores son estandares de blastp")
    logger.info("https://biopython.readthedocs.io/en/latest/chapter_blast.html")
    if(coverage.get() < 0 or coverage.get() > 100):
        messagebox.showerror("Error", "El pocentaje de coverage debe estar entre 0 y 100")
        return
    if(evalue.get() < 0 or evalue.get() > 0.5):
        messagebox.showerror("Error", "El evalue esperado debe estar entre 0 y 0.5")
        return
    fasta_seq = blast_service.blastp_query(pdb_id,evalue.get(),coverage.get(),data,sequence,input_path,logger)
    define_progress(25)
    align_and_generate_structures(pdb_id, fasta_seq, sequence, data,input_path,logger)
    
#Alineamiento de estructura primaria con proteinas homologas usando clustal
def align_and_generate_structures(pdb_id, fasta_seq, sequence, data,input_path,logger):
    logger.info("Se ejecuto un alineamiento multiple de la cadena problema junto con las proteinas homologas obtenidas previamente")
    logger.info("Esto se ejecuto con Clustal Omega, utilizando como parametro un archivo fasta con todas las cadenas y con el resto de valores por default")
    logger.info("http://www.clustal.org/omega/")
    output_path = clustal_service.run_clustal(pdb_id,fasta_seq,input_path,logger)
    if(output_path is None):
        return
    define_progress(50)
    generate_alignment_view(output_path,pdb_id,input_path,logger)
    generate_3structure(pdb_id,input_path,logger)

#Genera el alineamiento de las estructuras secundarias.
#Genera graficos con logomaker para mostrar los grados de conservacion de la estructura secundaria y primaria
def generate_alignment_view(outputPath,pdb_id,input_path,logger):
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
                    raw_seqs.append(seqText.replace('\n', ""))
                    seqText = "" 
                raw_seqs.append(sq)
            else:
                seqText += sq
    seqs = [seq.strip() for seq in raw_seqs if ('#' not in seq) and ('>') not in seq]
    pdbs_set = list(set(pdbs))
    pdbs_set.append(pdb_id)
    pdbs_to_process.extend(pdbs_set[-10:])
    second_structure_fasta = dssp_service.generate_2structures(pdbs_to_process,input_path,pdb_id,logger)
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
    logomaker_service.primary_structure_conservation(seqs,pdb_id,secondFrame,logger,lbls,input_path)
    logomaker_service.secondary_structure_conservation(seqs2,pdb_id,secondFrame,logger,lbls,input_path)
    define_progress(75)


#Genera alineamiento de las estructuras secundarias con pymol
def generate_3structure(pdb_id,input_path,logger):
    define_progress(90)

    logger.info("Se corre Pymol con las proteínas: " + str(pdbs_to_process))
    logger.info("Se utilizan los comandos fetch , alignto y el formato cartoon para el gráfico")
    logger.info("https://pymol.org/2/")
    pymol_service.generate_3structure_image(pdb_id,pdbs_to_process,input_path,logger)

    logger.info("Se generó un espacio de trabajo Pymol en la carpeta: " + input_path)
    logger.info("Puede utilizarlo para ver en detalle las estructuras alineadas")
    pymol_label = Label(secondFrame,text="Alineamiento de estructura terciaria: ")
    pymol_label.config(font=("Verdana",20))
    pymol_label_pse = Label(secondFrame,text="(Para verlo en detalle abrir el archivo .pse, guardado en la carpeta de ejecucion actual, con Pymol)")
    pymol_label_pse.config(font=("Verdana",15))
    lbls.append(pymol_label)
    lbls.append(pymol_label_pse)
    pymol_label.pack(anchor=CENTER)
    pymol_label_pse.pack(anchor=CENTER)
    loadPymol = Image.open(input_path + "/" + pdb_id + ".png")
    renderPymol = ImageTk.PhotoImage(loadPymol)
    imgPymol = Label(secondFrame,image=renderPymol)
    lbls.append(imgPymol)
    imgPymol.image = renderPymol
    imgPymol.pack(anchor=CENTER) 

    #Finalizacion de la busqueda
    delete_unused_files(input_path)
    search_label.config(text="Búsqueda finalizada")
    search_label.pack_forget()
    progress_bar.pack_forget()
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

titleLabel = Label(secondFrame, text="Visualizador de regiones conservadas de estructuras homólogas a distintos niveles")
titleLabel.pack(anchor=CENTER, pady = (0,15))
titleLabel.config(font=("Verdana",18)) 

pdbLabel = Label(secondFrame, text="Ingrese un código PDB: ")
pdbLabel.pack(anchor=CENTER)
pdbLabel.config(font=("Verdana",10))

pdbTextbox = Entry(secondFrame, width=30)
pdbTextbox.pack()

evalueLabel = Label(secondFrame, text="eValue: ")
evalueLabel.pack(anchor=CENTER)
evalueLabel.config(font=("Verdana",10))

evalueTextbox = Entry(secondFrame, width=30, textvariable=evalue)
evalueTextbox.pack()

coverageLabel = Label(secondFrame, text="Coverage: ")
coverageLabel.pack(anchor=CENTER)
coverageLabel.config(font=("Verdana",10))

coverageTextbox = Entry(secondFrame, width=30,textvariable=coverage)
coverageTextbox.pack()

myButton = Button(secondFrame, text="Procesar", command = getPDB)
myButton.pack(anchor=CENTER, pady = (0, 15))

root.mainloop()