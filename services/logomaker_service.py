import logomaker as lm
from tkinter import *
import numpy as np
import math
from tkinter import messagebox, ttk
import matplotlib.pyplot as plt
from PIL import Image,ImageTk
import logging

def primary_structure_conservation(seqs,pdb_id,frame,logger,lbls):
    counts_mat = lm.alignment_to_matrix(seqs)
    divider = len(seqs[0]) / 35
    logger.info("Se divide las cadenas de las estructuras primarias en " + str(math.ceil(divider)) + " partes para la generacion de graficos con LogoMaker")
    logger.info("Se utiliza el color scheme skylign_protein")
    logger.info("https://academic.oup.com/bioinformatics/article/36/7/2272/5671693")
    counts_mat_list = np.array_split(counts_mat, math.ceil(divider))
    alignment_label = Label(frame,text="Alineamiento de estructura primaria: ")
    alignment_label.config(font=("Verdana",20))
    lbls.append(alignment_label)
    alignment_label.pack(anchor=CENTER)
    for df in counts_mat_list:
        crp_logo = lm.Logo(df,color_scheme='skylign_protein')
        # style using Axes methods
        crp_logo.ax.xaxis.set_ticks_position('none')
        crp_logo.ax.xaxis.set_tick_params(pad=-1)
        plt.savefig( pdb_id + "_aln.png")
        load = Image.open(pdb_id + "_aln.png")
        render = ImageTk.PhotoImage(load)
        img = Label(frame,image=render)
        lbls.append(img)
        img.image = render
        img.pack(pady = (50, 0))

def secondary_structure_conservation(seqs2,pdb_id,frame,logger,lbls):
    counts_mat2 = lm.alignment_to_matrix(seqs2)
    divider2= len(seqs2[0]) / 40
    logger.info("Se divide las cadenas de las estructuras secundarias en " + str(math.ceil(divider2)) + " partes para la generacion de graficos con LogoMaker")
    logger.info("Se utiliza el color scheme skylign_protein")
    logger.info("https://academic.oup.com/bioinformatics/article/36/7/2272/5671693")

    counts_mat_list2 = np.array_split(counts_mat2, math.ceil(divider2))
    alignment_label2 = Label(frame,text="Alineamiento de estructura secundaria: ")
    alignment_label2.config(font=("Verdana",20))
    lbls.append(alignment_label2)
    alignment_label2.pack(pady=(0,30))
    for df in counts_mat_list2:
        crp_logo = lm.Logo(df,color_scheme='skylign_protein')
        # style using Axes methods
        crp_logo.ax.xaxis.set_ticks_position('none')
        crp_logo.ax.xaxis.set_tick_params(pad=-1)
        plt.savefig(pdb_id + "_aln_secondary.png")
        load = Image.open(pdb_id + "_aln_secondary.png")
        render = ImageTk.PhotoImage(load)
        img = Label(frame,image=render)
        lbls.append(img)
        img.image = render
        img.pack(pady = (30, 0))