from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import urllib.request
import os

secondary_map = {}

def generate_2structures(pdbs_to_process,output_path,pdb_id):
    print(pdbs_to_process)
    pdb_files = []
    for pdbs_id in pdbs_to_process:
        url = "https://files.rcsb.org/download/"+ pdbs_id +".pdb"
        try:
            urllib.request.urlretrieve(url, output_path +"/" + pdbs_id + ".pdb")
            pdb_files.append(output_path +"/" + pdbs_id + ".pdb")
        except Exception as e:
            print(str(e))
    all_seq_fasta = output_path + "/"+pdb_id+".fasta"
    for pdb_file in pdb_files:
        dssp_tuple = dssp_dict_from_pdb_file(pdb_file)
        dssp_dict = dssp_tuple[0]
        #EL PRIMER VALOR DE LA TUPLA ES UN DICCIONARIO (TUPLA KEY, DATA DE LA ESTRUCTURA)
        #EL SEGUNDO VALOR ES LA LISTA DE KEYS QUE SON DEL FORMATO ("CADENA",('',NRO DE RESIDUO,''))
        #LAS CADENAS ESTAN SEPARADAS , POR EJ LA CADENA A SON MUCHAS KEYS TODAS EMPEZANDO CON A PERO CON DISTINTO NUMERO DE RESIDUO
        AChain = []
        for key in dssp_tuple[1]:
            if(key[0] == 'A'):
                #OBTENGO TODAS LAS KEYS DE LA CADENA A
                AChain.append(key)
        seq = ""
        for chainPart in AChain:
            #OBTENGO LAS ESTRUCTURAS DE LA CADENA A
            if(dssp_dict[chainPart][0] == 'a'):
                seq += 'C'
            else:
                seq += dssp_dict[chainPart][0]
        secondary_map[pdb_file.split('/')[2].split('.')[0]] = seq
    return generate_secondary_fasta(get_primary_map(pdb_id),output_path)


def generate_secondary_fasta(primary_map, input_path):
    new_file = open(input_path + '/secondaryFasta.fasta', 'w')
    final_map = {}
    for key in primary_map:
        if key in secondary_map.keys():
            secondary_chain_with_gaps = compare_chains(primary_map[key], secondary_map[key])
            final_map[key] = secondary_chain_with_gaps

    max_len = len(max(final_map.values(),key=len))
    for key,value in final_map.items():
            new_value = value
            while(len(new_value) < max_len):
                new_value += '-'
            new_file.writelines('>' + key + '\n')
            new_file.writelines(new_value + '\n')
    new_file.close()
    return input_path + '/secondaryFasta.fasta'


def get_primary_map(pdb_id):
    seq_map = {}
    with open("./fasta/"+pdb_id+"_aln.fasta", "r") as f:
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

def compare_chains(primary_chain, secondary_chain):
    result = ''
    index = 0
    for char in primary_chain:
        if index < len(secondary_chain):
            if char == '-':
                result += '-'
            else:
                result += secondary_chain[index]
                index += 1

    if len(secondary_chain) > len(primary_chain):
        result += secondary_chain[index : len(secondary_chain)]
    
    return result
        