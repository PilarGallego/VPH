#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq
from collections import defaultdict
import re
import logging
from Bio.SeqUtils import GC, GC123
from Bio.Alphabet import generic_dna
from collections import Counter
from Bio.Seq import MutableSeq
import pandas as pd
from scipy import stats
from scipy import mean as sc_m
import itertools
import pandas as pd
import researchpy as rp
from scipy import stats
import seaborn as sns
import statsmodels.api as sm
from statsmodels.formula.api import ols
import matplotlib.pyplot as plt
from statsmodels.stats import anova 
from matplotlib import cm 
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.libqsturng import psturng
from tabulate import tabulate

get_ipython().run_line_magic('matplotlib', 'inline')


# In[6]:


logging.basicConfig(filename = "avisos_cambios.log", filemode = "w", 
                    format = "%(name)s - %(levelname)s - %(asctime)s - %(message)s",
                    level = logging.INFO)


# In[7]:


#ruta a la carpeta en la que se encuentran almacenados los archivos genbank
directorio = "./GenomasHPV/"
list_files = []
for root, directories, files in os.walk(directorio):
    for name in files:
        list_files.append(os.path.join(root, name))


# In[8]:


genes = ['E1', 'E2', 'E4', 'E5', 'E6', 'E7', 'L2', 'L1']


# In[14]:


#procesar los archivos genbank para obtener la secuencia genómica completa, 
#el concatenado y los genes por separado

def process_genbank(list_files):
    
    genome_genes = defaultdict(dict)
    whole_genome = defaultdict(dict)
    genome_genes_cat = defaultdict(dict)

    for file in list_files:
        with open(file, "r") as f_hpv:
            data = SeqIO.parse(f_hpv, "genbank")
        
            virus_genus = f_hpv.name.split("/")[-2]
            virus_type = (re.findall(r'\d+', f_hpv.name.split("/")[-1])[0])
        
            whole_genome_seq = []
            genome_genes_seq = []
            genes_cat_seq = []
        
            for seqrecord in data:
            
                whole_genome_seq.append(SeqRecord(seqrecord.seq, id = seqrecord.id, 
                                                  description = seqrecord.description))
    
                genes_cat = ""
            
                gen_present = {"E1" : False, "E2" : False, "E4" : False, "E5" : False, "E6" : False, 
                               "E7" : False, "L1" : False, "L2" : False}
            
                if seqrecord.features:
                    gene=""
            
                    for feature in seqrecord.features:

                        if feature.type == "gene":
                        
                            if feature.qualifiers == "gene":
                                gene = feature.qualifiers["gene"][0]
                                gene = gene.upper()
                            
                        elif feature.type == "CDS":
                            if "gene" in feature.qualifiers:
                                if feature.qualifiers["gene"][0].split()[0] in genes:
                                    cds = feature.qualifiers["gene"][0].split()[0]
                                    cds = cds.upper()
                                
                                    if not gen_present[cds]:
                                        
                                        gen_present[cds] = True
                                    
                                        if gene and (cds != gene): #error
                                            warn = seqrecord.description + ";" + cds
                                            logging.warning('%s raised a warning, gene and CDS                                             are different', warn)
                                    #ignoraré este CDS porque el feature tipo gene y el cds.gene no coinciden
                                        else: 
                                            seq_cds = feature.extract(seqrecord.seq)
                                            genome_genes_seq.append(SeqRecord(seq_cds, 
                                                                              id = seqrecord.id, 
                                                                              name = cds, 
                                                                        description = seqrecord.description))
                                        
                                            if gene != "E5":
                                                genes_cat = genes_cat + seq_cds
                                    else:
                                        warn = seqrecord.description + ";" + cds
                                        logging.warning('%s raised a warning, gene is                                             duplicated', warn)
                                
                                        
                            #algunos genes están anotados como ORF pero se corresponden con las mismas unidades,
                            # por lo que si coinciden las tenemos en cuenta igual
                                elif "ORF" in feature.qualifiers["gene"][0].split():
                                
                                    name_orf = feature.qualifiers["gene"][0].split()
                                    name_orf.remove("ORF")
                                
                                    name_orf = name_orf[0].upper()
                                
                                    warn = seqrecord.description + ";" + name_orf
                                    logging.warning('%s raised a warning, gene name included ORF', warn)
                                
                                    if name_orf in genes:
                                        cds = name_orf
                                    
                                        if not gen_present[cds]:
                                    
                                            gen_present[cds] = True

                                            if gene and (cds != gene): #error
                                                warn = seqrecord.description + ";" + cds
                                                logging.warning('%s raised a warning, gene and CDS are different', 
                                                                warn)
                                            else: 
                                                seq_cds = feature.extract(seqrecord.seq)
                                                genome_genes_seq.append(SeqRecord(seq_cds, id = seqrecord.id, 
                                                                                  name = cds, 
                                                                            description = seqrecord.description))
                                                if gene != "E5":
                                                    genes_cat = genes_cat + seq_cds
                                        else:
                                            warn = seqrecord.description + ";" + cds
                                            logging.warning('%s raised a warning, gene is duplicated', warn)
                                
                                    else:
                                        warn = seqrecord.description + ";" + name_orf
                                        logging.warning('%s raised a warning, gene name not in genes', warn)
                            
                            # Si el nombre del gen no está en nuestra lista y tampoco presenta ORF hay un 
                            # warning y lo ignoramos 
                                else:
                                    warn = seqrecord.description + ";" + feature.qualifiers["gene"][0].split()[0]
                                    logging.warning('%s raised a warning, gene name not in genes', warn)
                       
                        #si no hay cualificador gene en CDS miro si hay product
                            elif "product" in feature.qualifiers: 
                                if feature.qualifiers["product"][0].split()[0] in genes:
                                    cds = feature.qualifiers["product"][0].split()[0]
                                    cds = cds.upper()
                                
                                    if not gen_present[cds]:
                                        
                                        gen_present[cds] = True
    
                                        if gene and (cds != gene): #error
                                            warn = seqrecord.description + ";" + cds
                                            logging.warning('%s raised a warning, gene and CDS are different', 
                                                            warn)
                                        else: 
                                            seq_cds = feature.extract(seqrecord.seq)
                                            genome_genes_seq.append(SeqRecord(seq_cds, id = seqrecord.id, 
                                                                              name = cds, 
                                                                        description = seqrecord.description))
                                        
                                            if gene != "E5":
                                                genes_cat = genes_cat + seq_cds
                                    else:
                                        warn = seqrecord.description + ";" + cds
                                        logging.warning('%s raised a warning, gene is duplicated', warn)
                                
                                                
                                elif "ORF" in feature.qualifiers["product"][0].split():
                                
                                    name_orf = feature.qualifiers["product"][0].split()
                                    name_orf.remove("ORF")
                                
                                    name_orf = name_orf[0].upper()
                                
                                    warn = seqrecord.description + ";" + name_orf
                                    logging.warning('%s raised a warning, product name included ORF', warn)
                                
                                    if name_orf in genes:
                                        cds = name_orf
                                    
                                        if not gen_present[cds]:
                                    
                                            gen_present[cds] = True

                                            if gene and (cds != gene): #error
                                                warn = seqrecord.description + ";" + cds
                                                logging.warning('%s raised a warning, gene and CDS are different', 
                                                                warn)
                                        
                                            else: 
                                                seq_cds = feature.extract(seqrecord.seq)
                                                genome_genes_seq.append(SeqRecord(seq_cds, id = seqrecord.id, 
                                                                                  name = cds, 
                                                                            description = seqrecord.description))
                                                if gene != "E5":
                                                    genes_cat = genes_cat + seq_cds
                                        else:
                                            warn = seqrecord.description + ";" + cds
                                            logging.warning('%s raised a warning, gene is duplicated', warn)
                                
                                    else:
                                        warn = seqrecord.description + ";" + name_orf
                                        logging.warning('%s raised a warning, gene name not in genes', warn)
                            
                if len(genes_cat) > 0:
                    genes_cat_seq.append(SeqRecord(genes_cat, id = seqrecord.id, 
                                                   description = seqrecord.description))
        
            whole_genome[virus_genus][virus_type] = whole_genome_seq
        
            if len(genes_cat_seq) > 0:
                genome_genes_cat[virus_genus][virus_type] = genes_cat_seq
                genome_genes[virus_genus][virus_type] = genome_genes_seq
        
        

    return whole_genome, genome_genes_cat, genome_genes


# In[15]:


whole_genome, genome_genes_cat, genome_genes = process_genbank(list_files)


# In[ ]:


#obtener valores GC 
def gc123_genomes(genomes_dict):
    
    gc123_output_list = defaultdict(dict)
    
    for virus_genus in genomes_dict.keys():
        
        for virus_type in genomes_dict[virus_genus].keys():
        
            gc123_values = []
        
            for record in genomes_dict[virus_genus][virus_type]:
                
                gc123 = list(GC123(record.seq))
                gc12 = np.mean([gc123[1],gc123[2]])
                gc123.append(gc12)
            
                gc123_values.append(gc123)

            gc123_output_list[virus_genus][virus_type] = gc123_values
    
    return gc123_output_list


# In[ ]:


def df_genomes_gc_original(gc123_list):
    
    df_final = pd.DataFrame()
    
    for virus_genus in gc123_list.keys():
        
        for virus_type in gc123_list[virus_genus].keys():
        
            df_gc = pd.DataFrame.from_dict(gc123_list[virus_genus][virus_type])
            
            col_gen = [virus_genus] * len(gc123_list[virus_genus][virus_type])
            col_type = [virus_type] * len(gc123_list[virus_genus][virus_type])

            df_gc["Genero"] = col_gen
            df_gc["Tipo"] = col_type
            
            cols = df_gc.columns.tolist()
            cols = [cols[5], cols[6], cols[0], cols[4], cols[3]]
            df_gc = df_gc[cols]
            df_gc.rename(columns = {0: "Total_GC", 4:"GC12", 3:"GC3"}, inplace = True)
            
            
            
            df_final = pd.concat([df_final, df_gc], axis = 0)
        
    return df_final


# In[11]:


#eliminar codones seleccionados de la secuencia de los genes 

def gc3s_genes(genome_genes):
    genomes_dict_by_type_and_gen_gc3s = defaultdict(dict)

    for virus_genus in genome_genes.keys():
    
        for virus_type in genome_genes[virus_genus].keys():
    
            genomes_dict_by_gen_gc3s = defaultdict(list)
        
            for record in genome_genes[virus_genus][virus_type]:
    
                sequ = Seq(str(record.seq))
                codons = [sequ[i:i+3] for i in range(0, len(sequ), 3)]
                delete = ["ATG", "TGA", "TAG", "TAA", "TGG"]
    
                sequence = ""
        
                for triplete in codons:
                    if triplete not in delete:
                        sequence = sequence + triplete

                try: 
                    len(sequence) % 3 != 0
                
                except:
                    err = virus_genus + ";" + virus_type 
                    logging.error("%s raised an error, it cannot be divided by 3")
            
            
                genomes_dict_by_gen_gc3s[record.name].append(sequence)
            
            genomes_dict_by_type_and_gen_gc3s[virus_genus][virus_type] = genomes_dict_by_gen_gc3s


# In[13]:


# obtener los valores GC en el caso de los genes

def gc_genes(genomes_dict_by_type_and_gen_gc3s):
    gc123_by_type_and_gene_list_todas = defaultdict(dict)

    for virus_genus in genomes_dict_by_type_and_gen_gc3s.keys():
    
        for virus_type in genomes_dict_by_type_and_gen_gc3s[virus_genus]:
        
            medias_by_gene = {}
            desv_by_gene = {}
        
            available_genes = genomes_dict_by_type_and_gen_gc3s[virus_genus][virus_type].keys()
        
            gc123_by_gen_dict = defaultdict(list)
            
            for gene in available_genes:
            
                for sequence in genomes_dict_by_type_and_gen_gc3s[virus_genus][virus_type][gene]:
                
                    gc123 = list(GC123(sequence))
                    gc12 = np.mean([gc123[1],gc123[2]])
                    gc123.append(gc12)
                    gc123_by_gen_dict[gene].append(gc123)
                
                
                gc123_by_gen_gc3s_array = np.asarray(gc123_by_gen_dict[gene])
            
            gc123_by_type_and_gene_list_todas[virus_genus][virus_type] = gc123_by_gen_dict
        


# In[17]:


def df_genomes_genes_gc_original(gc123_genes_list):
    
    df_final = pd.DataFrame()
    
    for virus_genus in gc123_genes_list.keys():
        
        for virus_type in gc123_genes_list[virus_genus].keys():
            
            for gene in gc123_genes_list[virus_genus][virus_type].keys():
                
                df_gc = pd.DataFrame.from_dict(gc123_genes_list[virus_genus][virus_type][gene])
                
                
                col_genero = [virus_genus] * len(gc123_genes_list[virus_genus][virus_type][gene])
                col_type = [virus_type] * len(gc123_genes_list[virus_genus][virus_type][gene])
                col_gen = [gene] * len(gc123_genes_list[virus_genus][virus_type][gene])
                
                
                df_gc["Genero"] = col_genero
                df_gc["Tipo"] = col_type
                df_gc["Gen"] = col_gen
                
                cols = df_gc.columns.tolist()
                cols = [cols[5], cols[6], cols[7], cols[0], cols[4], cols[3]]
                df_gc = df_gc[cols]
                df_gc.rename(columns = {0: "Total_GC", 4:"GC12", 3:"GC3"}, inplace = True)
            
                df_final = pd.concat([df_final, df_gc], axis = 0)
        
    return df_final


# # Dinucleótidos

# In[ ]:


#obtener la lista de dinucleótidos
nucl = ["A", "C", "T", "G"]
dinucl = [p for p in itertools.product(nucl, repeat = 2)]
dinucleotidos = []

for pair in dinucl:
    dinucleotidos.append(pair[0]+pair[1])
#print(dinucleotidos)


# In[ ]:


#determinar que nucleótidos hay que tener en cuenta para el cálculo de
#frecuencias en función de los que componen el dinucleótido

def value_dinucleotidos(pair, C_total, G_total, T_total, A_total, din_dict):

        din_dict[pair]
        
        x = {}
        y = {}
        if "A" in pair:
            x = A_total
            if pair == "AA":
                y = A_total
            
        if "T" in pair:
            if len(x.keys()) == 0:
                x = T_total
                if pair == "TT":
                    y = T_total
            else: 
                y = T_total
     
        if "C" in pair:
            if len(x.keys()) == 0:
                x = C_total
                if pair == "CC":
                    y = C_total
            else: 
                y = C_total
 
        if "G" in pair:
            if len(x.keys()) == 0:
                x = G_total
                if pair == "GG":
                    y = G_total
            else: 
                y = G_total
                
        return(x, y)


# In[ ]:


#cálculo de la frecuencia de dinucleótidos cuando no se tiene en cuenta la posición
def freq_dinucleotidos(genome):
    
    C_total = defaultdict(dict)
    G_total = defaultdict(dict)

    T_total = defaultdict(dict)
    A_total = defaultdict(dict)

    dinucleotidos_por_tipo = defaultdict(dict)

    for virus_genus in genome.keys():
    
        for virus_type in genome[virus_genus].keys():
            
            din_dict = defaultdict(list)
        
            for pair in dinucleotidos:
                din_dict[pair]
        
            C_total_list = [] 
            G_total_list = [] 
            T_total_list = [] 
            A_total_list = [] 
        
            for record in genome[virus_genus][virus_type]:
            
                sequence = Seq(str(record.seq))
                counter = Counter(sequence)
    
    
                C_total_list.append((counter["C"])/len(sequence))
                G_total_list.append((counter["G"])/len(sequence))
                T_total_list.append((counter["T"])/len(sequence))
                A_total_list.append((counter["A"])/len(sequence))
    
        
    
                for pair in dinucleotidos:
            
                    frequency = freq_overlapping(sequence, pair)

                    din_dict[pair].append(frequency)
        
            C_total[virus_genus][virus_type] = C_total_list
            G_total[virus_genus][virus_type] = G_total_list
            T_total[virus_genus][virus_type] = T_total_list
            A_total[virus_genus][virus_type] = A_total_list
        
        
            dinucleotidos_por_tipo[virus_genus][virus_type] = din_dict
            
    return dinucleotidos_por_tipo, C_total, G_total, T_total, A_total



#cálculo de la frecuencia de dicucleótidos teniendo en cuenta la posición
def freq_dinucleotidos_123_genome(genome):
    
    dinucleotidos_12 = defaultdict(dict)
    dinucleotidos_23 = defaultdict(dict)
    
    for virus_genus in genome.keys():
                
        for virus_type in genome[virus_genus].keys():
            
            din_dict_12 = defaultdict(list)
            din_dict_23 = defaultdict(list)

            for sequence in genome[virus_genus][virus_type]:
                
                sequ = Seq(str(sequence.seq))
                
                codons = [sequ[i:i+3] for i in range(0, len(sequ), 3)]
                codons_12 = [codons[i][0:2] for i in range (len(codons))]
                codons_23 = [codons[i][1:3] for i in range (len(codons))]
    
                count_12 = Counter(codons_12)
                count_23 = Counter(codons_23)
    
                for pair in dinucleotidos:
            
                    freq_12 = count_12[pair]/len(codons)
                    din_dict_12[pair].append(freq_12)
                    
                    freq_23 = count_23[pair]/len(codons)
                    din_dict_23[pair].append(freq_23)
            
            dinucleotidos_12[virus_genus][virus_type] = din_dict_12
            dinucleotidos_23[virus_genus][virus_type] = din_dict_23
            
    
    return dinucleotidos_12, dinucleotidos_23


#cálculo de la frecuencia observados/esperados
def obs_esp_genoma(frequency_dinucleotidos, C_total, G_total, T_total, A_total):

    
    obs_esp_dinucl_por_tipo = defaultdict(dict)

    for virus_genus in frequency_dinucleotidos.keys():
    
        for virus_type in frequency_dinucleotidos[virus_genus].keys():
        
            din_dict = defaultdict(list)
        
            for pair in dinucleotidos:
            
                x, y = value_dinucleotidos(pair, C_total, G_total, T_total, A_total, din_dict)
              
                for value in range(len(frequency_dinucleotidos[virus_genus][virus_type][pair])):

                    obs_esp = frequency_dinucleotidos[virus_genus][virus_type][pair][value] /                                             (x[virus_genus][virus_type][value] * y[virus_genus][virus_type][value])
                    din_dict[pair].append(obs_esp)
            
            obs_esp_dinucl_por_tipo[virus_genus][virus_type] = din_dict
    return obs_esp_dinucl_por_tipo


# In[ ]:


def df_genomes_dinucleotidos_original(obs_esp):
    
    df_final = pd.DataFrame()
    
    for virus_genus in obs_esp.keys():
        
        for virus_type in obs_esp[virus_genus].keys():
            
            df_din = pd.DataFrame()
            
            for pair in obs_esp[virus_genus][virus_type].keys():
            
                df_din_pair = pd.DataFrame.from_dict(obs_esp[virus_genus][virus_type][pair])
            
                df_din_pair.rename(columns = {0: pair}, inplace = True)
                #print(df_din_pair)
                
                df_din = pd.concat([df_din, df_din_pair], axis = 1)
                
            col_gen = [virus_genus] * len(df_din)
            col_type = [virus_type] * len(df_din)
            
            df_din["Genero"] = col_gen
            df_din["Tipo"] = col_type
                
            df_final = pd.concat([df_final, df_din], axis = 0)
        
    return df_final


# In[ ]:


#cálculo de la frecuencia de los dinucleótidos en los genes, 
#sin tener en cuenta la posición

def freq_dinucleotidos_genes(genome_by_genes):
    
    C_total = defaultdict(dict)
    G_total = defaultdict(dict)

    T_total = defaultdict(dict)
    A_total = defaultdict(dict)

    dinucleotidos_por_tipo = defaultdict(dict)

    for virus_genus in genome_by_genes.keys():
    
        for virus_type in genome_by_genes[virus_genus].keys():
            
            din_dict_genes = defaultdict(dict)
            C_total_list = defaultdict(list) 
            G_total_list = defaultdict(list)
            T_total_list = defaultdict(list)
            A_total_list = defaultdict(list) 
            
            for gene in genome_by_genes[virus_genus][virus_type].keys():
                
                din_dict = defaultdict(list)
        
                for pair in dinucleotidos:
                    din_dict[pair]
        
                
                for record in genome_by_genes[virus_genus][virus_type][gene]:
                    
                    counter = Counter(record)
    
        
                    C_total_list[gene].append((counter["C"])/len(record))
                    G_total_list[gene].append((counter["G"])/len(record))
                    T_total_list[gene].append((counter["T"])/len(record))
                    A_total_list[gene].append((counter["A"])/len(record))
    
                    for pair in dinucleotidos:
            
                        frequency = freq_overlapping(record, pair)

                        din_dict[pair].append(frequency)
        
                
                din_dict_genes[gene] = din_dict
            C_total[virus_genus][virus_type] = C_total_list
            G_total[virus_genus][virus_type] = G_total_list
            T_total[virus_genus][virus_type] = T_total_list
            A_total[virus_genus][virus_type] = A_total_list
        
            dinucleotidos_por_tipo[virus_genus][virus_type] = din_dict_genes
            
    return dinucleotidos_por_tipo, C_total, G_total, T_total, A_total


# In[ ]:


#cálculo de la frecuencia de los dinucleótidos en los genes,
#teniendo en cuenta la posición

def freq_dinucleotidos_123_genes(genome):
    
    dinucleotidos_12 = defaultdict(dict)
    dinucleotidos_23 = defaultdict(dict)
    
    for virus_genus in genome.keys():
                
        for virus_type in genome[virus_genus].keys():
            
            din_dict_genes_12 = defaultdict(dict)
            din_dict_genes_23 = defaultdict(dict)
            
            for gen in genome[virus_genus][virus_type].keys():
                
                din_dict_12 = defaultdict(list)
                din_dict_23 = defaultdict(list)
            
                for sequence in genome[virus_genus][virus_type][gen]:
    
                    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
                    codons_12 = [codons[i][0:2] for i in range (len(codons))]
                    codons_23 = [codons[i][1:3] for i in range (len(codons))]
    
                    count_12 = Counter(codons_12)
                    count_23 = Counter(codons_23)
    
                    for pair in dinucleotidos:
                        freq_12 = count_12[pair]/len(codons)
                        din_dict_12[pair].append(freq_12)
                        freq_23 = count_23[pair]/len(codons)
                        din_dict_23[pair].append(freq_23)
                
                din_dict_genes_12[gen] = din_dict_12
                din_dict_genes_23[gen] = din_dict_23
                
            dinucleotidos_12[virus_genus][virus_type] = din_dict_genes_12
            dinucleotidos_23[virus_genus][virus_type] = din_dict_genes_23
    
    return dinucleotidos_12, dinucleotidos_23


# In[ ]:


#cálculo de la frecuencia observados/esperados en los genes

def obs_esp_genes(frequency_dinucleotidos, C_total, G_total, T_total, A_total):

    
    obs_esp_dinucl_por_tipo = defaultdict(dict)

    for virus_genus in frequency_dinucleotidos.keys():
    
        for virus_type in frequency_dinucleotidos[virus_genus].keys():
            
            din_dict_genes = defaultdict(dict)
        
            for gene in frequency_dinucleotidos[virus_genus][virus_type].keys():
                
                din_dict = defaultdict(list)
        
                for pair in dinucleotidos:
            
                    x, y = value_dinucleotidos(pair, C_total, G_total, T_total, A_total, din_dict)
                    
                    for value in range(len(frequency_dinucleotidos[virus_genus][virus_type][gene][pair])):

                        obs_esp = frequency_dinucleotidos[virus_genus][virus_type][gene][pair][value] /                                                 (x[virus_genus][virus_type][gene][value] * 
                                                 y[virus_genus][virus_type][gene][value])
                        din_dict[pair].append(obs_esp)
                
                din_dict_genes[gene] = din_dict
            
            obs_esp_dinucl_por_tipo[virus_genus][virus_type] = din_dict_genes
                
    return obs_esp_dinucl_por_tipo


# In[ ]:


def df_genomes_genes_gc_original(obs_esp_gen):
    
    df_final = pd.DataFrame()
    
    for virus_genus in obs_esp_gen.keys():
        
        for virus_type in obs_esp_gen[virus_genus].keys():
            
            for gene in obs_esp_gen[virus_genus][virus_type].keys():
                
                df_din = pd.DataFrame()
                
                for pair in obs_esp_gen[virus_genus][virus_type][gene].keys():
                    
                    df_din_pair = pd.DataFrame.from_dict(obs_esp_gen[virus_genus][virus_type][gene])
            
                    df_din_pair.rename(columns = {0: pair}, inplace = True)
                
                    df_din = pd.concat([df_din, df_din_pair], axis = 1)
                
                col_genero = [virus_genus] * len(df_din)
                col_type = [virus_type] * len(df_din)
                col_gen = [gene] * len(df_din)
                
                df_din["Genero"] = col_genero
                df_din["Tipo"] = col_type
                df_din["Gen"] = col_gen
                
                df_final = pd.concat([df_final, df_din], axis = 0)
        
    return df_final

