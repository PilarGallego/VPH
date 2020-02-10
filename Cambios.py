#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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
from matplotlib import cm #https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.libqsturng import psturng
from tabulate import tabulate

get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


logging.basicConfig(filename = "avisos_cambios.log", filemode = "w", 
                    format = "%(name)s - %(levelname)s - %(asctime)s - %(message)s",
                    level = logging.INFO)


# In[3]:


list_files = []
for root, directories, files in os.walk("./GenomasHPV/"):
    for name in files:
        list_files.append(os.path.join(root, name))


# In[4]:


genes = ['E1', 'E2', 'E4', 'E5', 'E6', 'E7', 'L2', 'L1']


# In[5]:


genome_genes = defaultdict(dict)
whole_genome = defaultdict(dict)
genome_genes_cat = defaultdict(dict)

for file in list_files:
    with open(file, "r") as f_hpv:
        data = SeqIO.parse(f_hpv, "genbank")
        
        virus_genus = f_hpv.name.split("/")[-2]
        virus_type = (re.findall(r'\d+', f_hpv.name.split("/")[-1])[0])
        
        
        #print(virus_genus, ";", virus_type)
        
        whole_genome_seq = []
        genome_genes_seq = []
        genes_cat_seq = []
        
        for seqrecord in data:
            
            whole_genome_seq.append(SeqRecord(seqrecord.seq, id = seqrecord.id, description = seqrecord.description))

            genes_cat = ""
            
            gen_present = {"E1" : False, "E2" : False, "E4" : False, "E5" : False, "E6" : False, "E7" : False,
                          "L1" : False, "L2" : False}
            
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
                                        logging.warning('%s raised a warning, gene and CDS are different', warn)
                                    #ignoraré este CDS porque el feature tipo gene y el cds.gene no coinciden
                                    else: 
                                        seq_cds = feature.extract(seqrecord.seq)
                                        genome_genes_seq.append(SeqRecord(seq_cds, id = seqrecord.id, name = cds, 
                                                              description = seqrecord.description))
                                        
                                        if gene != "E5":
                                            genes_cat = genes_cat + seq_cds
                                else:
                                    warn = seqrecord.description + ";" + cds
                                    logging.warning('%s raised a warning, gene is duplicated', warn)
                                
                                        
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
                                            logging.warning('%s raised a warning, gene and CDS are different', warn)
                                        else: 
                                            seq_cds = feature.extract(seqrecord.seq)
                                            genome_genes_seq.append(SeqRecord(seq_cds, id = seqrecord.id, name = cds, 
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
                                        logging.warning('%s raised a warning, gene and CDS are different', warn)
                                    else: 
                                        seq_cds = feature.extract(seqrecord.seq)
                                        genome_genes_seq.append(SeqRecord(seq_cds, id = seqrecord.id, name = cds, 
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
                                            logging.warning('%s raised a warning, gene and CDS are different', warn)
                                        
                                        else: 
                                            seq_cds = feature.extract(seqrecord.seq)
                                            genome_genes_seq.append(SeqRecord(seq_cds, id = seqrecord.id, name = cds, 
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
                genes_cat_seq.append(SeqRecord(genes_cat, id = seqrecord.id, description = seqrecord.description))
        
        whole_genome[virus_genus][virus_type] = whole_genome_seq
        
        if len(genes_cat_seq) > 0:
            genome_genes_cat[virus_genus][virus_type] = genes_cat_seq
            genome_genes[virus_genus][virus_type] = genome_genes_seq
        


# In[58]:


for virus_type in genome_genes_cat["Alfa"]:
    
    print(virus_type + " ; " + str(len(genome_genes_cat["Alfa"][virus_type])))


# In[10]:


for virus_type in genome_genes_cat["Alfa"]:
    
    print(virus_type + " ; " + str(len(genome_genes_cat["Alfa"][virus_type])))


# In[107]:


for virus_type in whole_genome["Alfa"]:
    
    print(virus_type + " ; " + str(len(whole_genome["Alfa"][virus_type])))


# tipo 67 y 58

# In[39]:


for virus_type in whole_genome["Beta"]:
    
    print(virus_type + " ; " + str(len(whole_genome["Beta"][virus_type])))


# In[40]:


for virus_type in whole_genome["Gamma"]:
    
    print(virus_type + " ; " + str(len(whole_genome["Gamma"][virus_type])))


# In[6]:


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


# In[7]:


gc123_by_type_cat_list = gc123_genomes(genome_genes_cat)


# In[8]:


gc123_by_type_whole_list = gc123_genomes(whole_genome)


# In[76]:


gc123_by_type_whole_list


# In[28]:


genome_genes_cat["Alfa"]["68"][0]


# In[9]:


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


# In[10]:


df_gc_cat_original = df_genomes_gc_original(gc123_by_type_cat_list)


# In[11]:


df_gc_whole_original = df_genomes_gc_original(gc123_by_type_whole_list)


# In[85]:


df_gc_cat_original


# In[14]:


df_gc_cat_original.to_csv("./results/genome_cat_gc_original_cambios.csv", sep  = "\t", index = False, header = True)


# In[12]:


df_gc_whole_original.to_csv("./results/whole_genome_gc_original_cambios.csv", sep  = "\t", index = False, header = True)


# In[15]:


cat_gc_alfa_beta_gamma = df_gc_cat_original.loc[df_gc_cat_original["Genero"].isin(["Alfa", "Beta", "Gamma"])]


# In[16]:


whole_gc_alfa_beta_gamma = df_gc_whole_original.loc[df_gc_whole_original["Genero"].isin(["Alfa", "Beta", "Gamma"])]


# In[17]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.3)

number = 1

plt.suptitle("Genoma concatenado", fontsize = 35)

for virus_genus in cat_gc_alfa_beta_gamma["Genero"].unique():
    
    gc3 = cat_gc_alfa_beta_gamma.loc[cat_gc_alfa_beta_gamma["Genero"].isin([virus_genus])]["GC3"]
    gc12 = cat_gc_alfa_beta_gamma.loc[cat_gc_alfa_beta_gamma["Genero"].isin([virus_genus])]["GC12"]
        
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(gc3, gc12)

    plt.subplot(1,3, number)
    plt.plot(gc3, gc12, "o")
    plt.plot(gc3, intercept + slope*gc3, "r", label = "Fitted line.\nb =" + str(round(slope, 2)) 
                 +", p = "+ str('{:0.3e}'.format(pvalue)))
    
    plt.xlabel("GC\u2083", fontsize = 25)
    plt.ylabel("GC\u2081"+"\u2082", fontsize = 25)
    #plt.xlim([21, 56])
    #plt.ylim([37, 52])
    plt.xlim([13, 70])
    plt.ylim([27, 65])
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(fontsize = 20)
    plt.title("Género "+virus_genus, fontsize = 30)
    
    number += 1
    
plt.savefig("./results/neutralidad_genero_genoma_concat_cambios.png")
plt.show()


# In[18]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.3)

number = 1

plt.suptitle("Genoma entero", fontsize = 35)

for virus_genus in whole_gc_alfa_beta_gamma["Genero"].unique():
    
    gc3 = whole_gc_alfa_beta_gamma.loc[whole_gc_alfa_beta_gamma["Genero"].isin([virus_genus])]["GC3"]
    gc12 = whole_gc_alfa_beta_gamma.loc[whole_gc_alfa_beta_gamma["Genero"].isin([virus_genus])]["GC12"]
        
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(gc3, gc12)

    plt.subplot(1,3, number)
    plt.plot(gc3, gc12, "o")
    plt.plot(gc3, intercept + slope*gc3, "r", label = "Fitted line.\nb =" + str(round(slope, 2)) 
                 +", p = "+ str('{:0.3e}'.format(pvalue)))
    
    plt.xlabel("GC\u2083", fontsize = 25)
    plt.ylabel("GC\u2081"+"\u2082", fontsize = 25)
    #plt.xlim([21, 56])
    #plt.ylim([37, 52])
    plt.xlim([13, 70])
    plt.ylim([27, 65])
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(fontsize = 20)
    plt.title("Género "+virus_genus, fontsize = 30)
    
    number += 1
    
plt.show()


# In[25]:


len_whole = defaultdict(list)

for virus_type in whole_genome["Alfa"]:
    
        for record in whole_genome["Alfa"][virus_type]:
            
            len_whole[virus_type].append(len(record.seq))


# In[26]:


len_whole


# In[28]:


len_cat = defaultdict(list)

for virus_type in genome_genes_cat["Alfa"]:
    
        for record in genome_genes_cat["Alfa"][virus_type]:
            
            len_cat[virus_type].append(len(record.seq))


# In[30]:


len_cat


# In[33]:


for virus_type in len_whole:
        
    print("TIPO:" + virus_type)
    
    for value in range(len(len_cat[virus_type])):
            
        print(len_whole[virus_type][value] - len_cat[virus_type][value])
    
    if len(len_whole[virus_type]) > len(len_cat[virus_type]):
        
        print("Genoma que falta")
            
        
    print("*****************")


# In[35]:


genome_genes["Alfa"]["43"]


# ## concatenado con los 7 genes principales (es decir, el concatenado que tenía al principio)
# 

# In[36]:


len_whole = defaultdict(list)

for virus_type in whole_genome["Gamma"]:
    
        for record in whole_genome["Gamma"][virus_type]:
            
            len_whole[virus_type].append(len(record.seq))
len_cat = defaultdict(list)

for virus_type in genome_genes_cat["Gamma"]:
    
        for record in genome_genes_cat["Gamma"][virus_type]:
            
            len_cat[virus_type].append(len(record.seq))
            
for virus_type in len_whole:
        
    print("TIPO:" + virus_type)
    
    for value in range(len(len_cat[virus_type])):
            
        print(len_whole[virus_type][value] - len_cat[virus_type][value])
    
    if len(len_whole[virus_type]) > len(len_cat[virus_type]):
        
        print("Genoma que falta")
            
        
    print("*****************")


# In[37]:


len_whole = defaultdict(list)

for virus_type in whole_genome["Beta"]:
    
        for record in whole_genome["Beta"][virus_type]:
            
            len_whole[virus_type].append(len(record.seq))
len_cat = defaultdict(list)

for virus_type in genome_genes_cat["Beta"]:
    
        for record in genome_genes_cat["Beta"][virus_type]:
            
            len_cat[virus_type].append(len(record.seq))
            
for virus_type in len_whole:
        
    print("TIPO:" + virus_type)
    
    for value in range(len(len_cat[virus_type])):
            
        print(len_whole[virus_type][value] - len_cat[virus_type][value])
    
    if len(len_whole[virus_type]) > len(len_cat[virus_type]):
        
        print("Genoma que falta")
            
        
    print("*****************")


# ### poner un diccionario con los 8 genes a falso, if not dict[gene] lo pasas a True y lo añades, si no ya pasas de él y  avisas de duplicado --> ***Hecho***
# ## dict = {val: False for val in dict}
# ### cuando no sale en el diccionario de genes, en vez de ponerle solo el else ponerle un elif de que si ORF está en nombre_gen incluirlo --> ***Hecho***
# 
# ##### si da mucho por culo se pasa de lo de los orfs y se incluyen solo los que están en esa lista de genes
# 
# ### Coger solo la primera secuencia de cada tipo para los 3 géneros y hacer todos los cálculos solo con esa para ver si hay una diferencia en los resultados (sobre todo en la neutralidad) --> se puede simplemente modificar el primer paso de extraer las secuencias (hacer un duplicado y ya está, modificando que solo me coja el primer seqrecord)
# 
# 
# ### hacer una versión en la que sí que se encuentre e5 en el concatenado
# 
# 
# #### si todo esto sigue dando diferente se podría ver el gc de las zonas intergénicas y discutir pero probablemente no sea necesario llegar a ese punto
# 
# 
# ### lo del riesgo y tal también repetirlo con un solo genoma

# # ***1 solo genoma por tipo***

# In[13]:


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


# In[17]:


def one_seq(genome_dict):
    
    genome_dict_1_seq = defaultdict(dict)
    
    for virus_genus in genome_dict:
        
         for virus_type in genome_dict[virus_genus]:
                
                genome_dict_1_seq[virus_genus][virus_type] = genome_dict[virus_genus][virus_type][0]
                
    return genome_dict_1_seq


# In[18]:


genome_genes_cat_one_seq = one_seq(genome_genes_cat)
whole_genome_one_seq = one_seq(whole_genome)


# In[19]:


def one_seq_genes(genes_dict):
    
    genes_dict_1_seq = defaultdict(dict)
    
    for virus_genus in genes_dict:
        
         for virus_type in genes_dict[virus_genus]:
                
                genes_list = {}
                
                for gene in genes_dict[virus_genus][virus_type]:
                        
                        genes_list[gene] = genes_dict[virus_genus][virus_type][gene][0]
                    
                genes_dict_1_seq[virus_genus][virus_type] = genes_list
                
    return genes_dict_1_seq


# In[20]:


genomes_dict_by_type_and_gen_gc3s_one_seq = one_seq_genes(genomes_dict_by_type_and_gen_gc3s)


# In[208]:


def gc123_genomes_one_seq(genomes_dict):
    
    gc123_output_list = defaultdict(dict)
    
    for virus_genus in genomes_dict.keys():
        
        
        for virus_type in genomes_dict[virus_genus].keys():
        
            gc123 = list(GC123(genomes_dict[virus_genus][virus_type].seq))
            gc12 = np.mean([gc123[1],gc123[2]])
            gc123.append(gc12)
            
            gc123_output_list[virus_genus][virus_type] = gc123
    
    return gc123_output_list


# In[19]:


gc123_by_type_cat_list_one_seq = gc123_genomes_one_seq(genome_genes_cat_one_seq)
gc123_by_type_whole_list_one_seq = gc123_genomes_one_seq(whole_genome_one_seq)


# In[209]:


gc123_by_type_and_gene_list_one_seq = defaultdict(dict)

for virus_genus in genomes_dict_by_type_and_gen_gc3s_one_seq.keys():
    
    for virus_type in genomes_dict_by_type_and_gen_gc3s_one_seq[virus_genus]:
        
        medias_by_gene = {}
        desv_by_gene = {}
        
        available_genes = genomes_dict_by_type_and_gen_gc3s_one_seq[virus_genus][virus_type].keys()
        
        gc123_by_gen_dict = defaultdict(list)
            
        for gene in available_genes:
            
            sequence = genomes_dict_by_type_and_gen_gc3s_one_seq[virus_genus][virus_type][gene]
                
            gc123 = list(GC123(sequence))
            gc12 = np.mean([gc123[1],gc123[2]])
            gc123.append(gc12)
            gc123_by_gen_dict[gene].append(gc123)
                
            
            gc123_by_gen_gc3s_array = np.asarray(gc123_by_gen_dict[gene])
            
        gc123_by_type_and_gene_list_one_seq[virus_genus][virus_type] = gc123_by_gen_dict
        


# #### Pequeño inciso para obtener todas las secuencias de los genes con los cambios por si acaso resultaran diferentes 

# In[9]:


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
        


# In[10]:


gc123_by_type_and_gene_list_todas


# In[21]:


def df_genomes_gc_original(gc123_list):
    
    df_final = pd.DataFrame()
    
    for virus_genus in gc123_list.keys():
        
        for virus_type in gc123_list[virus_genus].keys():
        
            df_gc = pd.DataFrame.from_dict(gc123_list[virus_genus][virus_type])
            
            df_gc = df_gc.transpose()
            col_gen = [virus_genus] 
            col_type = [virus_type] 
            df_gc["Genero"] = col_gen
            df_gc["Tipo"] = col_type
            
            cols = df_gc.columns.tolist()
            cols = [cols[5], cols[6], cols[0], cols[4], cols[3]]
            df_gc = df_gc[cols]
            df_gc.rename(columns = {0: "Total_GC", 4:"GC12", 3:"GC3"}, inplace = True)
            
            
            
            df_final = pd.concat([df_final, df_gc], axis = 0)
        
    return df_final


# In[22]:


df_gc_cat_original_one_seq = df_genomes_gc_original(gc123_by_type_cat_list_one_seq)


# In[23]:


df_gc_whole_original_one_seq = df_genomes_gc_original(gc123_by_type_whole_list_one_seq)


# In[77]:


df_gc_cat_original_one_seq


# In[80]:


df_gc_cat_original_one_seq.to_csv("./results/genome_cat_gc_original_one_seq.csv", sep  = "\t", index = False, header = True)
df_gc_whole_original_one_seq.to_csv("./results/whole_genome_gc_original_one_seq.csv", sep  = "\t", index = False, header = True)


# In[13]:


def df_genomes_genes_gc_original(gc123_genes_list):
    
    df_final = pd.DataFrame()
    
    for virus_genus in gc123_genes_list.keys():
        
        for virus_type in gc123_genes_list[virus_genus].keys():
            
            for gene in gc123_genes_list[virus_genus][virus_type].keys():
                
                df_gc = pd.DataFrame.from_dict(gc123_genes_list[virus_genus][virus_type][gene])
                
                #para obtener el df de todas se quitan los #
                
                col_genero = [virus_genus] #* len(gc123_genes_list[virus_genus][virus_type][gene])
                col_type = [virus_type] #* len(gc123_genes_list[virus_genus][virus_type][gene])
                col_gen = [gene] #* len(gc123_genes_list[virus_genus][virus_type][gene])
                
                
                df_gc["Genero"] = col_genero
                df_gc["Tipo"] = col_type
                df_gc["Gen"] = col_gen
                
                cols = df_gc.columns.tolist()
                cols = [cols[5], cols[6], cols[7], cols[0], cols[4], cols[3]]
                df_gc = df_gc[cols]
                df_gc.rename(columns = {0: "Total_GC", 4:"GC12", 3:"GC3"}, inplace = True)
            
                df_final = pd.concat([df_final, df_gc], axis = 0)
        
    return df_final


# In[211]:


df_gc_gen_original_one_seq = df_genomes_genes_gc_original(gc123_by_type_and_gene_list_one_seq)


# In[66]:


#df_gc_gen_original_one_seq.to_csv("./results/genome_genes_gc_original_one_seq.csv", sep  = "\t", index = False, header = True)
df_gc_gen_original_one_seq = pd.read_csv("./results/genome_genes_gc_original_one_seq.csv", sep  = "\t")


# In[14]:


df_gc_gen_original_todas = df_genomes_genes_gc_original(gc123_by_type_and_gene_list_todas)


# In[15]:


df_gc_gen_original_todas


# In[16]:


#df_gc_gen_original_todas.to_csv("./results/genome_genes_gc_original_todas.csv", sep  = "\t", index = False, header = True)


# In[4]:


df_gc_gen_original_one_seq


# In[26]:


cat_gc_alfa_beta_gamma_one_seq = df_gc_cat_original_one_seq.loc[df_gc_cat_original_one_seq["Genero"].isin(["Alfa", "Beta", "Gamma"])]
whole_gc_alfa_beta_gamma_one_seq = df_gc_whole_original_one_seq.loc[df_gc_whole_original_one_seq["Genero"].isin(["Alfa", "Beta", "Gamma"])]


# In[67]:


genes_gc_alfa_one_seq = df_gc_gen_original_one_seq.loc[df_gc_gen_original_one_seq["Genero"].isin(["Alfa"])]
genes_gc_beta_one_seq = df_gc_gen_original_one_seq.loc[df_gc_gen_original_one_seq["Genero"].isin(["Beta"])]
genes_gc_gamma_one_seq = df_gc_gen_original_one_seq.loc[df_gc_gen_original_one_seq["Genero"].isin(["Gamma"])]


# In[28]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.3)

number = 1

plt.suptitle("Genoma concatenado. 1 secuencia", fontsize = 35)

for virus_genus in cat_gc_alfa_beta_gamma_one_seq["Genero"].unique():
    
    gc3 = cat_gc_alfa_beta_gamma_one_seq.loc[cat_gc_alfa_beta_gamma_one_seq["Genero"].isin([virus_genus])]["GC3"]
    gc12 = cat_gc_alfa_beta_gamma_one_seq.loc[cat_gc_alfa_beta_gamma_one_seq["Genero"].isin([virus_genus])]["GC12"]
        
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(gc3, gc12)

    plt.subplot(1,3, number)
    plt.plot(gc3, gc12, "o")
    plt.plot(gc3, intercept + slope*gc3, "r", label = "Fitted line.\nb =" + str(round(slope, 2)) 
                 +", p = "+ str('{:0.3e}'.format(pvalue)))
    
    plt.xlabel("GC\u2083", fontsize = 25)
    plt.ylabel("GC\u2081"+"\u2082", fontsize = 25)
    plt.xlim([13, 70])
    plt.ylim([27, 65])
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(fontsize = 20)
    plt.title("Género "+virus_genus, fontsize = 30)
    
    number += 1
    
plt.savefig("./results/neutralidad_genero_genoma_concat_one_seq.png")
plt.show()


# In[29]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.3)

number = 1

plt.suptitle("Genoma entero. 1 secuencia", fontsize = 35)

for virus_genus in whole_gc_alfa_beta_gamma_one_seq["Genero"].unique():
    
    gc3 = whole_gc_alfa_beta_gamma_one_seq.loc[whole_gc_alfa_beta_gamma_one_seq["Genero"].isin([virus_genus])]["GC3"]
    gc12 = whole_gc_alfa_beta_gamma_one_seq.loc[whole_gc_alfa_beta_gamma_one_seq["Genero"].isin([virus_genus])]["GC12"]
        
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(gc3, gc12)

    plt.subplot(1,3, number)
    plt.plot(gc3, gc12, "o")
    plt.plot(gc3, intercept + slope*gc3, "r", label = "Fitted line.\nb =" + str(round(slope, 2)) 
                 +", p = "+ str('{:0.3e}'.format(pvalue)))
    
    plt.xlabel("GC\u2083", fontsize = 25)
    plt.ylabel("GC\u2081"+"\u2082", fontsize = 25)
    #plt.xlim([21, 56])
    #plt.ylim([37, 52])
    plt.xlim([13, 70])
    plt.ylim([27, 65])
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(fontsize = 20)
    plt.title("Género "+virus_genus, fontsize = 30)
    
    number += 1
    
plt.savefig("./results/neutralidad_genero_genoma_entero_one_seq.png")

plt.show()


# In[30]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 35
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.3, top = 0.93)

number = 1

plt.suptitle("Género Alfa. 1 secuencia", fontsize = 35)

for gen in genes:
    
    gc3 = genes_gc_alfa_one_seq.loc[genes_gc_alfa_one_seq["Gen"].isin([gen])]["GC3"]
    gc12 = genes_gc_alfa_one_seq.loc[genes_gc_alfa_one_seq["Gen"].isin([gen])]["GC12"]
        
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(gc3, gc12)

    plt.subplot(3,3, number)
    plt.plot(gc3, gc12, "o")
    plt.plot(gc3, intercept + slope*gc3, "r", label = "Fitted line.\nb =" + str(round(slope, 2)) 
                 +", p = "+ str('{:0.3e}'.format(pvalue)))
    
    plt.xlabel("GC\u2083", fontsize = 25)
    plt.ylabel("GC\u2081"+"\u2082", fontsize = 25)
    plt.xlim([13, 70])
    plt.ylim([27, 65])
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(fontsize = 20)
    plt.title("Gen "+gen, fontsize = 30)
    
    number += 1
    
plt.savefig("./results/neutralidad_genes_alfa_one_seq.png")
plt.show()


# In[69]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 35
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.3, top = 0.93)

number = 1

plt.suptitle("Género Beta. 1 secuencia", fontsize = 35)

for gen in genes:
    
    if gen in genes_gc_beta_one_seq["Gen"].unique() and gen != "E5":
        
        gc3 = genes_gc_beta_one_seq.loc[genes_gc_beta_one_seq["Gen"].isin([gen])]["GC3"]
        gc12 = genes_gc_beta_one_seq.loc[genes_gc_beta_one_seq["Gen"].isin([gen])]["GC12"]
        
        slope, intercept, rvalue, pvalue, stderr = stats.linregress(gc3, gc12)

        plt.subplot(3,3, number)
        plt.plot(gc3, gc12, "o")
        plt.plot(gc3, intercept + slope*gc3, "r", label = "Fitted line.\nb =" + str(round(slope, 2)) 
                     +", p = "+ str('{:0.3e}'.format(pvalue)))
    
        plt.xlabel("GC\u2083", fontsize = 25)
        plt.ylabel("GC\u2081"+"\u2082", fontsize = 25)
        plt.xlim([16, 88])
        plt.ylim([34, 69])
        plt.xticks(fontsize = 15)
        plt.yticks(fontsize = 15)
        plt.legend(fontsize = 20)
        plt.title("Gen "+gen, fontsize = 30)
        
        number += 1
    
plt.savefig("./results/neutralidad_genes_beta_one_seq.png")
plt.show()


# In[18]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 35
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.3, top = 0.93)

number = 1

plt.suptitle("Género Gamma. 1 secuencia", fontsize = 35)

for gen in genes:
    
    if gen in genes_gc_gamma_one_seq["Gen"].unique():
        
        gc3 = genes_gc_gamma_one_seq.loc[genes_gc_gamma_one_seq["Gen"].isin([gen])]["GC3"]
        gc12 = genes_gc_gamma_one_seq.loc[genes_gc_gamma_one_seq["Gen"].isin([gen])]["GC12"]
        
        slope, intercept, rvalue, pvalue, stderr = stats.linregress(gc3, gc12)

        plt.subplot(3,3, number)
        plt.plot(gc3, gc12, "o")
        plt.plot(gc3, intercept + slope*gc3, "r", label = "Fitted line.\nb =" + str(round(slope, 2)) 
                     +", p = "+ str('{:0.3e}'.format(pvalue)))
    
        plt.xlabel("GC\u2083", fontsize = 25)
        plt.ylabel("GC\u2081"+"\u2082", fontsize = 25)
        plt.xlim([13, 65])
        plt.ylim([29, 58])
        plt.xticks(fontsize = 15)
        plt.yticks(fontsize = 15)
        plt.legend(fontsize = 20)
        plt.title("Gen "+gen, fontsize = 30)
        
        number += 1
    
plt.savefig("./results/neutralidad_genes_gamma_one_seq.png")
plt.show()


# In[33]:


def anovas(data_df):
    mod_anova_total = ols("Total_GC ~ Genero", data = data_df).fit()
    mod_anova_12 = ols("GC12 ~ Genero", data = data_df).fit()
    mod_anova_3 = ols("GC3 ~ Genero", data = data_df).fit()
    
    anova_table_total = anova.anova_lm(mod_anova_total, typ = 2)
    anova_table_12 = anova.anova_lm(mod_anova_12, typ = 2)
    anova_table_3 = anova.anova_lm(mod_anova_3, typ = 2)
    
    return anova_table_total, anova_table_12, anova_table_3


# In[34]:


anova_whole_gc_total_one_seq, anova_whole_gc_12_one_seq, anova_whole_gc_3_one_seq = anovas(whole_gc_alfa_beta_gamma_one_seq)
anova_cat_gc_total_one_seq, anova_cat_gc_12_one_seq, anova_cat_gc_3_one_seq = anovas(cat_gc_alfa_beta_gamma_one_seq)


# In[35]:


def post_hoc(data_df):
    result_total = (pairwise_tukeyhsd(data_df["Total_GC"], data_df["Genero"], alpha = 0.05))
    result_12 = (pairwise_tukeyhsd(data_df["GC12"], data_df["Genero"], alpha = 0.05))
    result_3 = (pairwise_tukeyhsd(data_df["GC3"], data_df["Genero"], alpha = 0.05))
    
    return result_total, result_12, result_3


# In[36]:


def post_hoc_results_table(tukey_res): 
    p_value = psturng(np.abs(tukey_res.meandiffs / tukey_res.std_pairs),
           len(tukey_res.groupsunique), tukey_res.df_total)
    p_val_col = pd.DataFrame(p_value)
        
    table_pd = pd.DataFrame(data = tukey_res._results_table.data[1:],
                               columns = tukey_res._results_table.data[0])
        
    table_pd["P-value"] = p_val_col
    
    return table_pd


# In[37]:


tukey_whole_gc_total_one_seq, tukey_whole_gc_12_one_seq, tukey_whole_gc_3_one_seq = post_hoc(whole_gc_alfa_beta_gamma_one_seq)
tukey_table_whole_gc_total_one_seq = post_hoc_results_table(tukey_whole_gc_total_one_seq)
tukey_table_whole_gc_12_one_seq = post_hoc_results_table(tukey_whole_gc_12_one_seq)
tukey_table_whole_gc_3_one_seq = post_hoc_results_table(tukey_whole_gc_3_one_seq)
tukey_cat_gc_total_one_seq, tukey_cat_gc_12_one_seq, tukey_cat_gc_3_one_seq = post_hoc(cat_gc_alfa_beta_gamma_one_seq)
tukey_table_cat_gc_total_one_seq = post_hoc_results_table(tukey_cat_gc_total_one_seq)
tukey_table_cat_gc_12_one_seq = post_hoc_results_table(tukey_cat_gc_12_one_seq)
tukey_table_cat_gc_3_one_seq = post_hoc_results_table(tukey_cat_gc_3_one_seq)


# In[44]:


with open("./results/anova_tukey_gc_results_one_seq", "w") as file:
    
    file.write("Genoma completo. 1 secuencia \n")
    file.write("\nG+C total\n")
    file.write((anova_whole_gc_total_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_whole_gc_total_one_seq).to_string())
    
    file.write("\n\nGC 1-2\n")
    file.write((anova_whole_gc_12_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_whole_gc_12_one_seq).to_string())
    
    file.write("\n\nGC3\n")
    file.write((anova_whole_gc_3_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_whole_gc_3_one_seq).to_string())
    
    file.write("\n\n\nGenes concatenados. 1 secuencia\n")
    file.write("\nG+C total\n")
    file.write((anova_cat_gc_total_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_cat_gc_total_one_seq).to_string())
    
    file.write("\n\nGC 1-2\n")
    file.write((anova_cat_gc_12_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_cat_gc_12_one_seq).to_string())

    file.write("\n\nGC3\n")
    file.write((anova_cat_gc_3_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_cat_gc_3_one_seq).to_string())
    file.write("\n")


# In[42]:


tukey_table_cat_gc_3_one_seq


# ### Obtener la longitud total de las secuencias

# In[14]:


def count_len(genomes_dict):
    
    total_beta = 0
    total_gamma = 0
    
    for virus_genus in genomes_dict.keys():
        
        if virus_genus == "Beta":
            
            for virus_type in genomes_dict[virus_genus].keys():
                
                total_beta += len(genomes_dict[virus_genus][virus_type])
        
        if virus_genus == "Gamma":
            
            for virus_type in genomes_dict[virus_genus].keys():
                
                total_gamma += len(genomes_dict[virus_genus][virus_type])
    
    return(total_beta, total_gamma)


# In[15]:


total_beta_whole, total_gamma_whole = count_len(whole_genome_one_seq)
print("Total longitud genoma completo Beta: ", total_beta_whole)
print("Total longitud genoma completo Gamma: ", total_gamma_whole)


# In[17]:


total_beta_cat, total_gamma_cat = count_len(genome_genes_cat_one_seq)
print("Total longitud genoma concatenado Beta: ", total_beta_cat)
print("Total longitud genoma concatenado Gamma: ", total_gamma_cat)


# In[23]:


len(whole_genome_one_seq["Gamma"])


# In[8]:


total_e4_beta = 0
total_e5_beta = 0
for virus_type in genome_genes["Beta"].keys():
    
    for seqrecord in genome_genes["Beta"][virus_type]:
        
        if seqrecord.name == "E4":
            total_e4_beta += len(seqrecord.seq)
        if seqrecord.name == "E5":
            total_e5_beta += len(seqrecord.seq)
            

total_e4_gamma = 0
total_e5_gamma = 0

for virus_type in genome_genes["Gamma"].keys():
    
    for seqrecord in genome_genes["Gamma"][virus_type]:
        
        if seqrecord.name == "E4":
            total_e4_gamma += len(seqrecord.seq)
        if seqrecord.name == "E5":
            total_e5_gamma += len(seqrecord.seq)


# In[33]:


for seqrecord in genome_genes["Beta"][virus_type]:
    
    print(seqrecord.name)
    print(len(seqrecord.seq))


# In[10]:


print("Total longitud E4 Beta: ", total_e4_beta)
print("Total longitud E4 Gamma: ", total_e4_gamma)


# In[11]:


print("Total longitud E5 Beta: ", total_e5_beta)
print("Total longitud E5 Gamma: ", total_e5_gamma)


# In[45]:


for virus_type in whole_genome_one_seq["Beta"].keys():
    
    print("Tipo: " + virus_type)
    print(len(whole_genome_one_seq["Beta"][virus_type].seq))


# ## Dinucleótidos

# In[6]:


nucl = ["A", "C", "T", "G"]
dinucl = [p for p in itertools.product(nucl, repeat = 2)]
dinucleotidos = []

for pair in dinucl:
    dinucleotidos.append(pair[0]+pair[1])
print(dinucleotidos)


# In[7]:


def freq_overlapping(string, substring):
    
    total = 0
    start = 0
    frequency = 0
    
    while start < len(string):
        
        position = string.find(substring, start)
    
        if position != -1:
            
            start = position + 1
            
            total += 1
        
        else:
            
            frequency = total / (len(string)-1)
            
            return frequency
            


# In[8]:


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
        
            sequence = Seq(str((genome[virus_genus][virus_type]).seq))
            counter = Counter(sequence)
    
    
            #C_total_list.append((counter["C"])/len(sequence))
            #G_total_list.append((counter["G"])/len(sequence))
            #T_total_list.append((counter["T"])/len(sequence))
            #A_total_list.append((counter["A"])/len(sequence))
    
        
    
            for pair in dinucleotidos:
            
                frequency = freq_overlapping(sequence, pair)

                din_dict[pair].append(frequency)
        
            C_total[virus_genus][virus_type] = (counter["C"])/len(sequence)
            G_total[virus_genus][virus_type] = (counter["G"])/len(sequence)
            T_total[virus_genus][virus_type] = (counter["T"])/len(sequence)
            A_total[virus_genus][virus_type] = (counter["A"])/len(sequence)
        
        
            dinucleotidos_por_tipo[virus_genus][virus_type] = din_dict
            
    return dinucleotidos_por_tipo, C_total, G_total, T_total, A_total


# In[25]:


dinucleotidos_por_tipo_whole_one_seq, C_total_whole_one_seq, G_total_whole_one_seq, T_total_whole_one_seq, A_total_whole_one_seq = freq_dinucleotidos(whole_genome_one_seq)


# In[26]:


dinucleotidos_por_tipo_cat_one_seq, C_total_cat_one_seq, G_total_cat_one_seq, T_total_cat_one_seq, A_total_cat_one_seq = freq_dinucleotidos(genome_genes_cat_one_seq)


# In[30]:


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
        
                record = genome_by_genes[virus_genus][virus_type][gene]
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


# In[31]:


dinucleotidos_por_tipo_gen_one_seq, C_total_gen_one_seq, G_total_gen_one_seq, T_total_gen_one_seq,     A_total_gen_one_seq = freq_dinucleotidos_genes(genomes_dict_by_type_and_gen_gc3s_one_seq)


# In[57]:


def freq_dinucleotidos_123_genome(genome):
    
    dinucleotidos_12 = defaultdict(dict)
    dinucleotidos_23 = defaultdict(dict)
    
    for virus_genus in genome.keys():
                
        for virus_type in genome[virus_genus].keys():
            
            din_dict_12 = defaultdict(list)
            din_dict_23 = defaultdict(list)

            sequ = Seq(str((genome[virus_genus][virus_type]).seq))
                
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


# In[58]:


dinucleotidos_12_whole_one_seq, dinucleotidos_23_whole_one_seq = freq_dinucleotidos_123_genome(whole_genome_one_seq)
dinucleotidos_12_cat_one_seq, dinucleotidos_23_cat_one_seq = freq_dinucleotidos_123_genome(genome_genes_cat_one_seq)


# In[60]:


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
            
                sequence = genome[virus_genus][virus_type][gen]
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


# In[61]:


dinucleotidos_12_gen_one_seq, dinucleotidos_23_gen_one_seq = freq_dinucleotidos_123_genes(genomes_dict_by_type_and_gen_gc3s_one_seq)


# In[12]:


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


# In[46]:


def obs_esp_genoma(frequency_dinucleotidos, C_total, G_total, T_total, A_total):

    
    obs_esp_dinucl_por_tipo = defaultdict(dict)

    for virus_genus in frequency_dinucleotidos.keys():
    
        for virus_type in frequency_dinucleotidos[virus_genus].keys():
        
            din_dict = defaultdict(list)
        
            for pair in dinucleotidos:
            
                x, y = value_dinucleotidos(pair, C_total, G_total, T_total, A_total, din_dict)
              
                obs_esp = frequency_dinucleotidos[virus_genus][virus_type][pair][0] /                                             (x[virus_genus][virus_type] * y[virus_genus][virus_type])
                din_dict[pair] = obs_esp
            
            obs_esp_dinucl_por_tipo[virus_genus][virus_type] = din_dict
    return obs_esp_dinucl_por_tipo


# In[47]:


obs_esp_dinucleotidos_whole_one_seq = obs_esp_genoma(dinucleotidos_por_tipo_whole_one_seq, C_total_whole_one_seq, 
                                             G_total_whole_one_seq, T_total_whole_one_seq, A_total_whole_one_seq)


# In[48]:


obs_esp_dinucleotidos_cat_one_seq = obs_esp_genoma(dinucleotidos_por_tipo_cat_one_seq, C_total_cat_one_seq, 
                                             G_total_cat_one_seq, T_total_cat_one_seq, A_total_cat_one_seq)


# In[62]:


obs_esp_dinucleotidos_12_whole_one_seq = obs_esp_genoma(dinucleotidos_12_whole_one_seq, C_total_whole_one_seq, 
                                             G_total_whole_one_seq, T_total_whole_one_seq, A_total_whole_one_seq)
obs_esp_dinucleotidos_23_whole_one_seq = obs_esp_genoma(dinucleotidos_23_whole_one_seq, C_total_whole_one_seq, 
                                             G_total_whole_one_seq, T_total_whole_one_seq, A_total_whole_one_seq)


# In[63]:


obs_esp_dinucleotidos_12_cat_one_seq = obs_esp_genoma(dinucleotidos_12_cat_one_seq, C_total_cat_one_seq, 
                                             G_total_cat_one_seq, T_total_cat_one_seq, A_total_cat_one_seq)
obs_esp_dinucleotidos_23_cat_one_seq = obs_esp_genoma(dinucleotidos_23_cat_one_seq, C_total_cat_one_seq, 
                                             G_total_cat_one_seq, T_total_cat_one_seq, A_total_cat_one_seq)


# In[64]:


def obs_esp_genes(frequency_dinucleotidos, C_total, G_total, T_total, A_total):

    
    obs_esp_dinucl_por_tipo = defaultdict(dict)

    for virus_genus in frequency_dinucleotidos.keys():
    
        for virus_type in frequency_dinucleotidos[virus_genus].keys():
            
            din_dict_genes = defaultdict(dict)
        
            for gene in frequency_dinucleotidos[virus_genus][virus_type].keys():
                
                din_dict = defaultdict(list)
        
                for pair in dinucleotidos:
            
                    x, y = value_dinucleotidos(pair, C_total, G_total, T_total, A_total, din_dict)
                    
                    obs_esp = frequency_dinucleotidos[virus_genus][virus_type][gene][pair][0] /                                                 (x[virus_genus][virus_type][gene][0] * 
                                                 y[virus_genus][virus_type][gene][0])
                    din_dict[pair].append(obs_esp)
                
                din_dict_genes[gene] = din_dict
            
            obs_esp_dinucl_por_tipo[virus_genus][virus_type] = din_dict_genes
                
    return obs_esp_dinucl_por_tipo


# In[50]:


obs_esp_dinucleotidos_gen_one_seq = obs_esp_genes(dinucleotidos_por_tipo_gen_one_seq, C_total_gen_one_seq, 
                                             G_total_gen_one_seq, T_total_gen_one_seq, A_total_gen_one_seq)


# In[65]:


obs_esp_dinucleotidos_12_gen_one_seq = obs_esp_genes(dinucleotidos_12_gen_one_seq, C_total_gen_one_seq, 
                                             G_total_gen_one_seq, T_total_gen_one_seq, A_total_gen_one_seq)
obs_esp_dinucleotidos_23_gen_one_seq = obs_esp_genes(dinucleotidos_23_gen_one_seq, C_total_gen_one_seq, 
                                             G_total_gen_one_seq, T_total_gen_one_seq, A_total_gen_one_seq)


# In[22]:


def dataframe_dinucleotidos(valores):
    
    df_dinucleotidos = pd.DataFrame()
    
    for virus_genus in valores.keys():
        
        df_din = pd.DataFrame.from_dict(valores[virus_genus])
        df_din = df_din.transpose().reset_index()
        
        col_genus = [virus_genus] * len(valores[virus_genus].keys())
    
        df_din["Género"] = col_genus
        
        df_din.rename(columns = {"index": "Tipo"}, inplace = True)
        
        df_din = df_din[["Género", "Tipo"] + [col for col in df_din if col not in ["Género", "Tipo"]]]
                
        df_dinucleotidos = pd.concat([df_dinucleotidos, df_din], axis = 0)
        
    return(df_dinucleotidos)


# In[70]:


df_din_whole_one_seq = dataframe_dinucleotidos(obs_esp_dinucleotidos_whole_one_seq)
df_din_12_whole_one_seq = dataframe_dinucleotidos(obs_esp_dinucleotidos_12_whole_one_seq)
df_din_23_whole_one_seq = dataframe_dinucleotidos(obs_esp_dinucleotidos_23_whole_one_seq)

df_din_cat_one_seq = dataframe_dinucleotidos(obs_esp_dinucleotidos_cat_one_seq)
df_din_12_cat_one_seq = dataframe_dinucleotidos(obs_esp_dinucleotidos_12_cat_one_seq)
df_din_23_cat_one_seq = dataframe_dinucleotidos(obs_esp_dinucleotidos_23_cat_one_seq)


# In[71]:


df_din_23_whole_one_seq


# In[74]:


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


# In[78]:


df_din_12_gen_one_seq = df_genomes_genes_gc_original(obs_esp_dinucleotidos_12_gen_one_seq)


# In[79]:


df_din_12_gen_one_seq


# In[80]:


df_din_23_gen_one_seq = df_genomes_genes_gc_original(obs_esp_dinucleotidos_23_gen_one_seq)
df_din_gen_one_seq = df_genomes_genes_gc_original(obs_esp_dinucleotidos_gen_one_seq)


# In[81]:


df_din_12_whole_one_seq.to_csv("./results/whole_genome_din_12_one_seq.csv", sep  = "\t", index = False, header = True)
df_din_12_cat_one_seq.to_csv("./results/genome_cat_din_12_one_seq.csv", sep  = "\t", index = False, header = True)
df_din_23_whole_one_seq.to_csv("./results/whole_genome_din_23_one_seq.csv", sep  = "\t", index = False, header = True)
df_din_23_cat_one_seq.to_csv("./results/genome_cat_din_23_one_seq.csv", sep  = "\t", index = False, header = True)
df_din_whole_one_seq.to_csv("./results/whole_genome_din_one_seq.csv", sep = "\t", index = False, header = True)
df_din_cat_one_seq.to_csv("./results/genome_cat_din_one_seq.csv", sep = "\t", index = False, header = True)
df_din_gen_one_seq.to_csv("./results/genes_din_one_seq.csv", sep  = "\t", index = False, header = True)
df_din_12_gen_one_seq.to_csv("./results/genes_din_12_one_seq.csv", sep  = "\t", index = False, header = True)
df_din_23_gen_one_seq.to_csv("./results/genes_din_23_one_seq.csv", sep  = "\t", index = False, header = True)


# In[87]:


df_din_whole_one_seq = pd.read_csv("./results/whole_genome_din_one_seq.csv", sep = "\t")
alfa_beta_gamma_din_whole_one_seq = df_din_whole_one_seq[df_din_whole_one_seq["Género"].isin(["Alfa", "Beta", "Gamma"])]
df_din_12_whole_one_seq = pd.read_csv("./results/whole_genome_din_12_one_seq.csv", sep = "\t")
alfa_beta_gamma_din_12_whole_one_seq = df_din_12_whole_one_seq[df_din_12_whole_one_seq["Género"].isin(["Alfa", "Beta", "Gamma"])]
df_din_23_whole_one_seq = pd.read_csv("./results/whole_genome_din_23_one_seq.csv", sep = "\t")
alfa_beta_gamma_din_23_whole_one_seq = df_din_23_whole_one_seq[df_din_23_whole_one_seq["Género"].isin(["Alfa", "Beta", "Gamma"])]
df_din_cat_one_seq = pd.read_csv("./results/genome_cat_din_one_seq.csv", sep = "\t")
alfa_beta_gamma_din_cat_one_seq = df_din_cat_one_seq[df_din_cat_one_seq["Género"].isin(["Alfa", "Beta", "Gamma"])]
df_din_12_cat_one_seq = pd.read_csv("./results/genome_cat_din_12_one_seq.csv", sep = "\t")
alfa_beta_gamma_din_12_cat_one_seq = df_din_12_cat_one_seq[df_din_12_cat_one_seq["Género"].isin(["Alfa", "Beta", "Gamma"])]
df_din_23_cat_one_seq = pd.read_csv("./results/genome_cat_din_23_one_seq.csv", sep = "\t")
alfa_beta_gamma_din_23_cat_one_seq = df_din_23_cat_one_seq[df_din_23_cat_one_seq["Género"].isin(["Alfa", "Beta", "Gamma"])]


# In[83]:


df_din_whole_one_seq


# In[14]:


### Pequeño inciso para obtener los datos de todas las secuencias

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

dinucleotidos_por_tipo_whole, C_total_whole, G_total_whole, T_total_whole, A_total_whole = freq_dinucleotidos(whole_genome)

dinucleotidos_por_tipo_cat, C_total_cat, G_total_cat, T_total_cat, A_total_cat = freq_dinucleotidos(genome_genes_cat)

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

dinucleotidos_12_whole, dinucleotidos_23_whole = freq_dinucleotidos_123_genome(whole_genome)
dinucleotidos_12_cat, dinucleotidos_23_cat = freq_dinucleotidos_123_genome(genome_genes_cat)

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

obs_esp_dinucleotidos_whole = obs_esp_genoma(dinucleotidos_por_tipo_whole, C_total_whole, 
                                             G_total_whole, T_total_whole, A_total_whole)

obs_esp_dinucleotidos_cat = obs_esp_genoma(dinucleotidos_por_tipo_cat, C_total_cat, 
                                             G_total_cat, T_total_cat, A_total_cat)


# In[29]:


obs_esp_dinucleotidos_12_whole = obs_esp_genoma(dinucleotidos_12_whole, C_total_whole, 
                                             G_total_whole, T_total_whole, A_total_whole)
obs_esp_dinucleotidos_23_whole = obs_esp_genoma(dinucleotidos_23_whole, C_total_whole, 
                                             G_total_whole, T_total_whole, A_total_whole)
obs_esp_dinucleotidos_12_cat = obs_esp_genoma(dinucleotidos_12_cat, C_total_cat, 
                                             G_total_cat, T_total_cat, A_total_cat)
obs_esp_dinucleotidos_23_cat = obs_esp_genoma(dinucleotidos_23_cat, C_total_cat, 
                                             G_total_cat, T_total_cat, A_total_cat)

def stats_dinucleotidos_genome(obs_esp):
    
    medias_dinucleotidos = defaultdict(dict)
    desv_dinucleotidos = defaultdict(dict)

    for virus_genus in obs_esp.keys():
    
        for virus_type in obs_esp[virus_genus].keys():
        
            din_dict_med = {}
            din_dict_desv = {}
        
            for pair in dinucleotidos:
            
                media_pair = np.mean(obs_esp[virus_genus][virus_type][pair])
                desv_pair = np.std(obs_esp[virus_genus][virus_type][pair])
            
                din_dict_med[pair] = media_pair
                din_dict_desv[pair] = desv_pair
            
            medias_dinucleotidos[virus_genus][virus_type] = din_dict_med
            desv_dinucleotidos[virus_genus][virus_type] = din_dict_desv
    
    return medias_dinucleotidos, desv_dinucleotidos

def dataframe_dinucleotidos(valores):
    
    df_dinucleotidos = pd.DataFrame()
    
    for virus_genus in valores.keys():
        
        df_din = pd.DataFrame.from_dict(valores[virus_genus])
        df_din = df_din.transpose().reset_index()
        
        col_genus = [virus_genus] * len(valores[virus_genus].keys())
    
        df_din["Género"] = col_genus
        
        df_din.rename(columns = {"index": "Tipo"}, inplace = True)
        
        df_din = df_din[["Género", "Tipo"] + [col for col in df_din if col not in ["Género", "Tipo"]]]
                
        df_dinucleotidos = pd.concat([df_dinucleotidos, df_din], axis = 0)
        
    return(df_dinucleotidos)


# In[34]:


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

df_din_whole = df_genomes_dinucleotidos_original(obs_esp_dinucleotidos_whole)
df_din_12_whole = df_genomes_dinucleotidos_original(obs_esp_dinucleotidos_12_whole)
df_din_23_whole = df_genomes_dinucleotidos_original(obs_esp_dinucleotidos_23_whole)

df_din_cat = df_genomes_dinucleotidos_original(obs_esp_dinucleotidos_cat)
df_din_12_cat = df_genomes_dinucleotidos_original(obs_esp_dinucleotidos_12_cat)
df_din_23_cat = df_genomes_dinucleotidos_original(obs_esp_dinucleotidos_23_cat)

df_din_whole.to_csv("./results/whole_genome_din_original_cambios.csv", sep  = "\t", 
                                     index = False, header = True)
df_din_12_whole.to_csv("./results/whole_genome_din_12_original_cambios.csv", sep  = "\t", 
                                     index = False, header = True)
df_din_23_whole.to_csv("./results/whole_genome_din_23_original_cambios.csv", sep  = "\t", 
                                     index = False, header = True)

df_din_cat.to_csv("./results/genome_cat_din_original_cambios.csv", sep  = "\t", 
                                     index = False, header = True)
df_din_12_cat.to_csv("./results/genome_cat_din_12_original_cambios.csv", sep  = "\t", 
                                     index = False, header = True)
df_din_23_cat.to_csv("./results/genome_cat_din_23_original_cambios.csv", sep  = "\t", 
                                     index = False, header = True)


# In[21]:


medias_dinucleotidos_whole, desv_dinucleotidos_whole = stats_dinucleotidos_genome(obs_esp_dinucleotidos_whole)
medias_dinucleotidos_cat, desv_dinucleotidos_cat = stats_dinucleotidos_genome(obs_esp_dinucleotidos_cat)
medias_dinucleotidos_12_whole, desv_dinucleotidos_12_whole = stats_dinucleotidos_genome(obs_esp_dinucleotidos_12_whole)
medias_dinucleotidos_12_cat, desv_dinucleotidos_12_cat = stats_dinucleotidos_genome(obs_esp_dinucleotidos_12_cat)
medias_dinucleotidos_23_whole, desv_dinucleotidos_23_whole = stats_dinucleotidos_genome(obs_esp_dinucleotidos_23_whole)
medias_dinucleotidos_23_cat, desv_dinucleotidos_23_cat = stats_dinucleotidos_genome(obs_esp_dinucleotidos_23_cat)

df_medias_dinucleotidos_whole = dataframe_dinucleotidos(medias_dinucleotidos_whole)
df_desviaciones_dinucleotidos_whole = dataframe_dinucleotidos(desv_dinucleotidos_whole)
df_medias_dinucleotidos_12_whole = dataframe_dinucleotidos(medias_dinucleotidos_12_whole)
df_desviaciones_dinucleotidos_12_whole = dataframe_dinucleotidos(desv_dinucleotidos_12_whole)
df_medias_dinucleotidos_23_whole = dataframe_dinucleotidos(medias_dinucleotidos_23_whole)
df_desviaciones_dinucleotidos_23_whole = dataframe_dinucleotidos(desv_dinucleotidos_23_whole)

df_medias_dinucleotidos_cat = dataframe_dinucleotidos(medias_dinucleotidos_cat)
df_desviaciones_dinucleotidos_cat = dataframe_dinucleotidos(desv_dinucleotidos_cat)
df_medias_dinucleotidos_12_cat = dataframe_dinucleotidos(medias_dinucleotidos_12_cat)
df_desviaciones_dinucleotidos_12_cat = dataframe_dinucleotidos(desv_dinucleotidos_12_cat)
df_medias_dinucleotidos_23_cat = dataframe_dinucleotidos(medias_dinucleotidos_23_cat)
df_desviaciones_dinucleotidos_23_cat = dataframe_dinucleotidos(desv_dinucleotidos_23_cat)


# In[23]:


df_medias_dinucleotidos_whole = dataframe_dinucleotidos(medias_dinucleotidos_whole)
df_desviaciones_dinucleotidos_whole = dataframe_dinucleotidos(desv_dinucleotidos_whole)
df_medias_dinucleotidos_12_whole = dataframe_dinucleotidos(medias_dinucleotidos_12_whole)
df_desviaciones_dinucleotidos_12_whole = dataframe_dinucleotidos(desv_dinucleotidos_12_whole)
df_medias_dinucleotidos_23_whole = dataframe_dinucleotidos(medias_dinucleotidos_23_whole)
df_desviaciones_dinucleotidos_23_whole = dataframe_dinucleotidos(desv_dinucleotidos_23_whole)

df_medias_dinucleotidos_cat = dataframe_dinucleotidos(medias_dinucleotidos_cat)
df_desviaciones_dinucleotidos_cat = dataframe_dinucleotidos(desv_dinucleotidos_cat)
df_medias_dinucleotidos_12_cat = dataframe_dinucleotidos(medias_dinucleotidos_12_cat)
df_desviaciones_dinucleotidos_12_cat = dataframe_dinucleotidos(desv_dinucleotidos_12_cat)
df_medias_dinucleotidos_23_cat = dataframe_dinucleotidos(medias_dinucleotidos_23_cat)
df_desviaciones_dinucleotidos_23_cat = dataframe_dinucleotidos(desv_dinucleotidos_23_cat)


# In[24]:


df_medias_dinucleotidos_cat


# In[25]:


df_medias_dinucleotidos_whole.to_csv("./results/medias_dinucleotidos_whole_medias_tipos.csv", sep  = "\t", 
                                     index = False, header = True)
df_desviaciones_dinucleotidos_whole.to_csv("./results/desviaciones_dinucleotidos_whole_medias_tipos.csv", sep = "\t", 
                                     index = False, header = True)
df_medias_dinucleotidos_12_whole.to_csv("./results/medias_dinucleotidos_12_whole_medias_tipos.csv", sep  = "\t", 
                                     index = False, header = True)
df_desviaciones_dinucleotidos_12_whole.to_csv("./results/desviaciones_dinucleotidos_12_whole_medias_tipos.csv", sep = "\t", 
                                     index = False, header = True)
df_medias_dinucleotidos_23_whole.to_csv("./results/medias_dinucleotidos_23_whole_medias_tipos.csv", sep  = "\t", 
                                     index = False, header = True)
df_desviaciones_dinucleotidos_23_whole.to_csv("./results/desviaciones_dinucleotidos_23_whole_medias_tipos.csv", sep = "\t", 
                                     index = False, header = True)

df_medias_dinucleotidos_cat.to_csv("./results/medias_dinucleotidos_cat_medias_tipos.csv", sep  = "\t", 
                                     index = False, header = True)
df_desviaciones_dinucleotidos_cat.to_csv("./results/desviaciones_dinucleotidos_cat_medias_tipos.csv", sep = "\t", 
                                     index = False, header = True)
df_medias_dinucleotidos_12_cat.to_csv("./results/medias_dinucleotidos_12_cat_medias_tipos.csv", sep  = "\t", 
                                     index = False, header = True)
df_desviaciones_dinucleotidos_12_cat.to_csv("./results/desviaciones_dinucleotidos_12_cat_medias_tipos.csv", sep = "\t", 
                                     index = False, header = True)
df_medias_dinucleotidos_23_cat.to_csv("./results/medias_dinucleotidos_23_cat_medias_tipos.csv", sep  = "\t", 
                                     index = False, header = True)
df_desviaciones_dinucleotidos_23_cat.to_csv("./results/desviaciones_dinucleotidos_23_cat_medias_tipos.csv", sep = "\t", 
                                     index = False, header = True)


# In[89]:


def anovas_din(data_df, pair):
    
    comparison = pair + " ~ Género"
    
    mod_anova_pair = ols(comparison, data = data_df).fit()
 
    anova_table_pair = anova.anova_lm(mod_anova_pair, typ = 2)
 
    return anova_table_pair


# In[90]:


def post_hoc_din(data_df, pair):
    
    result_tukey_pair = (pairwise_tukeyhsd(data_df[pair], data_df["Género"], alpha = 0.05))

    return result_tukey_pair


# In[92]:


def post_hoc_results_table(tukey_res): 
    p_value = psturng(np.abs(tukey_res.meandiffs / tukey_res.std_pairs),
           len(tukey_res.groupsunique), tukey_res.df_total)
    p_val_col = pd.DataFrame(p_value)
        
    table_pd = pd.DataFrame(data = tukey_res._results_table.data[1:],
                               columns = tukey_res._results_table.data[0])
        
    table_pd["P-value"] = p_val_col
    
    return table_pd


# In[161]:


with open("./results/anova_tukey_din_whole_results_one_seq", "w") as file:
    
    file.write("Genoma entero 1 secuencia. Dinucleotidos sin importar la posicion\n")
    
    for pair in dinucleotidos:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_whole_one_seq, pair)

        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_whole_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")

                
    file.write("\n\n******************** \nGenoma entero 1 secuencia. Dinucleotidos en posicion 1-2\n")
    
    for pair in dinucleotidos:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_12_whole_one_seq, pair)
        
        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_12_whole_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")

    file.write("\n\n******************** \nGenoma entero 1 secuencia. Dinucleotidos en posicion 2-3\n")
    
    for pair in dinucleotidos:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_23_whole_one_seq, pair)
        
        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_23_whole_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)

            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")


# In[162]:


with open("./results/anova_tukey_din_cat_results_one_seq", "w") as file:
    
    file.write("Genes concatenados 1 secuencia. Dinucleotidos sin importar la posicion\n")
    
    for pair in dinucleotidos:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_cat_one_seq, pair)

        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_cat_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")

                
    file.write("\n\n******************** \nGenes concatenados 1 secuencia. Dinucleotidos en posicion 1-2\n")
    
    for pair in dinucleotidos:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_12_cat_one_seq, pair)
        
        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_12_cat_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")

    file.write("\n\n******************** \nGenes concatenados 1 secuencia. Dinucleotidos en posicion 2-3\n")
    
    for pair in dinucleotidos:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_23_cat_one_seq, pair)
        
        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_23_cat_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)

            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")


# In[116]:


def medias_dataframe_din(dataframe, dinucleotidos):
    
    means = pd.DataFrame()
    standard_error_din = pd.DataFrame()
    
    for genero in dataframe["Género"].unique():
        
        mean_results = pd.DataFrame.mean(dataframe[dataframe["Género"].isin([genero])], 0)
        
        mean_results = mean_results.to_frame()
        mean_results = mean_results.rename(columns = {0: genero})
        mean_results = mean_results.drop(mean_results.index[0])
        means = pd.concat([means, mean_results], axis = 1)
        
        std_err_din = pd.DataFrame()
        
        for pair in dinucleotidos:
            
            standard_error = stats.sem(dataframe[pair][dataframe["Género"].isin([genero])])
        
            std_err = pd.DataFrame({pair : [standard_error]}, index = [genero])
            
            std_err = std_err.transpose()

            std_err_din = pd.concat([std_err_din, std_err], axis = 0)

        standard_error_din = pd.concat([standard_error_din, std_err_din], axis = 1)
        
    return means, standard_error_din


# In[128]:


medias_din_whole_one_seq, error_din_whole_one_seq = medias_dataframe_din(alfa_beta_gamma_din_whole_one_seq, dinucleotidos)
medias_din_12_whole_one_seq, error_din_12_whole_one_seq = medias_dataframe_din(alfa_beta_gamma_din_12_whole_one_seq, dinucleotidos)
medias_din_23_whole_one_seq, error_din_23_whole_one_seq = medias_dataframe_din(alfa_beta_gamma_din_23_whole_one_seq, dinucleotidos)
medias_din_cat_one_seq, error_din_cat_one_seq = medias_dataframe_din(alfa_beta_gamma_din_cat_one_seq, dinucleotidos)
medias_din_12_cat_one_seq, error_din_12_cat_one_seq = medias_dataframe_din(alfa_beta_gamma_din_12_cat_one_seq, dinucleotidos)
medias_din_23_cat_one_seq, error_din_23_cat_one_seq = medias_dataframe_din(alfa_beta_gamma_din_23_cat_one_seq, dinucleotidos)


# In[129]:


error_din_whole_one_seq


# In[119]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 30
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.25, hspace = 0.25)

fig = plt.figure()
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

x = np.linspace(16, -1, 16)
y = [1] * 16

medias_din_whole_one_seq.plot.bar(y = ["Alfa", "Beta", "Gamma"], ax = ax1, yerr = error_din_whole_one_seq)
ax1.set_title("Genoma completo 1 secuencia. Dinucleótidos sin importar la posición", fontsize = 30)
ax1.tick_params(axis = "x", labelsize = 20, rotation = 0)
ax1.tick_params(axis = "y", labelsize = 20, rotation = 0)
ax1.legend(fontsize = 20)
ax1.plot(x, y, color = "black", linestyle = "dashed")

medias_din_12_whole_one_seq.plot.bar(y = ["Alfa", "Beta", "Gamma"], ax = ax2, yerr = error_din_12_whole_one_seq)
ax2.set_title("Genoma completo 1 secuencia. Dinucleótidos en posición 1-2", fontsize = 30)
ax2.tick_params(axis = "x", labelsize = 20, rotation = 0)
ax2.tick_params(axis = "y", labelsize = 20, rotation = 0)
ax2.legend(fontsize = 20)

ax2.plot(x, y, color = "black", linestyle = "dashed")

medias_din_23_whole_one_seq.plot.bar(y = ["Alfa", "Beta", "Gamma"], ax = ax3, yerr = error_din_23_whole_one_seq)
ax3.set_title("Genoma completo 1 secuencia. Dinucleótidos en posición 2-3", fontsize = 30)
ax3.tick_params(axis = "x", labelsize = 20, rotation = 0)
ax3.tick_params(axis = "y", labelsize = 20, rotation = 0)
ax3.legend(fontsize = 20)

ax3.plot(x, y, color = "black", linestyle = "dashed")

plt.savefig("./results/din_whole_one_seq.png")

plt.show()


# In[120]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 30
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.25, hspace = 0.25)

fig = plt.figure()
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

x = np.linspace(16, -1, 16)
y = [1] * 16

medias_din_cat_one_seq.plot.bar(y = ["Alfa", "Beta", "Gamma"], ax = ax1, yerr = error_din_cat_one_seq)
ax1.set_title("Genes concatenados 1 secuencia. Dinucleótidos sin importar la posición", fontsize = 30)
ax1.tick_params(axis = "x", labelsize = 20, rotation = 0)
ax1.tick_params(axis = "y", labelsize = 20, rotation = 0)
ax1.legend(fontsize = 20)

ax1.plot(x, y, color = "black", linestyle = "dashed")

medias_din_12_cat_one_seq.plot.bar(y = ["Alfa", "Beta", "Gamma"], ax = ax2, yerr = error_din_12_cat_one_seq)
ax2.set_title("Genes concatenados 1 secuencia. Dinucleótidos en posición 1-2", fontsize = 30)
ax2.tick_params(axis = "x", labelsize = 20, rotation = 0)
ax2.tick_params(axis = "y", labelsize = 20, rotation = 0)
ax2.legend(fontsize = 20)

ax2.plot(x, y, color = "black", linestyle = "dashed")

medias_din_23_cat_one_seq.plot.bar(y = ["Alfa", "Beta", "Gamma"], ax = ax3, yerr = error_din_23_cat_one_seq)
ax3.set_title("Genes concatenados 1 secuencia. Dinucleótidos en posición 2-3", fontsize = 30)
ax3.tick_params(axis = "x", labelsize = 20, rotation = 0)
ax3.tick_params(axis = "y", labelsize = 20, rotation = 0)
ax3.legend(fontsize = 20)

ax3.plot(x, y, color = "black", linestyle = "dashed")

plt.savefig("./results/din_cat_one_seq.png")

plt.show()


# ### Dinucleótidos seleccionados

# In[127]:


dinucleotidos_selected = ["CG", "TC", "GA"]


# In[163]:


with open("./results/anova_tukey_selected_din_whole_results_one_seq", "w") as file:
    
    file.write("Genoma entero 1 secuencia. Dinucleotidos sin importar la posicion\n")
    
    for pair in dinucleotidos_selected:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_whole_one_seq, pair)

        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_whole_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")

                
    file.write("\n\n******************** \nGenoma entero 1 secuencia. Dinucleotidos en posicion 1-2\n")
    
    for pair in dinucleotidos_selected:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_12_whole_one_seq, pair)
        
        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_12_whole_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")

    file.write("\n\n******************** \nGenoma entero 1 secuencia. Dinucleotidos en posicion 2-3\n")
    
    for pair in dinucleotidos_selected:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_23_whole_one_seq, pair)
        
        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_23_whole_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)

            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")


# In[164]:


with open("./results/anova_tukey_selected_din_cat_results_one_seq", "w") as file:
    
    file.write("Genes concatenados 1 secuencia. Dinucleotidos sin importar la posicion\n")
    
    for pair in dinucleotidos_selected:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_cat_one_seq, pair)

        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_cat_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")

                
    file.write("\n\n******************** \nGenes concatenados 1 secuencia. Dinucleotidos en posicion 1-2\n")
    
    for pair in dinucleotidos_selected:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_12_cat_one_seq, pair)
        
        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_12_cat_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")

    file.write("\n\n******************** \nGenes concatenados 1 secuencia. Dinucleotidos en posicion 2-3\n")
    
    for pair in dinucleotidos_selected:
        
        file.write("\n" + pair +"\n")
        anova_table_pair = anovas_din(alfa_beta_gamma_din_23_cat_one_seq, pair)
        
        file.write("\nAnova\n")
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns) + "\n")
        
        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din(alfa_beta_gamma_din_23_cat_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)

            file.write("\nTukey\n")
            file.write((res_p_val).to_string() + "\n")


# In[130]:


medias_din_whole_selected_one_seq = medias_din_whole_one_seq.loc[dinucleotidos_selected]
medias_din_12_whole_selected_one_seq = medias_din_12_whole_one_seq.loc[dinucleotidos_selected]
medias_din_23_whole_selected_one_seq = medias_din_23_whole_one_seq.loc[dinucleotidos_selected]

error_din_whole_selected_one_seq = error_din_whole_one_seq.loc[dinucleotidos_selected]
error_din_12_whole_selected_one_seq = error_din_12_whole_one_seq.loc[dinucleotidos_selected]
error_din_23_whole_selected_one_seq = error_din_23_whole_one_seq.loc[dinucleotidos_selected]

medias_din_cat_selected_one_seq = medias_din_cat_one_seq.loc[dinucleotidos_selected]
medias_din_12_cat_selected_one_seq = medias_din_12_cat_one_seq.loc[dinucleotidos_selected]
medias_din_23_cat_selected_one_seq = medias_din_23_cat_one_seq.loc[dinucleotidos_selected]

error_din_cat_selected_one_seq = error_din_cat_one_seq.loc[dinucleotidos_selected]
error_din_12_cat_selected_one_seq = error_din_12_cat_one_seq.loc[dinucleotidos_selected]
error_din_23_cat_selected_one_seq = error_din_23_cat_one_seq.loc[dinucleotidos_selected]


# In[131]:


medias_din_whole_selected_one_seq


# In[133]:


colores = plt.cm.get_cmap("YlGnBu")#.colors
colores


# In[147]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 60
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.35)

barwidth = 0.25
x = np.arange(3)
y = np.linspace(3, -0.25, 3)
z = [1] * 3

plt.subplot(131)
plt.title("Genoma entero 1 secuencia. Dinucleótidos sin importar la posición", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_whole_selected_one_seq["Alfa"], width = barwidth, 
       label = "Alfa", color = colores(40), 
        yerr = error_din_whole_selected_one_seq["Alfa"])
plt.bar(x + barwidth, medias_din_whole_selected_one_seq["Beta"], width = barwidth,
       label = "Beta", color = colores(70),
       yerr = error_din_whole_selected_one_seq["Beta"])
plt.bar(x + (2 * barwidth), medias_din_whole_selected_one_seq["Gamma"], width = barwidth,
       label = "Gamma", color = colores(120),
       yerr = error_din_whole_selected_one_seq["Gamma"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.2])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_whole_selected_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.subplot(132)
plt.title("Genoma entero 1 secuencia. Dinucleótidos en posición 1-2", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_12_whole_selected_one_seq["Alfa"], width = barwidth, 
       label = "Alfa", color = colores(40), 
        yerr = error_din_12_whole_selected_one_seq["Alfa"])
plt.bar(x + barwidth, medias_din_12_whole_selected_one_seq["Beta"], width = barwidth,
       label = "Beta", color = colores(70),
       yerr = error_din_12_whole_selected_one_seq["Beta"])
plt.bar(x + (2 * barwidth), medias_din_12_whole_selected_one_seq["Gamma"], width = barwidth,
       label = "Gamma", color = colores(120),
       yerr = error_din_12_whole_selected_one_seq["Gamma"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.2])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_12_whole_selected_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.subplot(133)
plt.title("Genoma entero 1 secuencia. Dinucleótidos en posición 2-3", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_23_whole_selected_one_seq["Alfa"], width = barwidth, 
       label = "Alfa", color = colores(40), 
        yerr = error_din_23_whole_selected_one_seq["Alfa"])
plt.bar(x + barwidth, medias_din_23_whole_selected_one_seq["Beta"], width = barwidth,
       label = "Beta", color = colores(70),
       yerr = error_din_23_whole_selected_one_seq["Beta"])
plt.bar(x + (2 * barwidth), medias_din_23_whole_selected_one_seq["Gamma"], width = barwidth,
       label = "Gamma", color = colores(120),
       yerr = error_din_23_whole_selected_one_seq["Gamma"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.2])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_23_whole_selected_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.savefig("./results/selected_din_whole_one_seq.png")

plt.show()


# In[148]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 60
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.35)

barwidth = 0.25
x = np.arange(3)
y = np.linspace(3, -0.25, 3)
z = [1] * 3

plt.subplot(131)
plt.title("Genes concatenados 1 secuencia. Dinucleótidos sin importar la posición", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_cat_selected_one_seq["Alfa"], width = barwidth, 
       label = "Alfa", color = colores(40), 
        yerr = error_din_cat_selected_one_seq["Alfa"])
plt.bar(x + barwidth, medias_din_cat_selected_one_seq["Beta"], width = barwidth,
       label = "Beta", color = colores(70),
       yerr = error_din_cat_selected_one_seq["Beta"])
plt.bar(x + (2 * barwidth), medias_din_cat_selected_one_seq["Gamma"], width = barwidth,
       label = "Gamma", color = colores(120),
       yerr = error_din_cat_selected_one_seq["Gamma"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 2])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_cat_selected_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.subplot(132)
plt.title("Genes concatenados 1 secuencia. Dinucleótidos en posición 1-2", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_12_cat_selected_one_seq["Alfa"], width = barwidth, 
       label = "Alfa", color = colores(40), 
        yerr = error_din_12_cat_selected_one_seq["Alfa"])
plt.bar(x + barwidth, medias_din_12_cat_selected_one_seq["Beta"], width = barwidth,
       label = "Beta", color = colores(70),
       yerr = error_din_12_cat_selected_one_seq["Beta"])
plt.bar(x + (2 * barwidth), medias_din_12_cat_selected_one_seq["Gamma"], width = barwidth,
       label = "Gamma", color = colores(120),
       yerr = error_din_12_cat_selected_one_seq["Gamma"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 2])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_12_cat_selected_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.subplot(133)
plt.title("Genes concatenados 1 secuencia. Dinucleótidos en posición 2-3", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_23_cat_selected_one_seq["Alfa"], width = barwidth, 
       label = "Alfa", color = colores(40), 
        yerr = error_din_23_cat_selected_one_seq["Alfa"])
plt.bar(x + barwidth, medias_din_23_cat_selected_one_seq["Beta"], width = barwidth,
       label = "Beta", color = colores(70),
       yerr = error_din_23_cat_selected_one_seq["Beta"])
plt.bar(x + (2 * barwidth), medias_din_23_cat_selected_one_seq["Gamma"], width = barwidth,
       label = "Gamma", color = colores(120),
       yerr = error_din_23_cat_selected_one_seq["Gamma"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 2])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_23_cat_selected_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.savefig("./results/selected_din_cat_one_seq.png")

plt.show()


# ## En función del riesgo en Alfa

# In[20]:


alto = [39, 45, 51, 52, 56, 58, 59, 68, 73, 82, 16, 18, 31, 33, 35]
bajo = [11, 40, 42, 43, 44, 53, 54, 61, 66, 6, 70, 72, 74, 81, 89]
indefinido = [102, 106, 10, 114, 117, 125, 13, 160, 177, 26, 27, 28, 29, 2,
             30, 32, 34, 3, 57, 62, 67, 69, 71, 77, 78, 7, 83, 84, 85, 86, 
             87, 90, 91, 94, 97]


# In[154]:


df_din_gen_one_seq = pd.read_csv("./results/genes_din_one_seq.csv", sep = "\t")
alfa_din_gen_one_seq = df_din_gen_one_seq[df_din_gen_one_seq["Genero"].isin(["Alfa"])]
df_din_12_gen_one_seq = pd.read_csv("./results/genes_din_12_one_seq.csv", sep = "\t")
alfa_din_12_gen_one_seq = df_din_12_gen_one_seq[df_din_12_gen_one_seq["Genero"].isin(["Alfa"])]
df_din_23_gen_one_seq = pd.read_csv("./results/genes_din_23_one_seq.csv", sep = "\t")
alfa_din_23_gen_one_seq = df_din_23_gen_one_seq[df_din_23_gen_one_seq["Genero"].isin(["Alfa"])]


# In[179]:


alfa_din_whole_one_seq = df_din_whole_one_seq[df_din_whole_one_seq["Género"].isin(["Alfa"])]
alfa_din_12_whole_one_seq = df_din_12_whole_one_seq[df_din_12_whole_one_seq["Género"].isin(["Alfa"])]
alfa_din_23_whole_one_seq = df_din_23_whole_one_seq[df_din_23_whole_one_seq["Género"].isin(["Alfa"])]

alfa_din_cat_one_seq = df_din_cat_one_seq[df_din_cat_one_seq["Género"].isin(["Alfa"])]
alfa_din_12_cat_one_seq = df_din_12_cat_one_seq[df_din_12_cat_one_seq["Género"].isin(["Alfa"])]
alfa_din_23_cat_one_seq = df_din_23_cat_one_seq[df_din_23_cat_one_seq["Género"].isin(["Alfa"])]


# In[183]:


columnas_extraer_genes = dinucleotidos_selected[:]
columnas_extraer_genes.append("Tipo")
columnas_extraer_genes.append("Gen")

columnas_extraer_genoma = dinucleotidos_selected[:]
columnas_extraer_genoma.append("Tipo")


# In[233]:


def select_add_riesgo(df, columnas_extraer):
    
    df_selected = df[columnas_extraer]
    
    list_riesgo = []

    for virus_type in df_selected["Tipo"]:
    
        if virus_type in alto:
        
            list_riesgo.append("Alto")

        elif virus_type in bajo:
        
            list_riesgo.append("Bajo")

        elif virus_type in indefinido:
        
            list_riesgo.append("Indefinido")

        else:
            print("El tipo" + virus_type + "no está clasificado")
    
    df_selected = df_selected.assign(Riesgo = list_riesgo)
    
    return df_selected


# In[185]:


alfa_selected_din_genes_one_seq = select_add_riesgo(alfa_din_gen_one_seq, columnas_extraer_genes)
alfa_selected_din_12_genes_one_seq = select_add_riesgo(alfa_din_12_gen_one_seq, columnas_extraer_genes)
alfa_selected_din_23_genes_one_seq = select_add_riesgo(alfa_din_23_gen_one_seq, columnas_extraer_genes)

alfa_selected_din_whole_one_seq = select_add_riesgo(alfa_din_whole_one_seq, columnas_extraer_genoma)
alfa_selected_din_12_whole_one_seq = select_add_riesgo(alfa_din_12_whole_one_seq, columnas_extraer_genoma)
alfa_selected_din_23_whole_one_seq = select_add_riesgo(alfa_din_23_whole_one_seq, columnas_extraer_genoma)

alfa_selected_din_cat_one_seq = select_add_riesgo(alfa_din_cat_one_seq, columnas_extraer_genoma)
alfa_selected_din_12_cat_one_seq = select_add_riesgo(alfa_din_12_cat_one_seq, columnas_extraer_genoma)
alfa_selected_din_23_cat_one_seq = select_add_riesgo(alfa_din_23_cat_one_seq, columnas_extraer_genoma)


# In[186]:


def anovas_din_alfa(data_df, pair):
    
    comparison = pair + " ~ Riesgo"
    
    mod_anova_pair = ols(comparison, data = data_df).fit()
 
    anova_table_pair = anova.anova_lm(mod_anova_pair, typ = 2)
 
    return anova_table_pair


# In[187]:


def post_hoc_din_alfa(data_df, pair):
    
    result_tukey_pair = (pairwise_tukeyhsd(data_df[pair], data_df["Riesgo"], alpha = 0.05))

    return result_tukey_pair


# In[188]:


with open("./results/anova_tukey_selected_din_riesgo_genes_results_one_seq", "w") as file:
    
    file.write("Género Alfa 1 secuencia. Dinucleótidos sin importar la posición\n")
    
    for pair in dinucleotidos_selected:
        
        anova_table_pair = anovas_din_alfa(alfa_selected_din_genes_one_seq, pair)
        
        file.write("\n" + pair +"\n")
        file.write("\n ANOVA \n")    
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns))

        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din_alfa(alfa_selected_din_genes_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\n\n TUKEY \n")
            file.write((res_p_val).to_string() + "\n")

                
    file.write("\n\n******************** \Género Alfa 1 secuencia. Dinucleótidos en posición 1-2\n")
    
    for pair in dinucleotidos_selected:
        
        anova_table_pair = anovas_din_alfa(alfa_selected_din_12_genes_one_seq, pair)
        
        file.write("\n" + pair +"\n")
        file.write("\n ANOVA \n")    
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns))

        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din_alfa(alfa_selected_din_12_genes_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\n\n TUKEY \n")
            file.write((res_p_val).to_string() + "\n")

    file.write("\n\n******************** \nGénero Alfa 1 secuencia. Dinucleótidos en posición 2-3\n")
    
    for pair in dinucleotidos_selected:
        
        anova_table_pair = anovas_din_alfa(alfa_selected_din_23_genes_one_seq, pair)
        
        file.write("\n" + pair +"\n")
        file.write("\n ANOVA \n")    
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns))

        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din_alfa(alfa_selected_din_23_genes_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\n\n TUKEY \n")
            file.write((res_p_val).to_string() + "\n")


# In[189]:


with open("./results/anova_tukey_selected_din_riesgo_whole_results_one_seq", "w") as file:
    
    file.write("Género Alfa, genoma completo 1 secuencia. Dinucleótidos sin importar la posición\n")
    
    for pair in dinucleotidos_selected:
        
        anova_table_pair = anovas_din_alfa(alfa_selected_din_whole_one_seq, pair)
        
        file.write("\n" + pair +"\n")
        file.write("\n ANOVA \n")    
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns))

        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din_alfa(alfa_selected_din_whole_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\n\n TUKEY \n")
            file.write((res_p_val).to_string() + "\n")

                
    file.write("\n\n******************** \Género Alfa, genoma completo 1 secuencia. Dinucleótidos en posición 1-2\n")
    
    for pair in dinucleotidos_selected:
        
        anova_table_pair = anovas_din_alfa(alfa_selected_din_12_whole_one_seq, pair)
        
        file.write("\n" + pair +"\n")
        file.write("\n ANOVA \n")    
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns))

        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din_alfa(alfa_selected_din_12_whole_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\n\n TUKEY \n")
            file.write((res_p_val).to_string() + "\n")

    file.write("\n\n******************** \nGénero Alfa, genoma completo 1 secuencia. Dinucleótidos en posición 2-3\n")
    
    for pair in dinucleotidos_selected:
        
        anova_table_pair = anovas_din_alfa(alfa_selected_din_23_whole_one_seq, pair)
        
        file.write("\n" + pair +"\n")
        file.write("\n ANOVA \n")    
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns))

        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din_alfa(alfa_selected_din_23_whole_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\n\n TUKEY \n")
            file.write((res_p_val).to_string() + "\n")


# In[191]:


with open("./results/anova_tukey_selected_din_riesgo_cat_results_one_seq", "w") as file:
    
    file.write("Género Alfa, genes concatenados 1 secuencia. Dinucleótidos sin importar la posición\n")
    
    for pair in dinucleotidos_selected:
        
        anova_table_pair = anovas_din_alfa(alfa_selected_din_cat_one_seq, pair)
        
        file.write("\n" + pair +"\n")
        file.write("\n ANOVA \n")    
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns))

        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din_alfa(alfa_selected_din_cat_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\n\n TUKEY \n")
            file.write((res_p_val).to_string() + "\n")

                
    file.write("\n\n******************** \Género Alfa, genes concatenados 1 secuencia. Dinucleótidos en posición 1-2\n")
    
    for pair in dinucleotidos_selected:
        
        anova_table_pair = anovas_din_alfa(alfa_selected_din_12_cat_one_seq, pair)
        
        file.write("\n" + pair +"\n")
        file.write("\n ANOVA \n")    
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns))

        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din_alfa(alfa_selected_din_12_cat_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\n\n TUKEY \n")
            file.write((res_p_val).to_string() + "\n")

    file.write("\n\n******************** \nGénero Alfa, genes concatenados 1 secuencia. Dinucleótidos en posición 2-3\n")
    
    for pair in dinucleotidos_selected:
        
        anova_table_pair = anovas_din_alfa(alfa_selected_din_23_cat_one_seq, pair)
        
        file.write("\n" + pair +"\n")
        file.write("\n ANOVA \n")    
        file.write(tabulate(anova_table_pair, headers = anova_table_pair.columns))

        if anova_table_pair["PR(>F)"][0] < 0.05:
            result_tukey_pair = post_hoc_din_alfa(alfa_selected_din_23_cat_one_seq, pair)
            
            res_p_val = post_hoc_results_table(result_tukey_pair)
            
            file.write("\n\n TUKEY \n")
            file.write((res_p_val).to_string() + "\n")


# In[192]:


def medias_dataframe_din_riesgo(dataframe, dinucleotidos):
    
    means = pd.DataFrame()
    standard_error_din = pd.DataFrame()
    
    for riesgo in dataframe["Riesgo"].unique():
        
        mean_results = pd.DataFrame.mean(dataframe[dataframe["Riesgo"].isin([riesgo])], 0)
        
        mean_results = mean_results.to_frame()
        mean_results = mean_results.rename(columns = {0: riesgo})
        mean_results = mean_results.drop(mean_results.index[-1])
                
        means = pd.concat([means, mean_results], axis = 1)
        
        std_err_din = pd.DataFrame()
        
        for pair in dinucleotidos:
            
            standard_error = stats.sem(dataframe[pair][dataframe["Riesgo"].isin([riesgo])])
        
            std_err = pd.DataFrame({pair : [standard_error]}, index = [riesgo])
            
            std_err = std_err.transpose()

            std_err_din = pd.concat([std_err_din, std_err], axis = 0)

        standard_error_din = pd.concat([standard_error_din, std_err_din], axis = 1)
        
    return means, standard_error_din


# In[171]:


medias_din_riesgo_one_seq, error_din_riesgo_one_seq = medias_dataframe_din_riesgo(alfa_selected_din_genes_one_seq, dinucleotidos_selected)
medias_din_12_riesgo_one_seq, error_din_12_riesgo_one_seq = medias_dataframe_din_riesgo(alfa_selected_din_12_genes_one_seq, dinucleotidos_selected)
medias_din_23_riesgo_one_seq, error_din_23_riesgo_one_seq = medias_dataframe_din_riesgo(alfa_selected_din_23_genes_one_seq, dinucleotidos_selected)


# In[201]:


medias_din_riesgo_whole_one_seq, error_din_riesgo_whole_one_seq = medias_dataframe_din_riesgo(alfa_selected_din_whole_one_seq, dinucleotidos_selected)
medias_din_12_riesgo_whole_one_seq, error_din_12_riesgo_whole_one_seq = medias_dataframe_din_riesgo(alfa_selected_din_12_whole_one_seq, dinucleotidos_selected)
medias_din_23_riesgo_whole_one_seq, error_din_23_riesgo_whole_one_seq = medias_dataframe_din_riesgo(alfa_selected_din_23_whole_one_seq, dinucleotidos_selected)

medias_din_riesgo_cat_one_seq, error_din_riesgo_cat_one_seq = medias_dataframe_din_riesgo(alfa_selected_din_cat_one_seq, dinucleotidos_selected)
medias_din_12_riesgo_cat_one_seq, error_din_12_riesgo_cat_one_seq = medias_dataframe_din_riesgo(alfa_selected_din_12_cat_one_seq, dinucleotidos_selected)
medias_din_23_riesgo_cat_one_seq, error_din_23_riesgo_cat_one_seq = medias_dataframe_din_riesgo(alfa_selected_din_23_cat_one_seq, dinucleotidos_selected)


# In[195]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 60
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.35)

barwidth = 0.25
x = np.arange(3)
y = np.linspace(3, -0.25, 3)
z = [1] * 3

plt.subplot(131)
plt.title("Género Alfa 1 secuencia. Dinucleótidos sin importar la posición", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_riesgo_one_seq["Alto"], width = barwidth, 
       label = "Alto", color = colores(40), 
        yerr = error_din_riesgo_one_seq["Alto"])
plt.bar(x + barwidth, medias_din_riesgo_one_seq["Bajo"], width = barwidth,
       label = "Bajo", color = colores(70),
       yerr = error_din_riesgo_one_seq["Bajo"])
plt.bar(x + (2 * barwidth), medias_din_riesgo_one_seq["Indefinido"], width = barwidth,
       label = "Indefinido", color = colores(120),
       yerr = error_din_riesgo_one_seq["Indefinido"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.8])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_riesgo_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.subplot(132)
plt.title("Género Alfa 1 secuencia. Dinucleótidos en posición 1-2", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_12_riesgo_one_seq["Alto"], width = barwidth, 
       label = "Alto", color = colores(40), 
        yerr = error_din_12_riesgo_one_seq["Alto"])
plt.bar(x + barwidth, medias_din_12_riesgo_one_seq["Bajo"], width = barwidth,
       label = "Bajo", color = colores(70),
       yerr = error_din_12_riesgo_one_seq["Bajo"])
plt.bar(x + (2 * barwidth), medias_din_12_riesgo_one_seq["Indefinido"], width = barwidth,
       label = "indefinido", color = colores(120),
       yerr = error_din_12_riesgo_one_seq["Indefinido"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.8])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_12_riesgo_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.subplot(133)
plt.title("Género Alfa 1 secuencia. Dinucleótidos en posición 2-3", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_23_riesgo_one_seq["Alto"], width = barwidth, 
       label = "Alto", color = colores(40), 
        yerr = error_din_23_riesgo_one_seq["Alto"])
plt.bar(x + barwidth, medias_din_23_riesgo_one_seq["Bajo"], width = barwidth,
       label = "Bajo", color = colores(70),
       yerr = error_din_23_riesgo_one_seq["Bajo"])
plt.bar(x + (2 * barwidth), medias_din_23_riesgo_one_seq["Indefinido"], width = barwidth,
       label = "Indefinido", color = colores(120),
       yerr = error_din_23_riesgo_one_seq["Indefinido"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.8])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_23_riesgo_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.savefig("./results/selected_din_riesgo_one_seq.png")

plt.show()


# In[198]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 50
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.30)

barwidth = 0.25
x = np.arange(3)
y = np.linspace(3, -0.25, 3)
z = [1] * 3

plt.subplot(131)
plt.title("Género Alfa, genoma completo 1 secuencia. \nDinucleótidos sin importar la posición", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_riesgo_whole_one_seq["Alto"], width = barwidth, 
       label = "Alto", color = colores(40), 
        yerr = error_din_riesgo_whole_one_seq["Alto"])
plt.bar(x + barwidth, medias_din_riesgo_whole_one_seq["Bajo"], width = barwidth,
       label = "Bajo", color = colores(70),
       yerr = error_din_riesgo_whole_one_seq["Bajo"])
plt.bar(x + (2 * barwidth), medias_din_riesgo_whole_one_seq["Indefinido"], width = barwidth,
       label = "Indefinido", color = colores(120),
       yerr = error_din_riesgo_whole_one_seq["Indefinido"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.1])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_riesgo_whole_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.subplot(132)
plt.title("Género Alfa, genoma completo 1 secuencia. \nDinucleótidos en posición 1-2", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_12_riesgo_whole_one_seq["Alto"], width = barwidth, 
       label = "Alto", color = colores(40), 
        yerr = error_din_12_riesgo_whole_one_seq["Alto"])
plt.bar(x + barwidth, medias_din_12_riesgo_whole_one_seq["Bajo"], width = barwidth,
       label = "Bajo", color = colores(70),
       yerr = error_din_12_riesgo_whole_one_seq["Bajo"])
plt.bar(x + (2 * barwidth), medias_din_12_riesgo_whole_one_seq["Indefinido"], width = barwidth,
       label = "indefinido", color = colores(120),
       yerr = error_din_12_riesgo_whole_one_seq["Indefinido"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.1])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_12_riesgo_whole_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.subplot(133)
plt.title("Género Alfa, genoma completo 1 secuencia. \nDinucleótidos en posición 2-3", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_23_riesgo_whole_one_seq["Alto"], width = barwidth, 
       label = "Alto", color = colores(40), 
        yerr = error_din_23_riesgo_whole_one_seq["Alto"])
plt.bar(x + barwidth, medias_din_23_riesgo_whole_one_seq["Bajo"], width = barwidth,
       label = "Bajo", color = colores(70),
       yerr = error_din_23_riesgo_whole_one_seq["Bajo"])
plt.bar(x + (2 * barwidth), medias_din_23_riesgo_whole_one_seq["Indefinido"], width = barwidth,
       label = "Indefinido", color = colores(120),
       yerr = error_din_23_riesgo_whole_one_seq["Indefinido"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.1])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_23_riesgo_whole_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.savefig("./results/selected_din_riesgo_whole_one_seq.png")

plt.show()


# In[205]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 50
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.30)

barwidth = 0.25
x = np.arange(3)
y = np.linspace(3, -0.25, 3)
z = [1] * 3

plt.subplot(131)
plt.title("Género Alfa, genes concatenados 1 secuencia. \nDinucleótidos sin importar la posición", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_riesgo_cat_one_seq["Alto"], width = barwidth, 
       label = "Alto", color = colores(40), 
        yerr = error_din_riesgo_cat_one_seq["Alto"])
plt.bar(x + barwidth, medias_din_riesgo_cat_one_seq["Bajo"], width = barwidth,
       label = "Bajo", color = colores(70),
       yerr = error_din_riesgo_cat_one_seq["Bajo"])
plt.bar(x + (2 * barwidth), medias_din_riesgo_cat_one_seq["Indefinido"], width = barwidth,
       label = "Indefinido", color = colores(120),
       yerr = error_din_riesgo_cat_one_seq["Indefinido"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.8])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_riesgo_cat_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.subplot(132)
plt.title("Género Alfa, genes concatenados 1 secuencia. \nDinucleótidos en posición 1-2", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_12_riesgo_cat_one_seq["Alto"], width = barwidth, 
       label = "Alto", color = colores(40), 
        yerr = error_din_12_riesgo_cat_one_seq["Alto"])
plt.bar(x + barwidth, medias_din_12_riesgo_cat_one_seq["Bajo"], width = barwidth,
       label = "Bajo", color = colores(70),
       yerr = error_din_12_riesgo_cat_one_seq["Bajo"])
plt.bar(x + (2 * barwidth), medias_din_12_riesgo_cat_one_seq["Indefinido"], width = barwidth,
       label = "indefinido", color = colores(120),
       yerr = error_din_12_riesgo_cat_one_seq["Indefinido"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.8])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_12_riesgo_cat_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.subplot(133)
plt.title("Género Alfa, genes concatenados 1 secuencia. \nDinucleótidos en posición 2-3", 
          fontsize = 30, fontweight = "bold")
plt.bar(x, medias_din_23_riesgo_cat_one_seq["Alto"], width = barwidth, 
       label = "Alto", color = colores(40), 
        yerr = error_din_23_riesgo_cat_one_seq["Alto"])
plt.bar(x + barwidth, medias_din_23_riesgo_cat_one_seq["Bajo"], width = barwidth,
       label = "Bajo", color = colores(70),
       yerr = error_din_23_riesgo_cat_one_seq["Bajo"])
plt.bar(x + (2 * barwidth), medias_din_23_riesgo_cat_one_seq["Indefinido"], width = barwidth,
       label = "Indefinido", color = colores(120),
       yerr = error_din_23_riesgo_cat_one_seq["Indefinido"])
plt.plot(y, z, color = "black", linestyle = "dashed")
plt.xlim([-0.25, 2.75])
plt.ylim([0, 1.8])
plt.xticks([r + barwidth for r in range(len(x))], list(medias_din_23_riesgo_cat_one_seq.index), fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylabel("Relación observados/esperados", fontsize = 25, fontweight = "bold")
plt.legend(fontsize = 25)

plt.savefig("./results/selected_din_riesgo_cat_one_seq.png")

plt.show()


# In[22]:


gc_list = ["Total_GC", "GC12", "GC3"]
columnas_extraer_gc = gc_list[:]
columnas_extraer_gc.append("Tipo")


# In[19]:


def select_add_riesgo_gc(df, columnas_extraer):
    
    df_selected = df[columnas_extraer]
    
    list_riesgo = []

    for virus_type in df_selected["Tipo"]:
    
        if virus_type in alto:
        
            list_riesgo.append("Alto")

        elif virus_type in bajo:
        
            list_riesgo.append("Bajo")

        elif virus_type in indefinido:
        
            list_riesgo.append("Indefinido")

        else:
            print("El tipo " + virus_type + " no está clasificado")
    
    df_selected = df_selected.assign(Riesgo = list_riesgo)
    
    return df_selected


# In[17]:


df_whole_gc_one_seq = pd.read_csv("./results/whole_genome_gc_original_one_seq.csv", sep = "\t")
df_cat_gc_one_seq = pd.read_csv("./results/genome_cat_gc_original_one_seq.csv", sep = "\t")


# In[18]:


alfa_whole_gc_one_seq = df_whole_gc_one_seq.loc[df_whole_gc_one_seq["Genero"].isin(["Alfa"])]
alfa_cat_gc_one_seq = df_cat_gc_one_seq.loc[df_cat_gc_one_seq["Genero"].isin(["Alfa"])]


# In[23]:


alfa_whole_gc_one_seq_riesgo = select_add_riesgo_gc(alfa_whole_gc_one_seq, columnas_extraer_gc)
alfa_cat_gc_one_seq_riesgo = select_add_riesgo_gc(alfa_cat_gc_one_seq, columnas_extraer_gc)


# In[11]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.3)

number = 1

plt.suptitle("Género Alfa. Genoma entero. 1 secuencia", fontsize = 35)

for risk in alfa_whole_gc_one_seq_riesgo["Riesgo"].unique():
    
    gc3 = alfa_whole_gc_one_seq_riesgo.loc[alfa_whole_gc_one_seq_riesgo["Riesgo"].isin([risk])]["GC3"]
    gc12 = alfa_whole_gc_one_seq_riesgo.loc[alfa_whole_gc_one_seq_riesgo["Riesgo"].isin([risk])]["GC12"]
        
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(gc3, gc12)

    plt.subplot(1,3, number)
    plt.plot(gc3, gc12, "o")
    plt.plot(gc3, intercept + slope*gc3, "r", label = "Fitted line.\nb =" + str(round(slope, 2)) 
                 +", p = "+ str('{:0.3e}'.format(pvalue)))
    
    plt.xlabel("GC\u2083", fontsize = 25)
    plt.ylabel("GC\u2081"+"\u2082", fontsize = 25)
    #plt.xlim([21, 56])
    #plt.ylim([37, 52])
    plt.xlim([25, 55])
    plt.ylim([30, 55])
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(fontsize = 20)
    plt.title("Riesgo "+risk, fontsize = 30)
    
    number += 1
    
plt.savefig("./results/neutralidad_riesgo_genoma_entero_one_seq.png")

plt.show()


# In[249]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.3)

number = 1
plt.suptitle("Género Alfa. Genes concatenados. 1 secuencia", fontsize = 35)

for risk in alfa_cat_gc_one_seq_riesgo["Riesgo"].unique():
    
    gc3 = alfa_cat_gc_one_seq_riesgo.loc[alfa_cat_gc_one_seq_riesgo["Riesgo"].isin([risk])]["GC3"]
    gc12 = alfa_cat_gc_one_seq_riesgo.loc[alfa_cat_gc_one_seq_riesgo["Riesgo"].isin([risk])]["GC12"]
        
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(gc3, gc12)

    plt.subplot(1,3, number)
    plt.plot(gc3, gc12, "o")
    plt.plot(gc3, intercept + slope*gc3, "r", label = "Fitted line.\nb =" + str(round(slope, 2)) 
                 +", p = "+ str('{:0.3e}'.format(pvalue)))
    
    plt.xlabel("GC\u2083", fontsize = 25)
    plt.ylabel("GC\u2081"+"\u2082", fontsize = 25)
    #plt.xlim([21, 56])
    #plt.ylim([37, 52])
    plt.xlim([20, 60])
    plt.ylim([38, 53])
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(fontsize = 20)
    plt.title("Riesgo "+risk, fontsize = 30)
    
    number += 1
    
plt.savefig("./results/neutralidad_riesgo_genoma_cat_one_seq.png")

plt.show()


# ### Mismas gráficas eliminando el outlier

# In[52]:


alfa_whole_gc_one_seq_riesgo.loc[alfa_whole_gc_one_seq_riesgo["Riesgo"].isin(["Alto"])]


# In[60]:


alfa_whole_gc_one_seq_riesgo.loc[47]


# In[24]:


alfa_cat_gc_one_seq_riesgo.loc[alfa_cat_gc_one_seq_riesgo["Riesgo"].isin(["Alto"])]


# In[62]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.3)

number = 1

plt.suptitle("Género Alfa. Genoma entero. 1 secuencia (sin outlier)", fontsize = 35)

for risk in alfa_whole_gc_one_seq_riesgo["Riesgo"].unique():
    
    gc3 = alfa_whole_gc_one_seq_riesgo.loc[alfa_whole_gc_one_seq_riesgo["Riesgo"].isin([risk])]["GC3"]
    gc12 = alfa_whole_gc_one_seq_riesgo.loc[alfa_whole_gc_one_seq_riesgo["Riesgo"].isin([risk])]["GC12"]
    
    if risk == "Alto": 
        
        for i in range(len(gc12.values) - 1):
        
            if gc12.values[i] < 40 and gc3.values[i] > 43:
                
                gc12.pop(gc12.index[i])
                gc3.pop(gc3.index[i])
                position = (alfa_whole_gc_one_seq_riesgo.loc[alfa_whole_gc_one_seq_riesgo                                                        ["Riesgo"].isin([risk])]["Tipo"].index[i])
                
                virus_type = alfa_whole_gc_one_seq_riesgo.loc[position]["Tipo"]
                print("El outlier es el tipo: ", virus_type)
                
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(gc3, gc12)

    plt.subplot(1,3, number)
    plt.plot(gc3, gc12, "o")
    plt.plot(gc3, intercept + slope*gc3, "r", label = "Fitted line.\nb =" + str(round(slope, 2)) 
                 +", p = "+ str('{:0.3e}'.format(pvalue)))
    
    plt.xlabel("GC\u2083", fontsize = 25)
    plt.ylabel("GC\u2081"+"\u2082", fontsize = 25)
    #plt.xlim([21, 56])
    #plt.ylim([37, 52])
    plt.xlim([25, 55])
    plt.ylim([30, 55])
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(fontsize = 20)
    plt.title("Riesgo "+risk, fontsize = 30)
    
    number += 1
    
plt.savefig("./results/neutralidad_riesgo_genoma_entero_one_seq_sin_outlier.png")

plt.show()


# In[40]:


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 30
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size
plt.subplots_adjust(wspace = 0.3)

number = 1
plt.suptitle("Género Alfa. Genes concatenados. 1 secuencia (sin outlier)", fontsize = 35)

for risk in alfa_cat_gc_one_seq_riesgo["Riesgo"].unique():
    
    gc3 = alfa_cat_gc_one_seq_riesgo.loc[alfa_cat_gc_one_seq_riesgo["Riesgo"].isin([risk])]["GC3"]
    gc12 = alfa_cat_gc_one_seq_riesgo.loc[alfa_cat_gc_one_seq_riesgo["Riesgo"].isin([risk])]["GC12"]
    
    if risk == "Alto": 
        
        for i in range(len(gc12.values) - 1):
        
            if gc12.values[i] < 40 and gc3.values[i] > 43:
                
                gc12.pop(gc12.index[i])
                gc3.pop(gc3.index[i])

    slope, intercept, rvalue, pvalue, stderr = stats.linregress(gc3, gc12)
    plt.subplot(1,3, number)
    plt.plot(gc3, gc12, "o")
    plt.plot(gc3, intercept + slope*gc3, "r", label = "Fitted line.\nb =" + str(round(slope, 2)) 
                 +", p = "+ str('{:0.3e}'.format(pvalue)))
    
    plt.xlabel("GC\u2083", fontsize = 25)
    plt.ylabel("GC\u2081"+"\u2082", fontsize = 25)
    #plt.xlim([21, 56])
    #plt.ylim([37, 52])
    plt.xlim([20, 60])
    plt.ylim([38, 53])
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(fontsize = 20)
    plt.title("Riesgo "+risk, fontsize = 30)
    
    number += 1
    
plt.savefig("./results/neutralidad_riesgo_genoma_cat_one_seq_sin_outlier.png")

plt.show()


# In[251]:


def anovas_riesgo(data_df):
    mod_anova_total = ols("Total_GC ~ Riesgo", data = data_df).fit()
    mod_anova_12 = ols("GC12 ~ Riesgo", data = data_df).fit()
    mod_anova_3 = ols("GC3 ~ Riesgo", data = data_df).fit()
    
    anova_table_total = anova.anova_lm(mod_anova_total, typ = 2)
    anova_table_12 = anova.anova_lm(mod_anova_12, typ = 2)
    anova_table_3 = anova.anova_lm(mod_anova_3, typ = 2)
    
    return anova_table_total, anova_table_12, anova_table_3


# In[252]:


anova_whole_gc_riesgo_total_one_seq, anova_whole_gc_riesgo_12_one_seq, anova_whole_gc_riesgo_3_one_seq = anovas_riesgo(alfa_whole_gc_one_seq_riesgo)
anova_cat_gc_riesgo_total_one_seq, anova_cat_gc_riesgo_12_one_seq, anova_cat_gc_riesgo_3_one_seq = anovas_riesgo(alfa_cat_gc_one_seq_riesgo)


# In[253]:


def post_hoc_riesgo(data_df):
    result_total = (pairwise_tukeyhsd(data_df["Total_GC"], data_df["Riesgo"], alpha = 0.05))
    result_12 = (pairwise_tukeyhsd(data_df["GC12"], data_df["Riesgo"], alpha = 0.05))
    result_3 = (pairwise_tukeyhsd(data_df["GC3"], data_df["Riesgo"], alpha = 0.05))
    
    return result_total, result_12, result_3


# In[254]:


def post_hoc_results_table(tukey_res): 
    p_value = psturng(np.abs(tukey_res.meandiffs / tukey_res.std_pairs),
           len(tukey_res.groupsunique), tukey_res.df_total)
    p_val_col = pd.DataFrame(p_value)
        
    table_pd = pd.DataFrame(data = tukey_res._results_table.data[1:],
                               columns = tukey_res._results_table.data[0])
        
    table_pd["P-value"] = p_val_col
    
    return table_pd


# In[256]:


tukey_whole_gc_riesgo_total_one_seq, tukey_whole_gc_riesgo_12_one_seq, tukey_whole_gc_riesgo_3_one_seq = post_hoc_riesgo(alfa_whole_gc_one_seq_riesgo)
tukey_table_whole_gc_riesgo_total_one_seq = post_hoc_results_table(tukey_whole_gc_riesgo_total_one_seq)
tukey_table_whole_gc_riesgo_12_one_seq = post_hoc_results_table(tukey_whole_gc_riesgo_12_one_seq)
tukey_table_whole_gc_riesgo_3_one_seq = post_hoc_results_table(tukey_whole_gc_riesgo_3_one_seq)
tukey_cat_gc_riesgo_total_one_seq, tukey_cat_gc_riesgo_12_one_seq, tukey_cat_gc_riesgo_3_one_seq = post_hoc_riesgo(alfa_cat_gc_one_seq_riesgo)
tukey_table_cat_gc_riesgo_total_one_seq = post_hoc_results_table(tukey_cat_gc_riesgo_total_one_seq)
tukey_table_cat_gc_riesgo_12_one_seq = post_hoc_results_table(tukey_cat_gc_riesgo_12_one_seq)
tukey_table_cat_gc_riesgo_3_one_seq = post_hoc_results_table(tukey_cat_gc_riesgo_3_one_seq)


# In[257]:


with open("./results/anova_tukey_gc_riesgo_results_one_seq", "w") as file:
    
    file.write("Genoma completo. 1 secuencia \n")
    file.write("\nG+C total\n")
    file.write((anova_whole_gc_riesgo_total_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_whole_gc_riesgo_total_one_seq).to_string())
    
    file.write("\n\nGC 1-2\n")
    file.write((anova_whole_gc_riesgo_12_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_whole_gc_riesgo_12_one_seq).to_string())
    
    file.write("\n\nGC3\n")
    file.write((anova_whole_gc_riesgo_3_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_whole_gc_riesgo_3_one_seq).to_string())
    
    file.write("\n\n\nGenes concatenados. 1 secuencia\n")
    file.write("\nG+C total\n")
    file.write((anova_cat_gc_riesgo_total_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_cat_gc_riesgo_total_one_seq).to_string())
    
    file.write("\n\nGC 1-2\n")
    file.write((anova_cat_gc_riesgo_12_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_cat_gc_riesgo_12_one_seq).to_string())

    file.write("\n\nGC3\n")
    file.write((anova_cat_gc_riesgo_3_one_seq).to_string())
    file.write("\n\n")
    file.write((tukey_table_cat_gc_riesgo_3_one_seq).to_string())
    file.write("\n")


# In[ ]:




