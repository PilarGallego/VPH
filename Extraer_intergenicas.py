#!/usr/bin/env python
# coding: utf-8


import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
import re
import logging
from Bio.SeqUtils import GC, GC123
from Bio.Alphabet import generic_dna
from collections import Counter
from Bio.Seq import MutableSeq


logging.basicConfig(filename = "avisos_integenicas.log", filemode = "w", 
                    format = "%(name)s - %(levelname)s - %(asctime)s - %(message)s",
                    level = logging.INFO)

list_files = []
for root, directories, files in os.walk("./GenomasHPV/"):
    for name in files:
        list_files.append(os.path.join(root, name))

genes = ['E1', 'E2', 'E4', 'E5', 'E6', 'E7', 'L2', 'L1']


intergenic_genome = defaultdict(list)
intergenic_genome_info = defaultdict(list)

for file in list_files:
    with open(file, "r") as f_hpv:
        data = SeqIO.parse(f_hpv, "genbank")
        
        virus_genus = f_hpv.name.split("/")[-2]
        virus_type = (re.findall(r'\d+', f_hpv.name.split("/")[-1])[0])
        
        for seqrecord in data:
            
            if virus_genus in ["Alfa", "Beta", "Gamma"]:# and virus_type not in intergenic_genome_info[virus_genus].keys():
                intergenic_region = []
                intergenic_record_info = []
                intergenic_record = ""
                intergenic_genome_type = defaultdict(list)
                intergenic_genome_info_type = defaultdict(list)
                
                #print(virus_genus, ";", virus_type)
        
                #gen_present = {"E1" : False, "E2" : False, "E4" : False, "E5" : False, "E6" : False, "E7" : False,
                 #         "L1" : False, "L2" : False}
                gen_present = {"E1" : False, "E2" : False, "E6" : False, "E7" : False,
                               "L1" : False, "L2" : False}
                
                if seqrecord.features:
                
                    gene = ""
                    
                    for feature in seqrecord.features:
                            
                        if feature.type == "source":
                            
                            inicio = feature.location._start
                            final = feature.location._end
                            
                        elif feature.type == "CDS":
                        
                            if "gene" in feature.qualifiers:
                            
                                if feature.qualifiers["gene"][0].split()[0] in gen_present.keys():
                                    cds = feature.qualifiers["gene"][0].split()[0]
                                    cds = cds.upper()
                                                                
                                    if not gen_present[cds]:
                                        
                                        try:
                                    
                                            feature.location._start
                                            gen_present[cds] = True

                                            start = feature.location._start#.location
                                            end = feature.location._end#.position
                                        
                                            if feature.strand == 1:
                                                intergenic_region.append((start, end, 1))    
                                        except:
                                            aviso = (virus_genus + virus_type + ": sin start")
                            elif "product" in feature.qualifiers:
                            
                                if feature.qualifiers["product"][0].split()[0] in gen_present.keys():
                                    cds = feature.qualifiers["product"][0].split()[0]
                                    cds = cds.upper()
                                
                                    #print(cds)
                                    if not gen_present[cds]:
                                    
                                        gen_present[cds] = True
                                        
                                        start = feature.location._start#.location
                                        end = feature.location._end#.position
                                        if feature.strand == 1:
                                            intergenic_region.append((start, end, 1))    
                    intergenic_region.append([end, final, 1])
                    
                    for i, pospair in enumerate(intergenic_region):
                                  
                        if i == 0:
                            
                            last_end = inicio
                            this_start = intergenic_region[i][0]
                            strand = pospair[2]
                            #print(last_end, this_start)
                            
                        elif i > 0 and i < (len(intergenic_region) - 1):
                            
                            last_end = intergenic_region[i-1][1]
                            this_start = pospair[0]
                            strand = pospair[2]
                            #print(last_end, this_start)
                                                
                        elif i == (len(intergenic_region) - 1):
                            last_end = intergenic_region[i][0]
                            this_start = final
                            strand = pospair[2]
                            #print(last_end, this_start)
                            
                        if (this_start - last_end) >= 1:
                        
                            intergenic_region_seq = seqrecord.seq[last_end:this_start]
                            strand_string = "+"
                            intergenic_record_info.append(SeqRecord(intergenic_region_seq, id = "%s-interg-%s" % (seqrecord.name, i),
                                        description = "%s %d-%d %s" % (seqrecord.name, last_end + 1, this_start, "+")))
                            intergenic_record += intergenic_region_seq 
                            
                    #print(len(intergenic_record))
                    gen_present = {"E1" : False, "E2" : False, "E6" : False, "E7" : False,
                               "L1" : False, "L2" : False}
                
                    if len(intergenic_record) > 0 :
                        intergenic_genome_info_type[virus_type].append(intergenic_record_info)
                        intergenic_genome_type[virus_type].append(intergenic_record)
            
                intergenic_genome[virus_genus].append(intergenic_genome_type)
                intergenic_genome_info[virus_genus].append(intergenic_genome_info_type)

