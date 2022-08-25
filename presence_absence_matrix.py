#!/usr/bin/env python3

# This script is designed to parse through the 
# Zoonomia BLAST results to extract useful 
# information

#================================================
'''Import modules'''

import os 
from collections import namedtuple
import pandas as pd

#================================================
'''Find BLAST output files'''

working_dir = '/Users/baileyfrancis/OneDrive/Documents/UoN/miRNA_project/Zoonomia_BLAST'
BLAST_files = os.listdir(working_dir)

#================================================
'''Find sequence IDs from Zoonomia file'''

Zoonomia_file = open('/Users/baileyfrancis/OneDrive/Documents/UoN/miRNA_project/data/seqids.txt', 'rt') #Obtained via halStats --branches

#Seperate sequences from Zoonomia file 
for seq in Zoonomia_file:
    Zoonomia_seqs = seq.split()

Zoonomia_IDs = []

for ID in Zoonomia_seqs:
    if ID.startswith('fullTreeAnc'): #Exclude all ancestral sequences from Zoonomia file 
        pass
    else: 
        Zoonomia_IDs.append(ID)

#================================================
'''Create output dataframe'''

df = pd.DataFrame(columns=['Species', 'Mir-127', 'Mir-185',
                           'Mir-188', 'Mir-28', 'Mir-324',
                           'Mir-378', 'Mir-433', 'Mir-331',
                           'Mir-340', 'Mir-505', 'Mir-542',
                           'Mir-423', 'Mir-671'])

#================================================
'''Create a dictionary for each species that will form each row in dataframe'''
species_dicts = []

for species in Zoonomia_IDs:
    exec(f"{species} = {{ 'Species': '{species}', 'Mir-127':0, 'Mir-185':0, 'Mir-188-P1':0, 'Mir-188-P2':0, 'Mir-188-P3':0, 'Mir-28-P1':0, 'Mir-28-P2':0, 'Mir-28-P3':0, 'Mir-324':0, 'Mir-378':0, 'Mir-433':0, 'Mir-331':0, 'Mir-340':0, 'Mir-505':0, 'Mir-542':0, 'Mir-423':0, 'Mir-671':0}}")
    

#================================================
'''Extract species and miRNA information from file'''

BLAST = namedtuple('BLAST', ['qseqid', 'sseqid', 'qlen', 'length', 'pident', 
                             'sstart', 'send', 'evalue']) #named tuple to store BLAST fields 

duplication_files = []

for file in BLAST_files:
    if file.endswith('_BLAST.txt'): #Find BLAST output files 
        
        #Three word species names that previously caused errors in code:
        Canis_lupus_familiaris_str = 'Canis_lupus_familiaris' 
        Ceratotherium_simum_cottoni_str = 'Ceratotherium_simum_cottoni'
        
        if (Canis_lupus_familiaris_str in file) or (Ceratotherium_simum_cottoni_str in file): #If species name has three words...
            species = '_'.join(file.split('_',3)[:3]) #Extract species information from file name
            mirna = '_'.join(file.split('_',3)[3:]) #Extract mirna information from file name
            mirna = mirna.replace('_pre_BLAST.txt', '') #Remove file extension
            mirna = mirna[4:]
        
        else: #If only two word species name...
            species = '_'.join(file.split('_',2)[:2]) #Extract species information from file name
            mirna = '_'.join(file.split('_',2)[2:]) #Extract mirna information from file name
            mirna = mirna.replace('_pre_BLAST.txt', '') #Remove file extension
            mirna = mirna[4:]
            print(mirna)
    
#================================================
        '''Calculate the number of significant hits for each miRNA'''
    
        #Check if file has no significant hits
        file_path = os.path.join(working_dir, file)
        filesize = os.path.getsize(file_path)
        
        if filesize == 0: #If file is empty...
            pass
        
        else: #If file contains hits...
            with open(file_path, 'rt') as file:
                lines = file.readlines() #Create a list containing each BLAST hit 
                
                count = 0 #Counter to track the number of BLAST hits
                
                #Lists to store the values of BLAST fields for each BLAST hit:
                file_pidents = [] 
                file_paligns = []
                file_evalues = []              
                
                for line in lines: #For each hit...
                    BLAST_fields = line.split() #Create a list containing each BLAST field 
                    blast = BLAST._make(BLAST_fields) #Create BLAST named tuple using BLAST field list (see above)
                    
                    palign = (int(blast.length)/int(blast.qlen))*100 #Calculate percentage alignment cover 
                    
                    count += 1 #Count the number of hits 
                    
                    #Add BLAST fields to appropriate list:
                    file_pidents.append(blast.pident)
                    file_paligns.append(palign)
                    file_evalues.append(blast.evalue)
                    
                    if (float(blast.pident) > 90) and (palign > 85): #If a hit satisfies our homolog criteria...
                        exec(f'{species}["{mirna}"] = 1') #Show the miRNA as being present
                    
                    if count > 1: #If there are multiple hits...
                    
                        #If below are true, then BLAST fields are identical for all hits in a given file:
                        identical_pidents = file_pidents.count(file_pidents[0]) == len(file_pidents)
                        identical_paligns = file_paligns.count(file_paligns[0]) == len(file_paligns)
                        identical_evalues = file_evalues.count(file_evalues[0]) == len (file_evalues)
                        
                        #if (identical_pidents) and (identical_paligns) and (identical_evalues): #If BLAST hits identical...
                            #print(f'{mirna} duplication likely an artefact of assembly in {species}')
                        
                        #else: #If BLAST hits differ...
                            #print(f'Possible duplication of {mirna} in {species}')
                            #duplication_files.append(file_path)
                            #break
                        

#================================================
'''Create miRNA presence/absence data frame'''
                       
for species in Zoonomia_IDs:
    exec(f'species_dicts.append({species})') #Create data frame using species dictionaries 
    
df = pd.DataFrame(species_dicts)
df.to_csv('BLAST_hits_summary_Sus.csv') #Convert data frame to csv 

#================================================
'''Create a file containing a list of potential duplications '''

duplication_output = open('potential_duplications.txt', 'w')

for file in duplication_files:
    duplication_output.write(f'{file}\n')
    
duplication_output.close()


