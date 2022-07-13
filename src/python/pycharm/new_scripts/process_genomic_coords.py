# Objective: This script extracts relevant genomic elements (gene locations, promoter sites)
# from an ENSEMBLE file for PhiX-174. It then calculates the 'linearized' coordinates.
# The terminator positions (which are relative to the positions of the linearized gene coordinates)
# are hardcoded. The sequence of relevant elements is as follows: 
## promoter A, terminator H, gene A, gene A*, promoter B1, gene B, 
## gene K, gene C, promoter D, gene D, gene E, gene J,
## terminator J, gene F, terminator F, gene G, terminator G. 
### NOTE: we exclude promoter B2, because it was not referenced in Logel 2020. 
# An output file called 'genomidddddc_coords.csv' is written to the working directory. 
# This file is read into the optimization script.  


# STEP 1: Import packages
from Bio import Entrez, SeqIO
import pandas as pd

# STEP 2: Parse genomic coordinates from Genebank File
Entrez.email = "tanvi.ingle@utexas.edu"
handle = Entrez.efetch(db="nuccore", id=["MN385565.1"], rettype="gb",retmode="text")
record = SeqIO.read(handle, "genbank")

df = pd.DataFrame(columns=['type', 'name', 'start1', \
                           'end1','start2','end2'])

for feature in record.features:
    name = ""
    feature_type = feature.type
    start1 = None
    end1 = None
    start2 = None
    end2 = None
    if feature_type == 'gene':
        name = feature.qualifiers.get('gene')[0]
        if feature.location_operator == 'join':
            first = feature.location.parts[0]
            second = feature.location.parts[1]
            start1 = first.start
            end1 = first.end
            start2 = second.start
            end2 = second.end
        else:
            start1 = feature.location.start
            end1 = feature.location.end

        new_row = {'type':feature_type, 'name':name, \
                   'start1':start1, 'end1':end1, \
                   'start2': start2, 'end2':end2}
        df = df.append(new_row, ignore_index=True)

    # we assume that regulatory (promoter) feature types do not
    # have CompoundLocation ``location'' fields
    if feature_type == 'regulatory':
        name = feature.qualifiers.get('note')[0]
        # remove Promoter from the name
        name = name.replace('Promoter ', '')
        if feature.location_operator == 'join':
            print('join not implemented for regulatory ',
                  'feature type!')
            continue
        else:
            start1 = feature.location.start
            end1 = feature.location.end

        new_row = {'type':'promoter', 'name':name, \
                   'start1':start1, 'end1':end1, \
                   'start2': start2, 'end2':end2}
        # append row to the dataframe (check efficiency)
        df = df.append(new_row, ignore_index=True)


# STEP 3: Calculate Lineraized genomic positions
i=0
df["new_start"] = ""
df["new_end"] = ""

index = df['start1'][12]
end_index =  df['end1'][0]

while i < len(df):
    row = df.loc[i]
    if(row['type'] == "promoter" and row['name'] == "A"):
        #print("pA")
        df.at[i,'new_start'] = 1
        df.at[i, 'new_end'] = row['end1'] - row['start1'] + 1
        
    elif(row['type'] == "promoter" and row['name'] in ["B1","B2"]):
        #print("Displaced promoters")
        df.at[i,'new_start'] = (row['start1'] - index)
        df.at[i, 'new_end'] = row['end1'] - row['start1'] + df.at[i,'new_start']
        
    elif(row['type'] == "gene" and row['name'] in ["A", "A*", "B"]):
        #print("looping genes")
        df.at[i,'new_start'] = (row['start1'] - index)
        df.at[i, 'new_end'] = (row['end1'] - row['start1']) + (row['end2'] - row['start2']) + df.at[i,'new_start']
    
    else:
        #print("remaining genes")
        df.at[i,'new_start'] = row['start1'] + (end_index-index)
        df.at[i, 'new_end'] = (row['end1'] - row['start1']) + df.at[i,'new_start'] 

    i = i+1

# STEP 4: Add terminators
# Positions are relative to the end of their corresponding genes.
# Position of Terminator G has been shifted further along, to avoid overlap with gene H
add_tH = {'type': "terminator", 
          'name' : "H", 
          'start1' : "NA",
          'end1' : "NA", 
          'start2' : "NA",
          'end2' : "NA",
          'new_start' : 55,
          'new_end' : 56}

add_tJ = { 'type': "terminator", 
          'name' : "J", 
          'start1' : "NA",
          'end1' : "NA", 
          'start2' : "NA",
          'end2' : "NA",
          'new_start' : 2436,
          'new_end' : 2437}
          
add_tF = {'type': "terminator", 
          'name' : "F", 
          'start1' : "NA",
          'end1' : "NA", 
          'start2' : "NA",
          'end2' : "NA",
          'new_start' : 3796,
          'new_end' : 3797}
          
add_tG = {'type': "terminator", 
          'name' : "G", 
          'start1' : "NA",
          'end1' : "NA", 
          'start2' : "NA",
          'end2' : "NA",
          'new_start' : 4400,
          'new_end' : 4401}
          
df = df.append(add_tH, ignore_index = True)
df = df.append(add_tJ, ignore_index = True)
df = df.append(add_tF, ignore_index = True)
df = df.append(add_tG, ignore_index = True)

# STEP 5: Write file to working directory
df.to_csv("genomic_coords.csv")















