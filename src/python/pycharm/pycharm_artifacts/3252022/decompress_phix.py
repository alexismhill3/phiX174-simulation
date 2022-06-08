# STEP 1: IMPORT PACKAGES
import pandas as pd
import numpy as np
import pinetree as pt
from sklearn.metrics import mean_squared_error
import sys

# STEP 2: RUN PINETREE FUNCTION
def main(pA, pB, pD, tJ, tF, tG, tH, base_dir, date):
    print(pA, pB, pD, tJ, tF, tG, tH)
    print("Defining PhiX-174 genome")

    # Create host cell & genome
    CELL_VOLUME = (1.1e-15)  # from T7
    model = pt.Model(cell_volume=CELL_VOLUME)
    phage = pt.Genome(name="phix_174", length=5400) # 5386 old genome length

    # Read genomic coordinates from csv into dataframe
    genomic_coords = pd.read_csv(base_dir + "src/python/pycharm/3252022/" + "decompress_genomic_coords.csv")
    
    # Change locations of gene G, terminator G, and gene H to have more spacing between them
     
    print(genomic_coords.at[0, "type"])

    n = 0
    while (n < len(genomic_coords)):

        if genomic_coords.at[n, "type"] == "gene":
            phage.add_gene(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                           start=genomic_coords.at[n, "new_start" ],
                           stop=genomic_coords.at[n, "new_end"],
                           rbs_start=genomic_coords.at[n, "new_start"] - 15,
                           rbs_stop=genomic_coords.at[n, "new_start"] - 1, rbs_strength=1e7) # change fro 1e7

        elif genomic_coords.at[n, "type"] == "promoter" and genomic_coords.at[n, "name"] == "A":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               interactions={"ecolipol": np.exp(pA)})

        elif genomic_coords.at[n, "type"] == "promoter" and genomic_coords.at[n, "name"] == "B1":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               interactions={"ecolipol": np.exp(pB)})

        elif genomic_coords.at[n, "type"] == "promoter" and genomic_coords.at[n, "name"] == "D":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               interactions={"ecolipol": np.exp(pD)})
                               
        elif genomic_coords.at[n, "type"] == "terminator" and genomic_coords.at[n, "name"] == "H":
            phage.add_terminator(name="terminator_H", 
                                start=genomic_coords.at[n, "new_start"], 
                                stop=genomic_coords.at[n, "new_end"],  
                                efficiency={"ecolipol": tH})
        
        elif genomic_coords.at[n, "type"] == "terminator" and genomic_coords.at[n, "name"] == "J":
            phage.add_terminator(name="terminator_J", 
                                start=genomic_coords.at[n, "new_start"], 
                                stop=genomic_coords.at[n, "new_end"],  
                                efficiency={"ecolipol": tJ})
            
        elif genomic_coords.at[n, "type"] == "terminator" and genomic_coords.at[n, "name"] == "F":
            phage.add_terminator(name="terminator_F", 
                                start=genomic_coords.at[n, "new_start"], 
                                stop=genomic_coords.at[n, "new_end"],  
                                efficiency={"ecolipol": tF})
                                
        elif genomic_coords.at[n, "type"] == "terminator" and genomic_coords.at[n, "name"] == "G":
            phage.add_terminator(name="terminator_G", 
                                start=genomic_coords.at[n, "new_start"], 
                                stop=genomic_coords.at[n, "new_end"],  
                                efficiency={"ecolipol": tG})

        n = n + 1

    print("all genes, promoters, and terminators added")

    # Register genome after promoters/terminators are added
    model.register_genome(phage)
    print("genome is registered")

    # Define interactions
    print("Defining Polymerases & Interactions")
    # Add polymerases & species
    model.add_polymerase(name="ecolipol", speed=35, footprint=35, copy_number=0)
    model.add_species("bound_ecolipol", 2000)  
    model.add_species("ecoli_genome", 0)
    model.add_species("ecoli_transcript", 0)
    model.add_reaction(1e7, ["ecolipol", "ecoli_genome"], ["bound_ecolipol"])  # 1e6
    model.add_reaction(0.04, ["bound_ecolipol"], ["ecolipol", "ecoli_genome", "ecoli_transcript"])
    model.add_ribosome(10, 30, 100) # Number of ribosomes changed from 100
    model.add_species("bound_ribosome", 100) # Number of ribosomes changed from 100
    model.seed(34)

    # Run simulation
    print("running simulation")
    model.simulate(time_limit=1200, time_step=5,output=base_dir + "src/python/pycharm/" + str(date) + "/" + str(test) + ".tsv")
    print("Simulation successful!")


if __name__ == '__main__':
    date = 3252022
    base_dir = "/stor/work/Wilke/tingle/phix174/"
    pA = float(sys.argv[1])
    pB = float(sys.argv[2])
    pD = float(sys.argv[3])
    tJ = float(sys.argv[4])
    tF = float(sys.argv[5])
    tG = float(sys.argv[6])
    tH = float(sys.argv[7])
    test = sys.argv[8]
    main(pA, pB, pD, tJ, tF, tG, tH, base_dir, date)
    
    
    
