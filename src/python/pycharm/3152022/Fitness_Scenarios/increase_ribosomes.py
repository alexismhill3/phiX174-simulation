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
    genomic_coords = pd.read_csv(base_dir + "output/" + "genomic_coords.csv")
    
    # Change locations of gene G, terminator G, and gene H to have more spacing between them
        # gH ends at 5390, starts at 4413
        # tG 4400 - 4401    
    genomic_coords.at[11, "new_start"] = 4413
    genomic_coords.at[11, "new_end"] = 5390
    
    print(genomic_coords.at[0, "type"])

    n = 0
    while (n < len(genomic_coords)):

        if genomic_coords.at[n, "type"] == "gene":
            phage.add_gene(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                           start=genomic_coords.at[n, "new_start" ],
                           stop=genomic_coords.at[n, "new_end"],
                           rbs_start=genomic_coords.at[n, "new_start"] - 15,
                           rbs_stop=genomic_coords.at[n, "new_start"] - 1, rbs_strength=1e10) # change fro 1e7

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

        else:
            print("ignoring pB2")

        n = n + 1

    print("all genes and promoters added")

    # Add terminators manually
    # phage.add_terminator(name="terminator_J", start=2434, stop=2435,  
    #                      efficiency={"ecolipol": tJ}) # After gene J, which ends at 2433
    # phage.add_terminator(name="terminator_F", start=3754, stop=3755,  
    #                      efficiency={"ecolipol": tF})  # After gene F, which ends at 3753
    # phage.add_terminator(name="terminator_G", start=4392, stop=4393,  
    #                      efficiency={"ecolipol": tG}) # After gene G, which ends at 4391
    # phage.add_terminator(name="terminator_H", start=46, stop=47,  
    #              efficiency={"ecolipol": tH})  # After pA, which ends at 45


    phage.add_terminator(name="terminator_J", start=2436, stop=2437,  
                         efficiency={"ecolipol": tJ}) # After gene J, which ends at 2433
    phage.add_terminator(name="terminator_F", start=3796, stop=3797,  
                         efficiency={"ecolipol": tF})  # After gene F, which ends at 3753
    phage.add_terminator(name="terminator_G", start=4400, stop=4401,
                       efficiency={"ecolipol": tG}) # After gene G, which ends at 4391
    phage.add_terminator(name="terminator_H", start=55, stop=56,
                 efficiency={"ecolipol": tH})  # After pA, which ends at 45
    #phage.add_terminator(name="terminator_X", start=5387, stop=5388,
     #            efficiency={"ecolipol": tH}) 
   
    print("all terminators added")

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
    model.simulate(time_limit=1200, time_step=5,output=base_dir + "output/" + str(date) + "/" + str(folder)+ "/" + str(test) + ".tsv")
    print("Simulation successful!")


if __name__ == '__main__':
    date = 3152022
    base_dir = "/stor/work/Wilke/tingle/phix174/"
    folder = "fix_geneH"
    test = "increase_ribosomes"
    # Values from 3152022 run, sim 7 gen 260
    pA = 13.04585
    pB = 15.222113
    pD = 18.342678
    tJ = 0.3129013
    tF = 0.6671101
    tG = 0.4105418
    tH = 0.7791773
    main(pA, pB, pD, tJ, tF, tG, tH, base_dir, date)
    
    
    
