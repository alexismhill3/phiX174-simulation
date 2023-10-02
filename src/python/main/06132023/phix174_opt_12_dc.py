# EMPTY!! - FINALIZE WT SCRIPT & COPY HERE. MODIFY GENOMIC_COORDS
# OBJECTIVE: This script will run the pinetree model for a decompressed genome. 
# REQUIREMENTS: User must specify the promoter and terminator strengths. Also need to specify the 
# decompressed genome coordinates. 
# OUTPUT: This script will write results to the working directory/output/decompressed
#----------------------------------------------------------------------------------------------------------------#

# STEP 1: IMPORT PACKAGES
import pandas as pd
import numpy as np
import pinetree as pt
import sys

# STEP 2: RUN PINETREE FUNCTION
def main(sim, gen, pA, pB, pD, tJ, tF, tG, tH, base_dir, genomic_coords_path, date, seed):
    
    # STEP 2.1: Create E. Coli host cell and phage genome
    CELL_VOLUME = 1.1e-15  
    model = pt.Model(cell_volume=CELL_VOLUME)
    phage = pt.Genome(name="phix_174", length=5386)

    # STEP 2.2: Read genomic coordinates from csv into a dataframe
    genomic_coords = pd.read_csv(base_dir + genomic_coords_path)
    
    # STEP 2.3: Add genomic features from dataframe to the genome
    n = 0
    while (n < len(genomic_coords)):

        if genomic_coords.at[n, "type"] == "gene":
            phage.add_gene(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                           start=genomic_coords.at[n, "new_start"],
                           stop=genomic_coords.at[n, "new_end"],
                           rbs_start=genomic_coords.at[n, "new_start"] - 15,
                           rbs_stop=genomic_coords.at[n, "new_start"] - 1, rbs_strength=1e7)

        elif genomic_coords.at[n, "type"] == "promoter" and genomic_coords.at[n, "name"] == "A":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               interactions={"ecolipol": np.multiply(pA, 1000000)})

        elif genomic_coords.at[n, "type"] == "promoter" and genomic_coords.at[n, "name"] == "B1":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               interactions={"ecolipol": np.multiply(pB, 1000000)})

        elif genomic_coords.at[n, "type"] == "promoter" and genomic_coords.at[n, "name"] == "D":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               interactions={"ecolipol": np.multiply(pD, 1000000)})
        
        elif genomic_coords.at[n, "type"] == "terminator" and genomic_coords.at[n, "name"] == "H":
            phage.add_terminator(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               efficiency={"ecolipol": tH})
                               
        elif genomic_coords.at[n, "type"] == "terminator" and genomic_coords.at[n, "name"] == "J":
            phage.add_terminator(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               efficiency={"ecolipol": tJ})
                               
        elif genomic_coords.at[n, "type"] == "terminator" and genomic_coords.at[n, "name"] == "F":
            phage.add_terminator(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               efficiency={"ecolipol": tF})
                               
        elif genomic_coords.at[n, "type"] == "terminator" and genomic_coords.at[n, "name"] == "G":
            phage.add_terminator(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               efficiency={"ecolipol": tG})
                               
        else:
            print(" ")

        n = n + 1

    # Register genome after promoters/terminators are added
    model.register_genome(phage)

    # Add polymerases & species
    #model.add_polymerase(name="ecolipol", speed=35, footprint=35, copy_number=0)
    model.add_polymerase_with_readthrough(name="ecolipol", speed=35, footprint=35, copy_number=0)

    model.add_species("bound_ecolipol", 1800) 
    model.add_species("ecoli_genome", 0)
    model.add_species("ecoli_transcript", 0)
    
    # Reactions
    model.add_reaction(1e7, ["ecolipol", "ecoli_genome"], ["bound_ecolipol"])  
    model.add_reaction(0.04, ["bound_ecolipol"], ["ecolipol", "ecoli_genome", "ecoli_transcript"])
    model.add_reaction(1e6, ["ecoli_transcript", "__ribosome"], ["bound_ribosome"])
    model.add_reaction(0.04, ["bound_ribosome"], ["__ribosome", "ecoli_transcript"])
    model.add_reaction(0.001925, ["ecoli_transcript"], ["degraded_transcript"])
    
    # Add ribosomes
    model.add_ribosome(10, 30, 30)
    model.add_species("bound_ribosome", 10000)
    model.seed(seed)

    # Run simulation
    model.simulate(time_limit=500, time_step=5,output=base_dir + "output/" + str(date) + "/" + str(sim) + "_gen_" + str(gen) + "_seed_" + str(seed) + ".tsv")

if __name__ == '__main__':
    sim = sys.argv[1]
    gen = sys.argv[2]
    pA = float(sys.argv[3])
    pB = float(sys.argv[4])
    pD = float(sys.argv[5])
    tJ = float(sys.argv[6])
    tF = float(sys.argv[7])
    tG = float(sys.argv[8])
    tH = float(sys.argv[9])
    base_dir = sys.argv[10]
    genomic_coords_path = sys.argv[11]
    date = sys.argv[12]
    seed = int(sys.argv[13])
    main(sim, gen, pA, pB, pD, tJ, tF, tG, tH, base_dir, genomic_coords_path, date, seed)

