import pandas as pd
import numpy as np
import pinetree as pt
import sys

def run_pt(sim, gen, pA, pB, pD, tJ, tF, tG, tH, base_dir, output_folder):
    print(gen, pA, pB, pD, tJ, tF, tG, tH)
    print("Defining PhiX-174 genome")

    # Create host cell & genome
    CELL_VOLUME = 1.1e-15  # from T7
    model = pt.Model(cell_volume=CELL_VOLUME)
    phage = pt.Genome(name="phix_174", length=5400)

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

        else:
            print("ignoring pB2")

        n = n + 1

    print("all genes and promoters added")

    # Add terminators manually (moved to intragenic regions, tH at start, and genome lengthend)
    phage.add_terminator(name="terminator_J", start=2436, stop=2437,  
                         efficiency={"ecolipol": tJ}) # After gene J
    phage.add_terminator(name="terminator_F", start=3796, stop=3797,  
                         efficiency={"ecolipol": tF})  # After gene F
    phage.add_terminator(name="terminator_G", start=4400, stop=4401,
                       efficiency={"ecolipol": tG}) # After gene G
    phage.add_terminator(name="terminator_H", start=55, stop=56,
                 efficiency={"ecolipol": tH})

    print("all terminators added")

    # Register genome after promoters/terminators are added
    model.register_genome(phage)
    print("genome is registered")

    # Define interactions
    print("Defining Polymerases & Interactions")
    # Add polymerases & species
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
    #model.seed(34)
    print("Not using seeds")

    # Run simulation
    print("running simulation")
    model.simulate(time_limit=500, time_step=5,output=base_dir + "output/" + output_folder + "/sim_" + sim + "_gen_" + str(gen) + ".tsv")
    print("Simulation successful!")


if __name__ == "__main__":
    gen = sys.argv[1]
    base_dir = sys.argv[2]
    output_folder = sys.argv[3]
    sim = "11072022_13"
    pA = 11.610219
    pB = 4.728556
    pD = 20.615451
    tJ = 0.3533596
    tF = 0.247685719
    tG = 0.5370027
    tH = 0.9807020
    run_pt(sim, gen, pA, pB, pD, tJ, tF, tG, tH, base_dir, output_folder)
