# STEP 1: IMPORT PACKAGES
import pandas as pd
import numpy as np
import pinetree as pt
from sklearn.metrics import mean_squared_error
import sys
from datetime import date

today = date.today()

# @TODO Update basedir file path
# @TODO Update pinetree output file path
# @TODO Update genomic_coords file path. This should like to the CSV with the gene locations
# @TODO Update promoter & terminator values --> optimized values included in the 3152022 analysis
# @TODO Note that this still has the tH location at the very end. This needs to be moved earlier so pA, tH, then gA..

# STEP 2: RUN PINETREE FUNCTION
def main(pA, pB, pD, tJ, tF, tG, tH, date_):
    print(pA, pB, pD, tJ, tF, tG, tH)
    print("Defining PhiX-174 genome")

    # Create host cell & genome
    CELL_VOLUME = (1.1e-15)  # from T7
    model = pt.Model(cell_volume=CELL_VOLUME)
    phage = pt.Genome(name="phix_174", length=5395)

    # Read genomic coordinates from csv into dataframe
    genomic_coords = pd.read_csv("/stor/work/Wilke/tingle/phix174/output/genomic_coords.csv")
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
    phage.add_terminator(name="terminator_J", start=2436, stop=2437,  # Right before gene F start=2404, stop=3687,
                         efficiency={"ecolipol": 0.7})  # 0.7
    phage.add_terminator(name="terminator_F", start=3796, stop=3797,  # Right before gene G start=3798, stop=4325
                         efficiency={"ecolipol": 0.8})  # 0.8
    #phage.add_terminator(name="terminator_G", start=4332, stop=4333,  # Right before gene H start=4334, stop=5320
    #                     efficiency={"ecolipol": tG})  # 0.6
    phage.add_terminator(name="terminator_G", start=4374, stop=4375,  # Right before gene H start=4334, stop=5320
                         efficiency={"ecolipol": 0.6})  # 4374 
    #phage.add_terminator(name="terminator_H", start=5321, stop=5322,  # Right after gene H
    #                     efficiency={"ecolipol": tH})  # 5321
    phage.add_terminator(name="terminator_H", start=5390, stop=5391,  # Right after gene H
                         efficiency={"ecolipol": 0.9})  # 5321

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
    model.add_ribosome(10, 30, 100)
    model.add_species("bound_ribosome", 100)
    model.seed(34)

    # Run simulation
    print("running simulation")
    model.simulate(time_limit=1200, time_step=5,
                   output=f"{date_}.tsv")
    print("Simulation successful!")

if __name__ == '__main__':
    date_ = today.strftime("%b-%d-%Y")
    #base_dir = "."
    #folder = "fix_geneH"
    #test = "increase_ribosomes"
    # promoter values from 3152022 run, sim 7 gen 260
    pA = 13.04585
    pB = 15.222113
    pD = 18.342678
    tJ = 0.01
    tF = 0.01
    tG = 0.01
    tH = 0.001
    main(pA, pB, pD, tJ, tF, tG, tH, date_)

