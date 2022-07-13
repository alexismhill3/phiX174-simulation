# STEP 1: IMPORT PACKAGES
import pandas as pd
import numpy as np
import pinetree as pt
from sklearn.metrics import mean_squared_error
import sys

# STEP 2: RUN PINETREE FUNCTION
def run_pt(gen, pA, pB, pD, tJ, tF, tG, tH, base_dir, genomic_coords_path, date):
    
    # STEP 2.1: Create E. Coli host cell and phage genome
    CELL_VOLUME = 1.1e-15  
    model = pt.Model(cell_volume=CELL_VOLUME)
    phage = pt.Genome(name="phix_174", length=5400) # length is lengthened to increase space between gene H and terminator G

    # STEP 2.2: Read genomic coordinates from csv into a dataframe
    genomic_coords = pd.read_csv(base_dir + genomic_coords_path + "genomic_coords.csv")
    # Change locations of gene G, terminator G, and gene H to have more spacing between them
        # gH ends at 5390, starts at 4413
        # tG 4400 - 4401    
    genomic_coords.at[11, "new_start"] = 4413
    genomic_coords.at[11, "new_end"] = 5390
    
    #print(genomic_coords.at[0, "type"])
    
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
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               efficiency={"ecolipol": tH})
                               
        elif genomic_coords.at[n, "type"] == "terminator" and genomic_coords.at[n, "name"] == "J":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               efficiency={"ecolipol": tJ})
                               
        elif genomic_coords.at[n, "type"] == "terminator" and genomic_coords.at[n, "name"] == "F":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               efficiency={"ecolipol": tF})
                               
        elif genomic_coords.at[n, "type"] == "terminator" and genomic_coords.at[n, "name"] == "G":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               efficiency={"ecolipol": tG})
                               
        else:
            print(" ")

        n = n + 1

    # Register genome after promoters/terminators are added
    model.register_genome(phage)

    # Add polymerases & species
    model.add_polymerase(name="ecolipol", speed=35, footprint=35, copy_number=0)

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
    model.seed(34)

    # Run simulation
    model.simulate(time_limit=1200, time_step=5,output=base_dir + "output/" + str(date) + "/sim_" + str(sim) + "_gen_" + str(gen) + ".tsv")
    

# STEP 3: GET ERROR FUNCTION
def get_error(file, gen):
    sim = pd.read_csv(file, sep="\t")
    sim = sim.round({'time': 0})
    sim = sim[sim['time'] == 1200.0]
    sim = sim[sim.species.str.match("gene_")]

    if (gen == 0):
        if (len(sim) != 11):
            error = -1
            return (error)
        else:
            # OLD ERROR CALC - RELATIVE ABUNDANCES
            # sim["norm"] = sim['transcript'] / (sim.iloc[0]["transcript"])
            # sim["exp"] = [1, 1, 6, 6, 17, 17, 11, 5, 1, 17, 6]
            
            # NEW ERROR CALC - RAW COUNTS
            sim["norm"] = sim['transcript'] # / (sim.iloc[0]["transcript"])
            sim["exp"] = [10, 10, 60, 60, 170, 170, 110, 50, 10, 170, 60] # scaled values based on output/4152022/sim_79_gen_468.tsv
            error = round(mean_squared_error(sim.exp, sim.norm, squared=False), 5)
            return (error)

    else:
        if (len(sim) != 11):
            error = np.exp(100)
            return (error)

        else:
            # sim["norm"] = sim['transcript'] / (sim.iloc[0]["transcript"])
            # sim["exp"] = [1, 1, 6, 6, 17, 17, 11, 5, 1, 17, 6]
            sim["norm"] = sim['transcript'] # / (sim.iloc[0]["transcript"])
            sim["exp"] = [10, 10, 60, 60, 170, 170, 110, 50, 10, 170, 60] # scaled values based on output/4152022/sim_79_gen_468.tsv
            error = round(mean_squared_error(sim.exp, sim.norm, squared=False), 5)
            return (error)

# STEP 4: OPTIMIZE 
def main(u, o, date, sim, base_dir, genomic_coords_path, total_gens):
    gen = 0
    pA = np.random.normal(u, o, 1)[0]
    pB = np.random.normal(u, o, 1)[0]
    pD = np.random.normal(u, o, 1)[0]
    tJ = np.random.random_sample(1, )[0]
    tF = np.random.random_sample(1, )[0]
    tG = np.random.random_sample(1, )[0]
    tH = np.random.random_sample(1, )[0]

    params = ["pA", "pB", "pD", "tJ", "tF", "tG", 'tH']
    additive_step = [1, -1]
    report_df = pd.DataFrame(columns=['gen', 'pA', 'pB', 'pD', 'tJ', 'tF', 'tG', 'tH', 'error', 'min_error'])
    min_error = None
    min_gen = None  # generation which produced min_error

    while (gen < total_gens+1):
        
        while (gen == 0):
            try:
                run_pt(gen, pA, pB, pD, tJ, tF, tG, tH, base_dir, date)
                new_error = get_error(file=base_dir + "output/" + str(date) + "/sim_" +str(sim) +"_gen_" + str(gen) + ".tsv", gen=gen)
                
                if (new_error == -1):
                    pA = np.random.normal(u, o, 1)[0]
                    pB = np.random.normal(u, o, 1)[0]
                    pD = np.random.normal(u, o, 1)[0]
                    #step = (np.random.choice(additive_step) * 0.01)
                    step = np.random.normal(0, 0.01, 1)[0]
                    new_tJ = tJ + step
                    if(new_tJ > 1 or new_tJ < 0):
                        tJ = tJ
                    else:
                        tJ = new_tJ

                    #step = (np.random.choice(additive_step) * 0.01)
                    step = np.random.normal(0, 0.01, 1)[0]
                    new_tF = tF + step
                    if(new_tF > 1 or new_tF < 0):
                        tF = tF
                    else:
                        tF = new_tF

                    #step = (np.random.choice(additive_step) * 0.01)
                    step = np.random.normal(0, 0.01, 1)[0]
                    new_tG = tG + step
                    if(new_tG > 1 or new_tG < 0):
                        tG = tG
                    else:
                        tG = new_tG

                    #step = (np.random.choice(additive_step) * 0.01)
                    step = np.random.normal(0, 0.01, 1)[0]
                    new_tH = tH + step
                    if(new_tH > 1 or new_tH < 0):
                        tH = tH
                    else:
                        tH = new_tH

                else:
                    new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'tJ': tJ, 'tF': tF, 'tG': tG, 'tH': tH,
                               'error': new_error, 'min_error': new_error}
                    report_df = report_df.append(new_run, ignore_index=True)
                    min_gen = gen
                    min_error = report_df.at[min_gen, "min_error"]  # gen 0 has the min error so far
                    gen = 1
            except:
                pA = np.random.normal(u, o, 1)[0]
                pB = np.random.normal(u, o, 1)[0]
                pD = np.random.normal(u, o, 1)[0]
                tJ = np.random.random_sample(1, )[0]
                tF = np.random.random_sample(1, )[0]
                tG = np.random.random_sample(1, )[0]
                tH = np.random.random_sample(1, )[0]
                gen = 0

        rX = np.random.choice(params)
        
        if (rX == "pA"):
            step = np.random.normal(1, 0.1, 1)[0]
            pA = pA * (step)  
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH, base_dir = base_dir, date = date)

        elif (rX == "pB"):
            step = np.random.normal(1, 0.1, 1)[0]
            pB = pB * (step)  
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH,base_dir = base_dir, date = date)

        elif (rX == "pD"):
            step = np.random.normal(1, 0.1, 1)[0]
            pD = pD * (step)  
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH, base_dir = base_dir, date = date)

        elif (rX == "tJ"):
            step = np.random.normal(0, 0.01, 1)[0]
            new_tJ = tJ + step
            if(new_tJ > 1 or new_tJ < 0):
                tJ = tJ
            else:
                tJ = new_tJ
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH, base_dir = base_dir, date = date)

        elif (rX == "tF"):
            step = np.random.normal(0, 0.01, 1)[0]
            new_tF = tF + step
            if(new_tF > 1 or new_tF < 0):
                tF = tF
            else:
                tF = new_tF
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH, base_dir = base_dir, date = date)

        elif (rX == "tG"):
            step = np.random.normal(0, 0.01, 1)[0]
            new_tG = tG + step
            if(new_tG > 1 or new_tG < 0):
                tG = tG
            else:
                tG = new_tG
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH, base_dir = base_dir, date = date)

        elif (rX == "tH"):
            step = np.random.normal(0, 0.01, 1)[0]
            new_tH = tH + step
            if(new_tH > 1 or new_tH < 0):
                tH = tH
            else:
                tH = new_tH
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH, base_dir = base_dir, date = date)

        new_error = get_error(file=base_dir + "output/" + str(date) + "/sim_" + str(sim) + "_gen_" + str(gen) + ".tsv", gen=gen)

        if (new_error > min_error):
            # Save poor run
            new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'tJ': tJ, 'tF': tF, 'tG': tG, 'tH': tH,
                       'error': new_error, 'min_error': min_error}
            report_df = report_df.append(new_run, ignore_index=True)

            # reset all promoter & terminator values to min_gen parameters
            pA = report_df.at[min_gen, "pA"]
            pB = report_df.at[min_gen, "pB"]
            pD = report_df.at[min_gen, "pD"]
            tJ = report_df.at[min_gen, "tJ"]
            tF = report_df.at[min_gen, "tF"]
            tG = report_df.at[min_gen, "tG"]
            tH = report_df.at[min_gen, "tH"]

            # Increment generation
            gen = gen + 1

        else:
            # save error & values.
            min_error = new_error
            min_gen = gen  # update generation with the min error
            new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'tJ': tJ, 'tF': tF, 'tG': tG, 'tH': tH,
                       'error': new_error, 'min_error': min_error}
            report_df = report_df.append(new_run, ignore_index=True)
            gen = gen + 1
            

    report_df.to_csv(base_dir + "output/" + str(date) + "/sim_" +str(sim) +"_report.csv")


if __name__ == '__main__':
    # Code to read in arguments
    # uA, oA, uB, oB, uD, oD, sim, date, base_dir
    date = sys.argv[1]
    sim = sys.argv[2]
    base_dir = sys.argv[3]
    u = float(sys.argv[4])
    o = float(sys.argv[5])
    genomic_coords_path = sys.argv[6]
    total_gens = int(sys.argv[7])
    main(u, o, date, sim, base_dir, genomic_coords_path, total_gens)
    #base_dir = "/Users/tanviingle/Documents/Wilke/phix174/"
