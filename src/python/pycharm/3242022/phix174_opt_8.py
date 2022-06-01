# STEP 1: IMPORT PACKAGES
import pandas as pd
import numpy as np
import pinetree as pt
from sklearn.metrics import mean_squared_error
import sys

# STEP 2: RUN PINETREE FUNCTION
def run_pt(gen, pA, pB, pD, tJ, tF, tG, tH, base_dir, date):
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
    model.add_polymerase(name="ecolipol", speed=35, footprint=35, copy_number=0)
    model.add_species("bound_ecolipol", 2000)  # initialization, 1800
    model.add_species("ecoli_genome", 0)
    model.add_species("ecoli_transcript", 0)
    model.add_reaction(1e7, ["ecolipol", "ecoli_genome"], ["bound_ecolipol"])  # 1e6
    model.add_reaction(0.04, ["bound_ecolipol"], ["ecolipol", "ecoli_genome", "ecoli_transcript"])
    model.add_ribosome(10, 30, 100)
    model.add_species("bound_ribosome", 100)
    model.seed(34)

    # Run simulation
    print("running simulation")
    model.simulate(time_limit=1200, time_step=5,output=base_dir + "output/" + str(date) + "/sim_" + str(sim) + "_gen_" + str(gen) + ".tsv")
    print("Simulation successful!")

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
            sim["norm"] = sim['transcript'] / (sim.iloc[0]["transcript"])
            sim["exp"] = [1, 1, 6, 6, 17, 17, 11, 5, 1, 17, 6]
            error = round(mean_squared_error(sim.exp, sim.norm, squared=False), 5)
            return (error)

    else:
        if (len(sim) != 11):
            error = np.exp(100)
            print("len != 11, report bogus error")
            return (error)

        else:
            sim["norm"] = sim['transcript'] / (sim.iloc[0]["transcript"])
            sim["exp"] = [1, 1, 6, 6, 17, 17, 11, 5, 1, 17, 6]
            error = round(mean_squared_error(sim.exp, sim.norm, squared=False), 5)
            return (error)

# STEP 4:
def main(u, o, date, sim, base_dir):
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

    while (gen < 501):
        print("gen = ", gen)

        while (gen == 0):
            try:
                print(f"\n")
                print("trying")
                print(f"\n")
                run_pt(gen, pA, pB, pD, tJ, tF, tG, tH, base_dir, date)
                new_error = get_error(file=base_dir + "output/" + str(date) + "/sim_" +str(sim) +"_gen_" + str(gen) + ".tsv", gen=gen)
                print("new_error = ", new_error)

                if (new_error == -1):
                    print("poor starting values chosen, retrying")
                    print(f"\n")
                    pA = np.random.normal(u, o, 1)[0]
                    pB = np.random.normal(u, o, 1)[0]
                    pD = np.random.normal(u, o, 1)[0]
                    step = (np.random.choice(additive_step) * 0.01)
                    tJ = tJ + step
                    if (tJ > 1):
                        tJ = tJ - (step)
                    if (tJ < 0):
                        tJ = tJ + (step)

                    step = (np.random.choice(additive_step) * 0.01)
                    tF = tF + step
                    if (tF > 1):
                        tF = tF - (step)
                    if (tF < 0):
                        tF = tF + (step)

                    step = (np.random.choice(additive_step) * 0.01)
                    tG = tG + step
                    if (tG > 1):
                        tG = tG - (step)
                    if (tG < 0):
                        tG = tG + (step)

                    step = (np.random.choice(additive_step) * 0.01)
                    tH = tH + step
                    if (tH > 1):
                        tH = tH - (step)
                    if (tH < 0):
                        tH = tH + (step)

                else:
                    print(f"\n")
                    print("good starting values")
                    print(f"\n")
                    print(gen, pA, pB, pD, tJ, tF, tG, tH)
                    new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'tJ': tJ, 'tF': tF, 'tG': tG, 'tH': tH,
                               'error': new_error, 'min_error': new_error}
                    report_df = report_df.append(new_run, ignore_index=True)
                    min_gen = gen
                    min_error = report_df.at[min_gen, "min_error"]  # gen 0 has the min error so far
                    gen = 1
            except:
                print(f"\n")
                print("except, retry")
                print(f"\n")
                pA = np.random.normal(u, o, 1)[0]
                pB = np.random.normal(u, o, 1)[0]
                pD = np.random.normal(u, o, 1)[0]
                tJ = np.random.random_sample(1, )[0]
                tF = np.random.random_sample(1, )[0]
                tG = np.random.random_sample(1, )[0]
                tH = np.random.random_sample(1, )[0]
                gen = 0

        rX = np.random.choice(params)
        print(f"\n")
        print("regulatory element to be optimized: ", rX)
        print(f"\n")

        if (rX == "pA"):
            step = np.random.normal(0, 0.1, 1)[0]
            pA = pA * (step + 1)  #
            print("new pA = ", pA)
            print(pA, pB, pD, tJ, tF, tG, tH)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH, base_dir = base_dir, date = date)

        elif (rX == "pB"):
            step = np.random.normal(0, 0.1, 1)[0]
            pB = pB * (step + 1)  #
            print("new pB = ", pB)
            print(pA, pB, pD, tJ, tF, tG, tH)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH,base_dir = base_dir, date = date)

        elif (rX == "pD"):
            step = np.random.normal(0, 0.1, 1)[0]
            pD = pD * (step + 1)  #
            print("new pD = ", pD)
            print(pA, pB, pD, tJ, tF, tG, tH)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH, base_dir = base_dir, date = date)

        elif (rX == "tJ"):
            step = 0.01
            direction = np.random.choice(additive_step)
            tJ = tJ + (direction * step)
            if (tJ > 1):
                tJ = tJ - (step)
            if (tJ < 0):
                tJ = tJ + (step)
            print("new tJ = ", tJ)
            print(pA, pB, pD, tJ, tF, tG, tH)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH, base_dir = base_dir, date = date)

        elif (rX == "tF"):
            step = 0.01
            direction = np.random.choice(additive_step)
            tF = tF + (direction * step)
            if (tF > 1):
                tF = tF - (step)
            if (tF < 0):
                tF = tF + (step)
            print("new tF = ", tF)
            print(pA, pB, pD, tJ, tF, tG, tH)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH, base_dir = base_dir, date = date)

        elif (rX == "tG"):
            step = 0.01
            direction = np.random.choice(additive_step)
            tG = tG + (direction * step)
            if (tG > 1):
                tG = tG - (step)
            if (tG < 0):
                tG = tG + (step)
            print("new tG = ", tG)
            print(pA, pB, pD, tJ, tF, tG, tH)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, tJ=tJ, tF=tF, tG=tG, tH=tH, base_dir = base_dir, date = date)

        elif (rX == "tH"):
            step = 0.01
            direction = np.random.choice(additive_step)
            tH = tH + (direction * step)
            if (tH > 1):
                tH = tH - (step)
            if (tH < 0):
                tH = tH + (step)
            print("new tH = ", tH)
            print(pA, pB, pD, tJ, tF, tG, tH)
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
            print(f"\n")

    report_df.to_csv(base_dir + "output/" + str(date) + "/sim_" +str(sim) +"_report.csv")


if __name__ == '__main__':
    # Code to read in arguments
    # uA, oA, uB, oB, uD, oD, sim, date, base_dir
    date = sys.argv[1]
    sim = sys.argv[2]
    base_dir = sys.argv[3]
    u = float(sys.argv[4])
    o = float(sys.argv[5])
    main(u, o, date, sim, base_dir)
    #base_dir = "/Users/tanviingle/Documents/Wilke/phix174/"
