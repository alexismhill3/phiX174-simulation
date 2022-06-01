# STEP 1: IMPORT PACKAGES
import pandas as pd
import numpy as np
import pinetree as pt
from sklearn.metrics import mean_squared_error
import sys

# STEP 2: RUN PINETREE FUNCTION
def run_pt(gen, pA, pB, pD, base_dir, date):
    print(gen, pA, pB, pD)
    print("Defining PhiX-174 genome")

    # Create host cell & genome
    CELL_VOLUME = 1.1e-15  # from T7
    model = pt.Model(cell_volume=CELL_VOLUME)
    phage = pt.Genome(name="phix_174", length=5386)

    # Read genomic coordinates from csv into dataframe
    genomic_coords = pd.read_csv(base_dir + "output/" + "genomic_coords.csv")
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
    phage.add_terminator(name="terminator_J", start=2402, stop=2403,  # Right before gene F start=2404, stop=3687,
                         efficiency={"ecolipol": 0.3})  # 0.7
    phage.add_terminator(name="terminator_F", start=3796, stop=3797,  # Right before gene G start=3798, stop=4325
                         efficiency={"ecolipol": 0.2})  # 0.8
    phage.add_terminator(name="terminator_G", start=4332, stop=4333,  # Right before gene H start=4334, stop=5320
                         efficiency={"ecolipol": 0.3})  # 0.6
    phage.add_terminator(name="terminator_H", start=5321, stop=5322,  # Right after gene H
                         efficiency={"ecolipol": 0.7})  # 0.3

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
    model.simulate(time_limit=1200, time_step=5,
                   output= base_dir + "output/" + str(date) + "/sim_" +str(sim) +"_gen_" + str(gen) + ".tsv")  # TODO change FILE PATHS
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
    pA = np.exp(np.random.normal(u, o, 1)[0])
    pB = np.exp(np.random.normal(u, o, 1)[0])
    pD = np.exp(np.random.normal(u, o, 1)[0])
    params = ["pA", "pB", "pD"]
    report_df = pd.DataFrame(columns=['gen', 'pA', 'pB', 'pD', 'error', 'min_error'])

    min_error = None
    min_gen = None  # generation which produced min_error
    new_error = None

    while (gen < 501):
        print("gen = ", gen)

        while (gen == 0):
            try:
                run_pt(gen, pA, pB, pD, base_dir, date)
                new_error = get_error(file=base_dir + "output/" + str(date) + "/sim_" +str(sim) +"_gen_" + str(gen) + ".tsv", gen=gen)
                print("new_error = ", new_error)

                if (new_error == -1):
                    print("poor starting values chosen, retrying")
                    print(f"\n")
                    pA = np.random.normal(12, 3, 1)[0]
                    pB = np.random.normal(12, 3, 1)[0]
                    pD = np.random.normal(12, 3, 1)[0]

                else:
                    print(f"\n")
                    print("good starting values")
                    print(f"\n")
                    print(gen, pA, pB, pD)
                    new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'error': new_error, 'min_error': new_error}
                    report_df = report_df.append(new_run, ignore_index=True)
                    min_gen = gen
                    min_error = report_df.at[min_gen, "min_error"]  # gen 0 has the min error so far
                    gen = 1
            except:
                pA = np.random.normal(12, 3, 1)[0]
                pB = np.random.normal(12, 3, 1)[0]
                pD = np.random.normal(12, 3, 1)[0]
                gen = 0

        pX = np.random.choice(params)
        print(f"\n")
        print("Promoter to be optimized: ", pX)
        print(f"\n")

        if (pX == "pA"):
            step = np.random.normal(0, 0.1, 1)[0]
            pA = pA * (step + 1)  #
            print("new pA = ", pA)
            print(pA, pB, pD)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, base_dir = base_dir, date = date)

        elif (pX == "pB"):
            step = np.random.normal(0, 0.1, 1)[0]
            pB = pB * (step + 1)  #
            print("new pB = ", pB)
            print(pA, pB, pD)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, base_dir = base_dir, date = date)

        elif (pX == "pD"):
            step = np.random.normal(0, 0.1, 1)[0]
            pD = pD * (step + 1)  #
            print("new pD = ", pB)
            print(pA, pB, pD)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD, base_dir = base_dir, date = date)

        new_error = get_error(file=base_dir + "output/" + str(date) + "/sim_" +str(sim) +"_gen_" + str(gen) + ".tsv", gen=gen)

        if (new_error > min_error):
            # Save poor run
            new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'error': new_error, 'min_error': min_error}
            report_df = report_df.append(new_run, ignore_index=True)

            # reset all promoter values to min_gen parameters
            pA = report_df.at[min_gen, "pA"]
            pB = report_df.at[min_gen, "pB"]
            pD = report_df.at[min_gen, "pD"]

            # Increment generation
            gen = gen + 1

        else:
            # save error & values.
            min_error = new_error
            min_gen = gen  # update generation with the min error
            new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'error': new_error, 'min_error': min_error}
            report_df = report_df.append(new_run, ignore_index=True)
            gen = gen + 1
            print(f"\n")

    report_df.to_csv(base_dir + "output/" + str(date) + "/sim_" +str(sim) +"_report" + ".csv")



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
