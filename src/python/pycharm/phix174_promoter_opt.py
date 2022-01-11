# STEP 1: Import packages & Define base directory
import pandas as pd
import numpy as np
import pinetree as pt
from sklearn.metrics import mean_squared_error

# STEP 2: run_pt() function runs the pinetree simulation
def run_pt(gen, pA, pB, pD):
    print(gen, pA, pB, pD)
    print("Defining PhiX-174 genome")

    # Create host cell & genome
    CELL_VOLUME = 1.1e-15  # from T7
    model = pt.Model(cell_volume=CELL_VOLUME)
    phage = pt.Genome(name="phix_174", length=5386)

    # Read genomic coordinates from csv into dataframe
    genomic_coords = pd.read_csv(base_dir + "output/" + "genomic_coords.csv")
    print(genomic_coords.at[0, "type"])

    # Add genomic ns (loop through ^ df)
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
                               interactions={"ecolipol": pA})

        elif genomic_coords.at[n, "type"] == "promoter" and genomic_coords.at[n, "name"] == "B1":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               interactions={"ecolipol": pB})

        elif genomic_coords.at[n, "type"] == "promoter" and genomic_coords.at[n, "name"] == "D":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               interactions={"ecolipol": pD})

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
                   output=base_dir + "output/" + str(date) + "/sim_" + str(sim) + "_gen_" + str(gen) + "_ptrun.tsv")
    print("Simulation successful!")


# STEP 3: get_error() calculates the RMSE between experimental & simulated transcript abundances
# Experimental abundances are hardcoded
def get_error(file):
    sim = pd.read_csv(file, sep="\t")
    sim = sim.round({'time': 0})
    sim = sim[sim['time'] == 1200.0]
    sim = sim[sim.species.str.match("gene_")]
    if (len(sim) != 11):
        error = -1
        return (error)
    else:
        sim["norm"] = sim['transcript'] / (sim.iloc[0]["transcript"])
        sim["exp"] = [1, 1, 6, 6, 17, 17, 11, 5, 1, 17, 6]
        error = mean_squared_error(sim.exp, sim.norm, squared=False)
        return (error)


if __name__ == '__main__':
    # Code to read in arguments
    # uA, oA, uB, oB, uD, oD, sim, date, base_dir

    base_dir = "/Users/tanviingle/Documents/Wilke/phix174/"

    gen = 0
    pA = np.exp(np.random.normal(uA, oA, 1)[0])
    pB = np.exp(np.random.normal(uB, oB, 1)[0])
    pD = np.exp(np.random.normal(uD, oD, 1)[0])
    j = 0
    error = 1e10
    step = 1e2
    params = ["pA", "pB", "pD"]

    report_df = pd.DataFrame(columns=['gen', 'pA', 'pB', 'pD', 'error'])
    first_run = {'gen': "NA", 'pA': pA, 'pB': pB, 'pD': pD, 'error': error}
    report_df = report_df.append(first_run, ignore_index=True)

    old_error = error
    new_error = error

    no_change = 0

    while (gen < 100):
        # STEP 1: Randomly select promoter to optimize
        print(gen)

        while (gen == 0 and j == 0):
            try:
                run_pt(gen, pA, pB, pD)
                j = 1
                gen = 1
                print("good starting values")
                print(gen, pA, pB, pD)
                print(f"\n")
            except:
                print("poor starting values chosen, retrying")
                pA = np.exp(np.random.normal(uA, oA, 1)[0])
                pB = np.exp(np.random.normal(uB, oB, 1)[0])
                pD = np.exp(np.random.normal(uD, oD, 1)[0])
                j = 0

        if (no_change == 0):
            print('pick random promoter')
            pX = np.random.choice(params)
            print(pX)

        # STEP 2: Run pinetree with selected promoter & value
        # Increment pX by oX
        if (pX == "pA"):
            print(pA)
            print("increment pA")
            pA = pA + np.exp(oA)
            print("new pA = ", pA)
            print(pA, pB, pD)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD)

        if (pX == "pB"):
            print(pB)
            print("increment pB")
            pB = pB + np.exp(oB)
            print("new pB = ", pB)
            print(pA, pB, pD)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD)

        if (pX == "pD"):
            print(pD)
            print("increment pD")
            pD = pD + np.exp(oD)
            print("new pD = ", pD)
            print(pA, pB, pD)
            run_pt(gen=gen, pA=pA, pB=pB, pD=pD)

        # STEP 3: Calculate Error
        old_error = new_error
        new_error = get_error(file=base_dir + "output/" + str(date) + "/sim_" + str(sim) + "_gen_" + str(gen) + "_ptrun.tsv")

        # STEP 4: Compare Old Error to New Error;
        if ((new_error == -1 or new_error >= old_error)):
            if (no_change >= 4):
                no_change = 0
                continue
            print("entering error loop")
            no_change = no_change + 1
            print("new_error = " + str(new_error) + ", old_error = " + str(old_error) + "\n")
            new_error = old_error
            old_error = report_df.at[gen - 1, "error"]

            if (pX == "pA"):
                step = oA * (np.random.normal(0.5, 0.2, 1)[0])
                pA = report_df.at[gen - 1, "pA"] + np.exp(step)
                print(step)
                # print(pA, pB, pD)
                new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'error': new_error}
                report_df = report_df.append(new_run, ignore_index=True)
                gen = gen + 1
                print(report_df)
                continue

            if (pX == "pB"):
                step = oB * (np.random.normal(0.5, 0.2, 1)[0])
                pB = report_df.at[gen - 1, "pB"] + np.exp(step)
                print(step)
                # print(pA, pB, pD)
                new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'error': new_error}
                report_df = report_df.append(new_run, ignore_index=True)
                gen = gen + 1
                print(report_df)
                continue

            if (pX == "pD"):
                step = oD * (np.random.normal(0.5, 0.2, 1)[0])
                pD = report_df.at[gen - 1, "pD"] + np.exp(step)
                print(step)
                # print(pA, pB, pD)
                new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'error': new_error}
                report_df = report_df.append(new_run, ignore_index=True)
                gen = gen + 1
                print(report_df)
                continue

        # print (step)
        new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'error': new_error}
        report_df = report_df.append(new_run, ignore_index=True)

        print(pA, pB, pD)
        gen = gen + 1
        no_change = 0
        print(f"\n")

    report_df.to_csv(base_dir + "output/" + str(date) + "/sim_" + str(sim) + "_report.csv")