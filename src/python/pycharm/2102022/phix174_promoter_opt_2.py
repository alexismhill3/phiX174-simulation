import pandas as pd
import numpy as np
import pinetree as pt
from sklearn.metrics import mean_squared_error
import sys

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

    # Add genomic ns (loop through ^ df); hardcode necessary strengths according to preomtimized_model
    ## for length of genomic_coords, add elements
    # for n in genomic_coords:{
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
                   output=base_dir + "output/" + str(date) + "/scenario_" + str(scenario) + "/sim_" + str(sim) + "_gen_" + str(gen) + "_ptrun.tsv")  # TODO change limit
    print("Simulation successful!")


# STEP 3: get_error() calculates the RMSE between experimental & simulated transcript abundances
# Experimental Q-PCR abundances are hardcoded
def get_error(file, input_report_df, input_gen):
    sim = pd.read_csv(file, sep="\t")
    sim = sim.round({'time': 0})
    sim = sim[sim['time'] == 1200.0]
    sim = sim[sim.species.str.match("gene_")]

    if (input_gen == 0):
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
            error = input_report_df.at[input_gen - 1, "error"]
            print("len != 11, reporting previous error")
            # error = report_df.at[gen-1, "error"]
            # while (error == -1):
            # error = report_df.at[gen-1, "error"]
            # serror = -1
            return (error)

        else:
            sim["norm"] = sim['transcript'] / (sim.iloc[0]["transcript"])
            sim["exp"] = [1, 1, 6, 6, 17, 17, 11, 5, 1, 17, 6]
            error = round(mean_squared_error(sim.exp, sim.norm, squared=False), 5)
            return (error)


def main(uA, oA, uB, oB, uD, oD, scenario, sim, date, base_dir):
    gen = 0
    pA = np.exp(np.random.normal(uA, oA, 1)[0])
    pB = np.exp(np.random.normal(uB, oB, 1)[0])
    pD = np.exp(np.random.normal(uD, oD, 1)[0])
    j = 0

    params = ["pA", "pB", "pD"]
    report_df = pd.DataFrame(columns=['gen', 'pA', 'pB', 'pD', 'error'])

    old_error = None
    new_error = None

    no_change = 0

    # while (old_error >= new_error and no_change <= 10):
    while (gen < 301):
        # STEP 1: Randomly select promoter to optimize; if gen == 10, quit.
        print(gen)
        while (gen == 0 and j == 0):
            try:
                run_pt(gen, pA, pB, pD)
                new_error = get_error(file=base_dir + "output/" + str(date) + "/scenario_" + str(scenario) + "/sim_" + str(sim) + "_gen_" + str(gen) + "_ptrun.tsv",
                              input_report_df=report_df,
                              input_gen=gen)
                if (new_error == -1):
                    print("poor starting values chosen, retrying")
                    pA = np.exp(np.random.normal(uA, oA, 1)[0])
                    pB = np.exp(np.random.normal(uB, oB, 1)[0])
                    pD = np.exp(np.random.normal(uD, oD, 1)[0])
                    j = 0

                else:
                    print("good starting values")
                    print(gen, pA, pB, pD)
                    new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'error': new_error}
                    report_df = report_df.append(new_run, ignore_index=True)
                    j = 1
                    gen = 1
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
        new_error = get_error(file=base_dir + "output/" + str(date) + "/scenario_" + str(scenario) + "/sim_" + str(sim) + "_gen_" + str(gen) + "_ptrun.tsv",
                              input_report_df=report_df,
                              input_gen=gen)

        # STEP 4: Compare Old Error to New Error;
        if ((new_error >= old_error)):  # new_error == -1 or
            if (no_change >= 5):
                no_change = 0
                continue
            print("entering error loop")
            no_change = no_change + 1
            print("new_error = " + str(new_error) + ", old_error = " + str(old_error) + "\n")
            new_error = old_error
            old_error = report_df.at[gen - 1, "error"]

            if (pX == "pA"):
                step = (np.random.normal(0, 0.2, 1)[0])
                while (abs(step) > 1 or step == 0):
                    step = np.random.normal(0, 0.2, 1)[0]
                print("step = ", step)
                print("increment by = ", step * np.exp(oA))
                pA = report_df.at[gen - 1, "pA"] + (step * np.exp(oA))
                new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'error': new_error}
                report_df = report_df.append(new_run, ignore_index=True)
                gen = gen + 1
                print(report_df)
                continue

            if (pX == "pB"):
                step = (np.random.normal(0, 0.2, 1)[0])
                while (abs(step) > 1 or step == 0):
                    step = np.random.normal(0, 0.2, 1)[0]
                print("step = ", step)
                print("increment by = ", step * np.exp(oB))
                pB = report_df.at[gen - 1, "pB"] + (step * np.exp(oB))
                new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'error': new_error}
                report_df = report_df.append(new_run, ignore_index=True)
                gen = gen + 1
                print(report_df)
                continue

            if (pX == "pD"):
                step = (np.random.normal(0, 0.2, 1)[0])
                while (abs(step) > 1 or step == 0):
                    step = np.random.normal(0, 0.2, 1)[0]
                print("step = ", step)
                print("increment by = ", step * np.exp(oD))
                pD = report_df.at[gen - 1, "pD"] + (step * np.exp(oD))
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

    report_df.to_csv(base_dir + "output/" + str(date) + "/scenario_" + str(scenario) + "/sim_" + str(sim) + "_report.csv")



if __name__ == '__main__':
    # Code to read in arguments
    # uA, oA, uB, oB, uD, oD, sim, date, base_dir
    scenario = sys.argv[1]
    sim = sys.argv[2]
    date = sys.argv[3]
    base_dir = sys.argv[4]
    uA = float(sys.argv[5])
    oA = float(sys.argv[6])
    uB = float(sys.argv[7])
    oB = float(sys.argv[8])
    uD = float(sys.argv[9])
    oD = float(sys.argv[10])
    main(uA, oA, uB, oB, uD, oD, scenario, sim, date, base_dir)
    #base_dir = "/Users/tanviingle/Documents/Wilke/phix174/"









