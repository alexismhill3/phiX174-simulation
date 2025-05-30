# STEP 1: IMPORT PACKAGES
import pandas as pd
import numpy as np
import pinetree as pt
from sklearn.metrics import mean_squared_error
import sys
import multiprocessing

# Adapted from https://www.geeksforgeeks.org/parallel-processing-in-python/
class Process(multiprocessing.Process):
    def __init__(self, sim, gen, ID, pA, pB, pD, btss49, tJ, tF, tG, tH, RUT3, base_dir, date):
        super(Process, self).__init__()
        self.sim = sim
        self.gen = gen
        self.ID = ID
        self.pA = pA
        self.pB = pB
        self.pD = pD
        self.btss49 = btss49
        self.tJ = tJ
        self.tF = tF
        self.tG = tG
        self.tH = tH
        self.RUT3 = RUT3
        self.base_dir = base_dir
        self.date = date
                 
    def run(self):
        print("Running pinetree simulation with ID: {}".format(self.ID))
        run_pt(self.sim, self.gen, self.ID, self.pA, self.pB, self.pD, self.btss49, self.tJ, self.tF, self.tG, self.tH, self.RUT3, self.base_dir, self.date)

# STEP 2: RUN PINETREE FUNCTION
def run_pt(sim, gen, ID, pA, pB, pD, btss49, tJ, tF, tG, tH, RUT3, base_dir, date):
    print(gen, pA, pB, pD, btss49, tJ, tF, tG, tH, RUT3)
    print("Defining PhiX-174 genome")

    # Create host cell & genome
    CELL_VOLUME = 1.1e-15  # from T7
    model = pt.Model(cell_volume=CELL_VOLUME)
    phage = pt.Genome(name="phix_174", length=5400)

    # Read genomic coordinates from csv into dataframe -> WITH CRYPTIC PROMOTER/TERMINATOR
    genomic_coords = pd.read_csv(base_dir + "data/" + "genomic_coords_with_predicted_v2.csv")
    # Change locations of gene G, terminator G, and gene H to have more spacing between them
        # gH ends at 5390, starts at 4413
        # tG 4400 - 4401    
    genomic_coords.at[11, "new_start"] = 4413
    genomic_coords.at[11, "new_end"] = 5390

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
        
        elif genomic_coords.at[n, "type"] == "promoter" and genomic_coords.at[n, "name"] == "btss49":
            phage.add_promoter(name=genomic_coords.at[n, "type"] + "_" + genomic_coords.at[n, "name"],
                               start=genomic_coords.at[n, "new_start"],
                               stop=genomic_coords.at[n, "new_end"],
                               interactions={"ecolipol": np.multiply(btss49, 1000000)})

        else:
            print("ignoring pB2")

        n = n + 1

    # Add terminators manually (moved to intragenic regions, tH at start, and genome lengthend)
    phage.add_terminator(name="terminator_J", start=2436, stop=2437,  
                         efficiency={"ecolipol": tJ}) # After gene J
    phage.add_terminator(name="terminator_F", start=3796, stop=3797,  
                         efficiency={"ecolipol": tF})  # After gene F
    phage.add_terminator(name="terminator_G", start=4400, stop=4401,
                       efficiency={"ecolipol": tG}) # After gene G
    phage.add_terminator(name="terminator_H", start=55, stop=56,
                 efficiency={"ecolipol": tH})
    phage.add_terminator(name="terminator_RUT3", start=1229, stop=1230,
                 efficiency={"ecolipol": RUT3})

    # Register genome after promoters/terminators are added
    model.register_genome(phage)

    # Define interactions
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
    model.simulate(time_limit=500, time_step=5,output=base_dir + "output/" + str(date) + "/sim_" + str(sim) + "_gen_" + str(gen) + "_ID_" + str(ID) + ".tsv")
    print("Simulation successful!")
    
def run_parallel(sim, gen, pA, pB, pD, btss49, tJ, tF, tG, tH, RUT3, base_dir, date, trials):
    processes = []
    for i in range(trials):
      p = Process(sim, gen, i, pA, pB, pD, btss49, tJ, tF, tG, tH, RUT3, base_dir, date)
      p.start()
      processes.append(p)
      
    for p in processes:
      p.join()
    
    print("Finished running all sims")

# STEP 3: GET ERROR FUNCTION
def get_error(sim, gen, trials, base_dir, date):
    results = pd.DataFrame()
    for i in range(trials): # Average over simulation trials
        file = base_dir + "output/" + str(date) + "/sim_" + str(sim) + "_gen_" + str(gen) + "_ID_" + str(i) + ".tsv"
        tmp = pd.read_csv(file, sep="\t")
        tmp = tmp.round({'time': 0})
        tmp = tmp[tmp['time'] == 500.0]
        tmp = tmp[tmp.species.str.match("gene_")]
        results = results.append(tmp, ignore_index=True)
    results = results.groupby("species").mean()
    
    if (gen == 0):
        if (len(results) != 11):
            error = -1
            return (error)
        else:
            #AA*BCKDEJFGH
            #results["exp"] = [17, 17, 102, 102, 289, 289, 187, 85, 17, 289, 102]
            results["exp"] = [27, 27, 120, 81, 228, 201, 231, 118, 76, 264, 123] # use RNAseq data values
            error = round(mean_squared_error(results.exp, results.transcript, squared=False), 5)
            return (error)

    else:
        if (len(results) != 11):
            error = np.exp(100)
            print("len != 11, report bogus error")
            return (error)

        else:
            #results["exp"] = [17, 17, 102, 102, 289, 289, 187, 85, 17, 289, 102]
            results["exp"] = [27, 27, 120, 81, 228, 201, 231, 118, 76, 264, 123]
            error = round(mean_squared_error(results.exp, results.transcript, squared=False), 5)
            return (error)

# STEP 4:
def main(u, o, o_prom, o_term, date, sim, base_dir, trials, iterations):
    gen = 0
    pA = np.random.normal(u, o, 1)[0]
    pB = np.random.normal(u, o, 1)[0]
    pD = np.random.normal(u, o, 1)[0]
    btss49 = np.random.normal(u, o, 1)[0]
    tJ = np.random.random_sample(1, )[0]
    tF = np.random.random_sample(1, )[0]
    tG = np.random.random_sample(1, )[0]
    tH = np.random.random_sample(1, )[0]
    RUT3 = np.random.random_sample(1, )[0]

    params = ["pA", "pB", "pD", "btss49", "tJ", "tF", "tG", 'tH', "RUT3"]
    additive_step = [1, -1]
    report_df = pd.DataFrame(columns=['gen', 'pA', 'pB', 'pD', "btss49", 'tJ', 'tF', 'tG', 'tH', "RUT3", 'error', 'min_error'])
    min_error = None
    min_gen = None  # generation which produced min_error

    while (gen < iterations):
        print("gen = ", gen)

        while (gen == 0):
            try:
                print(f"\n")
                print("trying")
                print(f"\n")
                
                ##----- PARALLELIZE HERE -----##
                run_parallel(sim, gen, pA, pB, pD, btss49, tJ, tF, tG, tH, RUT3, base_dir, date, trials)
                new_error = get_error(sim, gen, trials, base_dir, date)
                print("new_error = ", new_error)
                ##----------------------------##
                
                if (new_error == -1):
                    print("poor starting values chosen, retrying")
                    print(f"\n")
                    pA = np.random.normal(u, o, 1)[0]
                    pB = np.random.normal(u, o, 1)[0]
                    pD = np.random.normal(u, o, 1)[0]
                    btss49 = np.random.normal(u, o, 1)[0]
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
                        
                    step = np.random.normal(0, 0.01, 1)[0]
                    new_RUT3 = RUT3 + step
                    if(new_RUT3 > 1 or new_RUT3 < 0):
                        RUT3 = RUT3
                    else:
                        RUT3 = new_RUT3

                else:
                    print(f"\n")
                    print("good starting values")
                    print(f"\n")
                    print(gen, pA, pB, pD, btss49, tJ, tF, tG, tH, RUT3)
                    new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'btss49': btss49, 'tJ': tJ, 'tF': tF, 'tG': tG, 'tH': tH, 'RUT3': RUT3,
                               'error': new_error, 'min_error': new_error}
                    report_df = report_df.append(new_run, ignore_index=True)
                    min_gen = gen
                    min_error = report_df.at[min_gen, "min_error"]  # gen 0 has the min error so far
                    gen = 1
            except Exception as e:
                print(f"\n")
                print(e)
                print("except, retry")
                print(f"\n")
                pA = np.random.normal(u, o, 1)[0]
                pB = np.random.normal(u, o, 1)[0]
                pD = np.random.normal(u, o, 1)[0]
                btss49 = np.random.normal(u, o, 1)[0]
                tJ = np.random.random_sample(1, )[0]
                tF = np.random.random_sample(1, )[0]
                tG = np.random.random_sample(1, )[0]
                tH = np.random.random_sample(1, )[0]
                RUT3 = np.random.random_sample(1, )[0]
                gen = 0

        rX = np.random.choice(params)
        print(f"\n")
        print("regulatory element to be optimized: ", rX)
        print(f"\n")
  
        print(pA, pB, pD, btss49, tJ, tF, tG, tH, RUT3)
        if (rX == "pA"):
            step = np.random.normal(0, o_prom, 1)[0]
            pA = pA * (2 ** step)  
            print("new pA = ", pA)

        elif (rX == "pB"):
            step = np.random.normal(0, o_prom, 1)[0]
            pB = pB * (2 ** step)  
            print("new pB = ", pB)

        elif (rX == "pD"):
            step = np.random.normal(0, o_prom, 1)[0]
            pD = pD * (2 ** step)  
            print("new pD = ", pD)
        
        elif (rX == "btss49"):
            step = np.random.normal(0, o_prom, 1)[0]
            btss49 = btss49 * (2 ** step)  
            print("new btss49 = ", btss49)

        elif (rX == "tJ"):
            step = np.random.normal(0, o_term, 1)[0]
            new_tJ = tJ + step
            if(new_tJ > 1 or new_tJ < 0):
                tJ = tJ
            else:
                tJ = new_tJ
            print("new tJ = ", tJ)

        elif (rX == "tF"):
            step = np.random.normal(0, o_term, 1)[0]
            new_tF = tF + step
            if(new_tF > 1 or new_tF < 0):
                tF = tF
            else:
                tF = new_tF
            print("new tF = ", tF)

        elif (rX == "tG"):
            step = np.random.normal(0, o_term, 1)[0]
            new_tG = tG + step
            if(new_tG > 1 or new_tG < 0):
                tG = tG
            else:
                tG = new_tG
            print("new tG = ", tG)

        elif (rX == "tH"):
            step = np.random.normal(0, o_term, 1)[0]
            new_tH = tH + step
            if(new_tH > 1 or new_tH < 0):
                tH = tH
            else:
                tH = new_tH
            print("new tH = ", tH)
        
        elif (rX == "RUT3"):
            step = np.random.normal(0, o_term, 1)[0]
            new_RUT3 = RUT3 + step
            if(new_RUT3 > 1 or new_RUT3 < 0):
                RUT3 = RUT3
            else:
                RUT3 = new_RUT3
            print("new RUT3 = ", RUT3)
            
        ##----- PARALLELIZE HERE -----##
        run_parallel(sim, gen, pA, pB, pD, btss49, tJ, tF, tG, tH, RUT3, base_dir, date, trials)
        new_error = get_error(sim, gen, trials, base_dir, date)
        print("new_error = ", new_error)
        ##----------------------------##

        if (new_error > min_error):
            print(f"new error: {new_error}, min error: {min_error} - resetting {rX}")
            # Save poor run
            new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'btss49': btss49, 'tJ': tJ, 'tF': tF, 'tG': tG, 'tH': tH, 'RUT3': RUT3,
                       'error': new_error, 'min_error': min_error}
            report_df = report_df.append(new_run, ignore_index=True)

            # reset all promoter & terminator values to min_gen parameters
            pA = report_df.at[min_gen, "pA"]
            pB = report_df.at[min_gen, "pB"]
            pD = report_df.at[min_gen, "pD"]
            btss49 = report_df.at[min_gen, "btss49"]
            tJ = report_df.at[min_gen, "tJ"]
            tF = report_df.at[min_gen, "tF"]
            tG = report_df.at[min_gen, "tG"]
            tH = report_df.at[min_gen, "tH"]
            RUT3 = report_df.at[min_gen, "RUT3"]

            # Increment generation
            gen = gen + 1
            print(f"\n")

        else:
            print(f"new error: {new_error}, min error: {min_error} - keeping {rX}")
            # save error & values.
            min_error = new_error
            min_gen = gen  # update generation with the min error
            #new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'tJ': tJ, 'tF': tF, 'tG': tG, 'tH': tH,
            #           'error': new_error, 'min_error': min_error}
            new_run = {'gen': gen, 'pA': pA, 'pB': pB, 'pD': pD, 'btss49': btss49, 'tJ': tJ, 'tF': tF, 'tG': tG, 'tH': tH, 'RUT3': RUT3,
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
    o_prom = float(sys.argv[6])
    o_term = float(sys.argv[7])
    trials = 5
    iterations = 4000
        
    main(u, o, o_prom, o_term, date, sim, base_dir, trials, iterations)
