import pandas as pd
import numpy as np
import pinetree as pt
import sys
import multiprocessing

def run_pt(sim, ID, pA, pB, pD, tJ, tF, tG, tH, base_dir, output_dir):
    print(ID, pA, pB, pD, tJ, tF, tG, tH)
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

    # Add terminators manually (moved to intragenic regions, tH at start, and genome lengthend)
    phage.add_terminator(name="terminator_J", start=2436, stop=2437,  
                         efficiency={"ecolipol": tJ}) # After gene J
    phage.add_terminator(name="terminator_F", start=3796, stop=3797,  
                         efficiency={"ecolipol": tF})  # After gene F
    phage.add_terminator(name="terminator_G", start=4400, stop=4401,
                       efficiency={"ecolipol": tG}) # After gene G
    phage.add_terminator(name="terminator_H", start=55, stop=56,
                 efficiency={"ecolipol": tH})

    # Register genome after promoters/terminators are added
    model.register_genome(phage)

    # Define interactions
    # Add polymerases & species ## No readthrough!
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
    #model.seed(34)
    print("Not using seeds")

    # Run simulation
    print("running simulation")
    model.simulate(time_limit=500, time_step=5,output=output_dir + "/linear_" + sim + "_ID_" + str(ID) + ".tsv")
    print("Simulation successful!")
    
    
# Adapted from https://www.geeksforgeeks.org/parallel-processing-in-python/
class Process(multiprocessing.Process):
    def __init__(self, sim, ID, pA, pB, pD, tJ, tF, tG, tH, base_dir, output_dir):
        super(Process, self).__init__()
        self.sim = sim
        self.ID = ID
        self.pA = pA
        self.pB = pB
        self.pD = pD
        self.tJ = tJ
        self.tF = tF
        self.tG = tG
        self.tH = tH
        self.base_dir = base_dir
        self.output_dir = output_dir
                 
    def run(self):
        print("Running pinetree simulation with ID: {}".format(self.ID))
        run_pt(self.sim, self.ID, self.pA, self.pB, self.pD, self.tJ, self.tF, self.tG, self.tH, self.base_dir, self.output_dir)
  

if __name__ == "__main__":
    sim = sys.argv[1]
    base_dir = sys.argv[2]
    output_dir = sys.argv[3]
    
    # Best set of parameters from fitting to qPCR data (11162022)
    #pA = 8.448189	
    pA = 0.1
    pB = 3.230365	
    pD = 34.62136	
    tJ = 0.3467087	
    tF = 0.5057597	
    #tG = 0.8051313
    tG = float(sys.argv[4])
    #tH = 0.9287955
    tH = float(sys.argv[5])
    
    processes = []
    for i in range(5):
      p = Process(sim, i, pA, pB, pD, tJ, tF, tG, tH, base_dir, output_dir)
      p.start()
      processes.append(p)
      
    for p in processes:
      p.join()
    
    print("Finished running all sims")
