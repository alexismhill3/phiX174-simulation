# OBJECTIVE: This python script creates a list of jobs to be run in parallel. 
# ARGUMENTS: Described below. 
# OUTPUT: jobs_wt.txt in the working directory. The 'wt' means wildtype PhiX174. 
# WARNING: Make sure you delete the existing jobs_wt.txt before reruning this script. 
## Otherwise it will append everything to the end of the existing jobs_wt.txt file. 
#----------------------------------------------------------------------------------------------------------------#

# STEP 1: IMPORT LIBRARIES
import numpy as np
import pandas as pd

# STEP 2: SET PARAMETERS; Arguments match the calls in the script's main() function
python_version = "python3" 
script = "phix174_opt_12_wt.py" # the name of the optimization script
date = "7132022" # the date the optimization is run; the working directory should have a folder output/date/. 
sim = 1 # the first simulation number
base_dir = "/stor/work/Wilke/tingle/phix174/src/python/pycharm/new_scripts/" # the working directory. Contains folder output
u = "12" # mean 
o = "3" # st. deviation 
genomic_coords_path = "/stor/work/Wilke/tingle/phix174/src/python/pycharm/new_scripts/genomic_coords.csv" # stored in working directory
total_gens = "1000" # total number of generations per simulation 

df = pd.DataFrame()
total_sims = 100 # total number of simulations you want to run


# STEP 3: LOOP TO ADD ROWS WITH ARGUMENTS TO A DATAFRAME
while sim <= total_sims:
  new_row = {'1_python_version' : python_version,
             '2_script': script,
             '3_date' : date,
             '4_sim' : str(int(sim)),
             '5_base_dir' : base_dir,
             '6_u' : u,
             '7_o' : o,
             '8_genomic_coords_path' : genomic_coords_path,
             '9_total_gens': total_gens}
             
  df = df.append(new_row, ignore_index = True) 
  print(df)
  sim = sim + 1

# STEP 4: WRITE THE DATAFRAME (WITHOUT COLUMN HEADERS) TO jobs_wt.txt
with open("jobs_wt.txt", 'a') as f:
    dfAsString = df.to_string(header=False, index=False)
    f.write(dfAsString)

  
  
  
  
  
  
  
  
  
  
  
  
 
