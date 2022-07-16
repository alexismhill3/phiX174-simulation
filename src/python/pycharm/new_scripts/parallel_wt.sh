# OBJECTIVE: This shell script runs the jobs described in jobs_wt.txt 
# ARGUMENTS: The number of cores (80, meaning 80 simulations running at once); the list of jobs_wt.txt
# OUTPUT: Runs simulations described in jobs_wt.txt; output are files in working directory/date/output
#------------------------------------------------------------------------------------------------------#

parallel --jobs 80 < jobs_wt.txt