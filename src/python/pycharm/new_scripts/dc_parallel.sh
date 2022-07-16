# OBJECTIVE: This shell script runs the jobs described in jobs_dc.txt 
# ARGUMENTS: The number of cores (3, meaning 3 simulations running at once); the list of jobs_dc.txt
# OUTPUT: Runs simulations described in jobs_dc.txt; output are files in working directory/date/output
#------------------------------------------------------------------------------------------------------#

parallel --jobs 3 < jobs_dc.txt