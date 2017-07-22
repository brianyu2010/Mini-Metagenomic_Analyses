######################################
# This script is used to submit
# snakemake files. 
# 
# 1. loads modules
# 2. calls python to generate txt files inside biosample folders
# 3. calls snakemake
# 
# Must first log into singlecell-submit
#######################################

# THESE PARAMETERS NEED TO BE CHANGED EACH RUN
# root_folder=$1

# Generate text files
# source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda
# python generate_snakemake_seed.py $root_folder $subsamples_file"/code_analysis/subsamples.txt"

# Loading modules
module load use.singlecell
module load python/3.4.2

# Calling snakemake files --forcerun split_fastq use -k when actually running
root_folder=/datastore/brianyu/2015.03.11_SygnisTruePrime/
snakemake -j 30 -w 600 -k --config location=$root_folder --cluster "sbatch --job-name={params.name} --ntasks=1 --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}_%j.log" --rerun-incomplete -s Snakefile.toplevel 

# Actually do this
snakemake -j 5 -w 300 -k --config location=$root_folder --cluster "sbatch --job-name={params.name} --ntasks=1 --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} --time={params.time} -o {params.name}_%j.log" --rerun-incomplete -s Snakefile_toplevel_miniMeta.py --unlock -n


echo "Process Completed"
exit

