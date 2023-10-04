#!/bin/bash

# This file is batch script used to run commands on the BioHPC cluster.
# The script is submitted to the cluster using the SLURM `sbatch` command.
# Lines starting with # are comments, and will not be run.
# Lines starting with #SBATCH specify options for the scheduler.
# Lines that do not start with # or #SBATCH are commands that will run.

# Name for the job that will be visible in the job queue and accounting tools.
#SBATCH --job-name process_pipeline_RP

# select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH -p 256GB                                             

# Number of nodes required to run this job
#SBATCH -N 1

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 0-5:0:0

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Send an email when the job status changes, to the specfied address.
#SBATCH --mail-type ALL
#SBATCH --mail-user wan-chen.li@utsouthwestern.edu

module purge
module add shared slurm
module load python/3.7.x-anaconda
source activate /project/GCRB/JChen_lab/shared/conda/py36_jupyter_diffexpr/
# module add igvtools/2.3.71 samtools/gcc/1.8 fastx-toolkit/0.0.13.2 bowtie/1.1.2 tophat/gcc/2.1.2 fastqc/0.11.8


# ==============================================================
# EDIT HERE IF NEEDED

path_to_paths_file="paths.txt"
rna_thres=32
ribo_thres=32

path_to_appris_transcript_id="/project/GCRB/JChen_lab/shared/ref/human/hg19/others/plastid/gencode.v24lift37.annotation_ucsc_appris_merged.txt"
path_to_gene_id_name="/project/GCRB/JChen_lab/shared/ref/human/hg19/gencode_transcriptome/gencode.v24lift37.annotation_genenames.txt"

# ==============================================================

echo "Running differential expression"

/project/GCRB/JChen_lab/shared/scripts/runRiboSeq_DiffExp.py -p $path_to_paths_file -r $rna_thres -s $ribo_thres -a $path_to_appris_transcript_id -g $path_to_gene_id_name

echo "Done"




