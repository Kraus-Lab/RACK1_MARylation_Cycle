#!/bin/bash

# This file is batch script used to run commands on the BioHPC cluster.
# The script is submitted to the cluster using the SLURM `sbatch` command.
# Lines starting with # are comments, and will not be run.
# Lines starting with #SBATCH specify options for the scheduler.
# Lines that do not start with # or #SBATCH are commands that will run.

# Name for the job that will be visible in the job queue and accounting tools.
#SBATCH --job-name process_pipeline_RNA

# select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH -p 256GB                                             

# Number of nodes required to run this job
#SBATCH -N 1

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 0-20:0:0

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Send an email when the job status changes, to the specfied address.
#SBATCH --mail-type ALL
#SBATCH --mail-user wan-chen.li@utsouthwestern.edu


module purge
module add shared slurm
module load python/2.7.x-anaconda
source activate /project/GCRB/JChen_lab/shared/conda/py27riboseq/
module add igvtools/2.3.71 samtools/gcc/1.8 fastx-toolkit/0.0.13.2 bowtie/1.1.2 tophat/gcc/2.1.2 fastqc/0.11.8

# ==============================================================
# EDIT HERE IF NEEDED
adapter_sequence="CTGTAGGCACCATCAAT"  # Usually the adapters are the universal Illumina adapters.

fasta_path="/project/GCRB/JChen_lab/shared/ref/human/hg19/fasta/hg19"
rRNA_path="/project/GCRB/JChen_lab/shared/ref/human/rRNA/hs-rrna"
contaminant_path="/project/GCRB/JChen_lab/shared/ref/contaminants/contaminants"
transcriptome_index_path="/project/GCRB/JChen_lab/shared/ref/human/hg19/lnc_gencode_transcriptome/lncs_gencode24/lncs_gencode24"
cs_count_positions_file="/project/GCRB/JChen_lab/shared/ref/human/hg19/others/plastid/gencode.v24lift37.annotation_ucsc_appris_gene.positions"
gtf_file="/project/GCRB/JChen_lab/shared/ref/human/hg19/lnc_gencode_transcriptome/lncs_gencode24.gtf"
nproc=6  # number of processes to use for multiprocessing



# ==============================================================


# ==============================================================
# Demultiplex raw fastq files according to linker barcodes

trim() {
    local fqgz=$1
    echo "Running $fqgz"
    fstem=$(basename "$fqgz" _R1_001.fastq.gz)
    # zcat $fqgz | fastx_clipper -a $adapter_sequence -ncv \
    # | fastx_barcode_splitter.pl --eol --bcfile $bcfile --prefix ${fstem}_ --suffix "_disamb.fastq"
    # zcat $fqgz |  fastx_trimmer -l 50 -Q33 | fastx_clipper -a $adapter_sequence -ncv -Q33 \
    # | fastx_barcode_splitter.pl --eol --bcfile $bcfile --prefix ${fstem}_ --suffix "_demult.fastq"

    fq=${fstem}_trimmed.fastq
    zcat $fqgz | fastx_trimmer -l 50 -Q33 | fastx_reverse_complement -Q33 | fastx_clipper -a $adapter_sequence -nv -Q33 -o $fq


}

echo "Trimming..."
for fqgz in *.fastq.gz
do
    trim $fqgz & 
done

wait

# ==============================================================

mkdir nohup

echo "Running processRibo_RNAseq...."
for fq in *_trimmed.fastq
do
    fstem=$(basename "$fq" _trimmed.fastq)
    echo "Running $fq"
    nohup /project/GCRB/JChen_lab/shared/scripts/processRibo_RNAseq.sh $fq \
    -f $fasta_path -r $rRNA_path -c $contaminant_path -t $transcriptome_index_path -g $gtf_file -s $cs_count_positions_file -p $nproc \
    > nohup/${fstem}.out &
done

wait

gzip *.fastq

echo "Done!"
