#!/bin/bash

# workflow for processing HUMAN RNA-seq that are paired with ribosome profiling:
# 1. QC the reads with fastqc
# 2. Map reads with bowtie to rRNA and contaminants
# 3. Map reads with tophat to transcriptome (.gtf)
# 4. Sort and index bam files
# 5. Make files for IGV (.wig, .tdf)
# 6. Calculate sequencing stats (# aligned reads, etc)
# 7. Run plastid cs script for counts and rpkms
# """

if [ $# -eq 0 ]
then
	echo "Pipeline to process RNA-seq that are matched with Ribosome Profiling datasets, starting from fastq (NOT fastq.gz). Please modify the fastq file name first to get rid of unwanted parts of the name."
	echo "Suitable for RNA-seq libraries (such as the Illumina Tru-seq RNA-seq kit). Make sure to trim the adapters first." 
	echo "Usage: processRibo_RNAseq.sh [options] RNAseq.fastq [RNAseq.fastq]" 
	echo ""
	echo -e "\tOptions:"
	echo -e "\t-----------------"
	# echo -e "\t--trim)\tWhether to trim adapters."
	# echo -e "\t-a)\tAdapter sequence to trim"
	echo -e "\t-f)\tPath to indexed genome fasta file."
	echo -e "\t-r)\tPath to indexed rRNA fasta file."
	echo -e "\t-c)\tPath to indexed contaminant fasta file."
	echo -e "\t-t)\tPath to indexed transcriptome file"
	echo -e "\t-g)\tPath to gtf file"
	echo -e "\t-s)\tPath to .positions file for cs count"
	echo -e "\t-p)\tNumber of processes to use."

else
	##############################################
	# These are default parameters. Should manually change the paths to match your work environment
	# adapter_sequence="CTGTAGGCACCATCAAT" # Usually the adapters are the universal Illumina adapters.
	fasta_path="/project/GCRB/JChen_lab/shared/ref/human/hg19/fasta/hg19"
	rRNA_path="/project/GCRB/JChen_lab/shared/ref/human/rRNA/hs-rrna"
	contaminant_path="/project/GCRB/JChen_lab/shared/ref/contaminants/contaminants"
	transcriptome_index_path="/project/GCRB/JChen_lab/shared/ref/human/hg19/lnc_gencode_transcriptome/lncs_gencode24/lncs_gencode24"
	cs_count_positions_file="/project/GCRB/JChen_lab/shared/ref/human/hg19/others/plastid/gencode.v24lift37.annotation_ucsc_appris_gene.positions"
	gtf_file="/project/GCRB/JChen_lab/shared/ref/human/hg19/lnc_gencode_transcriptome/lncs_gencode24.gtf"
	nproc=6  # number of processes to use for multiprocessing
	##############################################

	# trim=0
	PARAMS=""

	while (( "$#" )); do
	  	case "$1" in
	   	# -a)
	  		# adapter_sequence=$2
	  		# shift 2
	    #   	;;
	   	-f)
	  		fasta_path=$2
	  		shift 2
	      	;;
	   	-r)
	  		rRNA_path=$2
	  		shift 2
	      	;;
	   	-c)
	  		contaminant_path=$2
	  		shift 2
	      	;;
	   	-t)
      		transcriptome_index_path=$2
      		shift 2
      		;;
	   	-p)
      		nproc=$2
      		shift 2
      		;;
      	-s)
      		cs_count_positions_file=$2
      		shift 2
      		;;
  	    -g)
      		gtf_file=$2
      		shift 2
      		;;
      	# --trim)
      	# 	trim=1
      	# 	shift 1
      	# 	;;
		--) # end argument parsing
	  		shift
	  		break
	      	;;
		-*|--*=) # unsupported flags
	  		echo "Error: Unsupported flag $1" >&2
	  		exit 1
	      	;;
		*) # preserve positional arguments
	  		PARAMS="$PARAMS $1"
	  		shift
	      	;;
	  	esac
	done

	echo "Processing Ribo-seq"
	echo "Using fasta: $fasta_path"
	echo "Using rRNA: $rRNA_path"
	echo "Using contaminant: $contaminant_path"
	echo "Using transcriptome index: $transcriptome_index_path"

	eval set -- "$PARAMS"
	FASTQ=$PARAMS
	echo "-----------------"
	echo 

	# Runs QC using FastQC
	mkdir fastqc
	for fq in $FASTQ
	do
		echo "Running QC using FastQC"	
		fastqc $fq -o fastqc &
	done
	wait
	echo 

	mkdir alignments
	mkdir counts
	mkdir igv 

	for fq in $FASTQ
	do
		echo "Running $fq"
		echo "======================================="
		fstem=$(basename "$fq" .fastq)
		# 	echo $fstem

		# if [ $trim -eq 1 ]
		# then
		# # Trim adapter 
		# 	echo "Trimming adapter"
		# 	fq=${fstem}_trimmed.fastq
		# # 	echo $fq
		# 	# zcat $fqgz | fastx_clipper -a $adapter_sequence -ncv -o $fq
		# 	zcat $fqgz | fastx_trimmer -l 50 -Q33 | fastx_reverse_complement -Q33 | fastx_clipper -a $adapter_sequence -nv -Q33 -o $fq

		# 	fastqc $fq -o fastqc # Run fastqc on the trimmed reads
		# else
		# 	echo "Dataset already trimmed..."
		# 	fq=${fstem}.fastq
		# 	zcat $fqgz > $fq
		# fi

		contamfq=${fstem}_contam_unalign.fastq
		echo "Aligning to library contaminants"
		bowtie -v3 -p $nproc $contaminant_path $fq --un alignments/$contamfq >/dev/null

		rrnafq=${fstem}_rrna_unalign.fastq
		echo "Aligning to rRNA"
		bowtie -v3 -p $nproc -S $rRNA_path alignments/$contamfq --un alignments/$rrnafq | samtools view -b - >alignments/${fstem}_hs-rrna_align.bam
		
		echo "Aligning to transcriptome"
		tophatdir=${fstem}_tophat
		tophat --bowtie1 -p $nproc --read-mismatches 2 -g 64 \
		    --no-novel-juncs -T \
		    --transcriptome-index $transcriptome_index_path \
		    -o alignments/$tophatdir \
		    $fasta_path alignments/$rrnafq >alignments/${tophatdir}.out

		# index bam files
		samtools index alignments/$tophatdir/accepted_hits.bam

		# Create soft link and change bam file names to organize all the bam files
		ln -s "$PWD"/alignments/$tophatdir/accepted_hits.bam alignments/${fstem}.bam 
		ln -s "$PWD"/alignments/$tophatdir/accepted_hits.bam.bai alignments/${fstem}.bam.bai

		echo "Making wig and tdf"
		# make .wig file to visualize in IGV
		make_wiggle --count_files alignments/${fstem}.bam  --min_length 25 --max_length 51 --center -o igv/$fstem > /dev/null

		# Convert .wig to .tdf for improved IGV performance
		# igvtools toTDF igv/${fstem}_fw.wig igv/${fstem}_fw.tdf hg19 > /dev/null # genome needs to be changed for other organisms
		# igvtools toTDF igv/${fstem}_rc.wig igv/${fstem}_rc.tdf hg19 > /dev/null # genome needs to be changed for other organisms


		echo "Running cs counts"
		# Use Plastid "cs count" command to get gene level counts and RPKMs
    	cs count --count_files alignments/${fstem}.bam --min_length 25 --max_length 51 --center $cs_count_positions_file counts/${fstem}_counts > /dev/null


    	# echo "Cleaning up temporary files..."
    	# rm alignments/$contamfq

    	echo "Done!"
    	echo "Total number of reads:"
    	echo $(cat $fq|wc -l)/4|bc
    	# echo "Number of reads after trim:"
    	# echo $(cat $fq|wc -l)/4|bc
    	echo "Number of reads after filtering for contaminants:"
    	echo $(cat alignments/$contamfq|wc -l)/4|bc
    	echo "Number of reads after filtering for rRNA:"
    	echo $(cat alignments/$rrnafq|wc -l)/4|bc
    	echo
    	echo

    	rm alignments/$contamfq
    	gzip alignments/$rrnafq
	done

fi

