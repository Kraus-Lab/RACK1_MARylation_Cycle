#!/usr/bin/env python
"""
Script for getting RPKM expression values and differentially expressed genes from Ribo-Seq and RNA-seq data.
Need a "Paths" file indicating the sample names, paths to cs count files, and paths to bam files. 
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
# from plastid import *
# import seaborn as sns
# import pysam
from collections import OrderedDict
import glob
from itertools import combinations
import subprocess, sys, os, datetime, argparse, string

PROGNAME=__file__

# def printer(inp):
#     now = datetime.datetime.now()
#     print >>sys.stdout, "%s [%s]: %s" % (PROGNAME, now, inp)
#     return

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-p","--path",metavar="path to the paths file",
                        help="paths file")
    parser.add_argument("-r","--rna_thres",metavar="threshold count for RNA-seq",
                        help="rna count threshold")
    parser.add_argument("-s","--ribo_thres",metavar="threshold count for Ribo-Seq",
                        help="ribo-seq count threshold")
    parser.add_argument("-a","--appris",metavar="path to APPRIS transcript ID",
                        help="path to APPRIS transcript ID (created by cs generate in Plastid)")
    parser.add_argument("-g","--gene_name",metavar="path to gene names file",
                        help="path to gene names file")

    args = parser.parse_args()
    path_to_paths_file = args.path  # '/work/GCRB/s195535/data/Mauro_CostaMattioli/1_eIF2A_MNase_cortical_RiboSeq/paths.txt'

    rna_thres = int(args.rna_thres)
    ribo_thres = int(args.ribo_thres)

    path_to_appris_transcript_id = args.appris # '/project/GCRB/JChen_lab/shared/ref/mouse/m27/plastid/gencode.vM27.annotation_appris_merged.txt'
    path_to_gene_id_name = args.gene_name # '/project/GCRB/JChen_lab/shared/ref/mouse/m27/gencode.vM27.annotation_genenames.txt'

    # print(path_to_paths_file)

    paths = pd.read_csv(path_to_paths_file, sep='\t', index_col=0)

    appris_transcript_id_conversion = pd.read_csv(path_to_appris_transcript_id, sep='\t', comment="#", header=None).set_index(1)
    gene_id_name = pd.read_csv(path_to_gene_id_name, sep='\t', comment="#", header=None).set_index(0)
    appris_transcript_id_conversion["gene_name"] = appris_transcript_id_conversion[0].apply(lambda x: gene_id_name.loc[x.split('_')[1]])
    appris_transcript_name = appris_transcript_id_conversion.loc[~appris_transcript_id_conversion.index.duplicated(keep='first')]


    rp_cds_rpkm = pd.DataFrame()
    rp_cds_reads = pd.DataFrame()
    rp_utr5_rpkm = pd.DataFrame()
    rp_utr5_reads = pd.DataFrame()
    rp_utr3_rpkm = pd.DataFrame()
    rp_utr3_reads = pd.DataFrame()

    rna_exon_rpkm = pd.DataFrame()
    rna_exon_reads = pd.DataFrame()


    for i, row in paths.iterrows():
        sample = i
        if row["assay"] == "ribo":
            filepath = row["count_path"]
            counts = pd.read_csv(filepath, comment="#", index_col=0, sep="\t")
            
            rp_cds_rpkm[sample] = counts.cds_rpkm
            rp_cds_reads[sample] = counts.cds_reads
            rp_utr5_rpkm[sample] = counts.utr5_rpkm
            rp_utr5_reads[sample] = counts.utr5_reads
            rp_utr3_rpkm[sample] = counts.utr3_rpkm
            rp_utr3_reads[sample] = counts.utr3_reads
            
            
        elif row["assay"] == "rna":
            filepath = row["count_path"]
            counts = pd.read_csv(filepath, comment="#", index_col=0, sep="\t")
            rna_exon_rpkm[sample] = counts.exon_rpkm
            rna_exon_reads[sample] = counts.exon_reads
            

    # We use a threshold of counts to get rid of genes with too low read count
    tmp = rna_exon_reads.replace(0, np.nan).dropna(how = "any", axis = 0) # Get rid of genes with 0 reads or with NaN reads
    rna_exon_reads_thres = tmp[(tmp > rna_thres)].dropna(how = "all", axis = 0) # threshold of 32 reads for RNA-seq

    tmp = rp_cds_reads.replace(0, np.nan).dropna(how = "any", axis = 0)
    rp_cds_reads_thres = tmp[(tmp > ribo_thres)].dropna(how = "all", axis = 0) # threshold of 64 reads for ribosome profiling


    genes_above_thresholdcount = rp_cds_reads_thres.merge(rna_exon_reads_thres, left_index=True, right_index=True, suffixes = ("_rp","_rna")).index

    print(len(genes_above_thresholdcount))

    rna_exon_reads_thres = rna_exon_reads.loc[genes_above_thresholdcount]
    rp_cds_reads_thres = rp_cds_reads.loc[genes_above_thresholdcount]

    rna_exon_rpkm_thres = rna_exon_rpkm.loc[genes_above_thresholdcount]
    rp_cds_rpkm_thres = rp_cds_rpkm.loc[genes_above_thresholdcount]

    rp_utr5_reads_thres = rp_utr5_reads.loc[genes_above_thresholdcount]
    rp_utr5_rpkm_thres = rp_utr5_rpkm.loc[genes_above_thresholdcount]

    rp_utr3_reads_thres = rp_utr3_reads.loc[genes_above_thresholdcount]
    rp_utr3_rpkm_thres = rp_utr3_rpkm.loc[genes_above_thresholdcount]


    rna_exon_reads_thres_genename = appris_transcript_name[["gene_name"]].merge(rna_exon_reads_thres, left_index=True, right_index=True).set_index('gene_name', drop=True)
    rp_cds_reads_thres_genename =  appris_transcript_name[["gene_name"]].merge(rp_cds_reads_thres, left_index=True, right_index=True).set_index('gene_name', drop=True)

    rna_exon_rpkm_thres_genename =  appris_transcript_name[["gene_name"]].merge(rna_exon_rpkm_thres, left_index=True, right_index=True).set_index('gene_name', drop=True)
    rp_cds_rpkm_thres_genename =  appris_transcript_name[["gene_name"]].merge(rp_cds_rpkm_thres, left_index=True, right_index=True).set_index('gene_name', drop=True)

    rp_utr5_reads_thres_genename =  appris_transcript_name[["gene_name"]].merge(rp_utr5_reads_thres, left_index=True, right_index=True).set_index('gene_name', drop=True)
    rp_utr5_rpkm_thres_genename =  appris_transcript_name[["gene_name"]].merge(rp_utr5_rpkm_thres, left_index=True, right_index=True).set_index('gene_name', drop=True)

    rp_utr3_reads_thres_genename =  appris_transcript_name[["gene_name"]].merge(rp_utr3_reads_thres, left_index=True, right_index=True).set_index('gene_name', drop=True)
    rp_utr3_rpkm_thres_genename =  appris_transcript_name[["gene_name"]].merge(rp_utr3_rpkm_thres, left_index=True, right_index=True).set_index('gene_name', drop=True)



    rna_exon_rpkm_thres_genename_avg = pd.DataFrame()
    rp_cds_rpkm_thres_genename_avg = pd.DataFrame()
    rp_utr5_rpkm_thres_genename_avg = pd.DataFrame()
    rp_utr3_rpkm_thres_genename_avg = pd.DataFrame()

    for i in set(paths.query('assay == "ribo"')["condition"].values):
        indices = paths.query('assay == "ribo" and condition == @i').index
        rp_cds_rpkm_thres_genename_avg["RP_"+i] = rp_cds_rpkm_thres_genename[indices].mean(axis=1)
        rp_utr5_rpkm_thres_genename_avg["RP_"+i] = rp_utr5_rpkm_thres_genename[indices].mean(axis=1)
        rp_utr3_rpkm_thres_genename_avg["RP_"+i] = rp_utr3_rpkm_thres_genename[indices].mean(axis=1)
        
    for i in set(paths.query('assay == "rna"')["condition"].values):
        indices = paths.query('assay == "rna" and condition == @i').index
        rna_exon_rpkm_thres_genename_avg["RNA_"+i] = rna_exon_rpkm_thres_genename[indices].mean(axis=1)


    conds = paths[["assay", "condition"]]
    conds['assaycondition'] = conds[['assay', 'condition']].apply(lambda x: '.'.join(x), axis=1)


    conditions = conds["condition"].drop_duplicates().values


    cds_counts_data_rep = pd.DataFrame()
    for sample in rp_cds_reads_thres.columns:
        cds_counts_data_rep[sample] = rp_cds_reads_thres_genename[sample].astype(int)
    for sample in rna_exon_reads_thres.columns:
        cds_counts_data_rep[sample] = rna_exon_reads_thres_genename[sample].astype(int)

    utr5_counts_data_rep = pd.DataFrame()
    for sample in rp_utr5_reads_thres.columns:
        utr5_counts_data_rep[sample] = rp_utr5_reads_thres_genename[sample].astype(int)
    for sample in rna_exon_reads_thres.columns:
        utr5_counts_data_rep[sample] = rna_exon_reads_thres_genename[sample].astype(int)

    utr3_counts_data_rep = pd.DataFrame()
    for sample in rp_utr3_reads_thres.columns:
        utr3_counts_data_rep[sample] = rp_utr3_reads_thres_genename[sample].astype(int)
    for sample in rna_exon_reads_thres.columns:
        utr3_counts_data_rep[sample] = rna_exon_reads_thres_genename[sample].astype(int)

    rna_counts_data_rep = pd.DataFrame()
    for sample in rna_exon_reads_thres.columns:
        rna_counts_data_rep[sample] = rna_exon_reads_thres_genename[sample].astype(int)

    print("Running DESeq2")

    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    import rpy2.robjects.packages as rpackages
    utils = rpackages.importr('utils')
    base = rpackages.importr('base')
    # print(base._libPaths())

    pandas2ri.activate()
    ro.r('library("DESeq2")')
    ro.globalenv['conds'] = conds
    ro.globalenv['data'] = cds_counts_data_rep

    ro.r('dds <- DESeqDataSetFromMatrix(countData=data, colData=conds, design = ~ assaycondition)')
    ro.r('dds <- DESeq(dds, fitType="parametric",  betaPrior=TRUE)')

    res_TE_cds = dict()
    for i, pair in enumerate(list(combinations(conditions, 2))):
        res_TE_cds[pair[1]+'_vs_'+pair[0]] = ro.r('data.frame(results(dds, contrast=(list(c("assayconditionribo.'+pair[1]+'", "assayconditionrna.'+pair[0]+'"), c("assayconditionribo.'+pair[0]+'", "assayconditionrna.'+pair[1]+'")))))')
        res_TE_cds[pair[1]+'_vs_'+pair[0]] = res_TE_cds[pair[1]+'_vs_'+pair[0]].reset_index(drop=True).merge(cds_counts_data_rep.reset_index(), left_index=True, right_index=True).set_index('gene_name')
        res_TE_cds[pair[1]+'_vs_'+pair[0]]["log10padj"] = -np.log10(res_TE_cds[pair[1]+'_vs_'+pair[0]]["padj"])


    pandas2ri.activate()
    ro.r('library("DESeq2")')
    ro.globalenv['conds'] = conds
    ro.globalenv['data'] = utr5_counts_data_rep

    ro.r('dds <- DESeqDataSetFromMatrix(countData=data, colData=conds, design = ~ assaycondition)')
    ro.r('dds <- DESeq(dds, fitType="parametric",  betaPrior=TRUE)')

    res_TE_utr5 = dict()
    for i, pair in enumerate(list(combinations(conditions, 2))):
        res_TE_utr5[pair[1]+'_vs_'+pair[0]] = ro.r('data.frame(results(dds, contrast=(list(c("assayconditionribo.'+pair[1]+'", "assayconditionrna.'+pair[0]+'"), c("assayconditionribo.'+pair[0]+'", "assayconditionrna.'+pair[1]+'")))))')
        res_TE_utr5[pair[1]+'_vs_'+pair[0]] = res_TE_utr5[pair[1]+'_vs_'+pair[0]].reset_index(drop=True).merge(utr5_counts_data_rep.reset_index(), left_index=True, right_index=True).set_index('gene_name')
        res_TE_utr5[pair[1]+'_vs_'+pair[0]]["log10padj"] = -np.log10(res_TE_utr5[pair[1]+'_vs_'+pair[0]]["padj"])

        
    pandas2ri.activate()
    ro.r('library("DESeq2")')
    ro.globalenv['conds'] = conds
    ro.globalenv['data'] = utr3_counts_data_rep

    ro.r('dds <- DESeqDataSetFromMatrix(countData=data, colData=conds, design = ~ assaycondition)')
    ro.r('dds <- DESeq(dds, fitType="parametric",  betaPrior=TRUE)')

    res_TE_utr3 = dict()
    for i, pair in enumerate(list(combinations(conditions, 2))):
        res_TE_utr3[pair[1]+'_vs_'+pair[0]] = ro.r('data.frame(results(dds, contrast=(list(c("assayconditionribo.'+pair[1]+'", "assayconditionrna.'+pair[0]+'"), c("assayconditionribo.'+pair[0]+'", "assayconditionrna.'+pair[1]+'")))))')
        res_TE_utr3[pair[1]+'_vs_'+pair[0]] = res_TE_utr3[pair[1]+'_vs_'+pair[0]].reset_index(drop=True).merge(utr3_counts_data_rep.reset_index(), left_index=True, right_index=True).set_index('gene_name')
        res_TE_utr3[pair[1]+'_vs_'+pair[0]]["log10padj"] = -np.log10(res_TE_utr3[pair[1]+'_vs_'+pair[0]]["padj"])


    conds_RNA = conds.query('assay == "rna"')[["condition"]]

    conditions_RNA = conds_RNA["condition"].drop_duplicates().values

    pandas2ri.activate()
    ro.r('library("DESeq2")')
    ro.globalenv['conds'] = conds_RNA
    ro.globalenv['data'] = rna_counts_data_rep

    ro.r('dds <- DESeqDataSetFromMatrix(countData=data, colData=conds, design = ~ condition)')
    ro.r('dds <- DESeq(dds, fitType="parametric",  betaPrior=TRUE)')

    res_RNA = dict()
    for i, pair in enumerate(list(combinations(conditions_RNA, 2))):
        res_RNA[pair[1]+'_vs_'+pair[0]] = ro.r('data.frame(results(dds, contrast=(c("condition", "' + pair[1] + '", "' + pair[0] + '"))))')
        res_RNA[pair[1]+'_vs_'+pair[0]] = res_RNA[pair[1]+'_vs_'+pair[0]].reset_index(drop=True).merge(rna_counts_data_rep.reset_index(), left_index=True, right_index=True).set_index('gene_name')
        res_RNA[pair[1]+'_vs_'+pair[0]]["log10padj"] = -np.log10(res_RNA[pair[1]+'_vs_'+pair[0]]["padj"])


    print("saving files...")

    writer = pd.ExcelWriter('TE_DE_genes.xlsx', engine='xlsxwriter')

    for pair in res_TE_cds:
        result = res_TE_cds[pair].sort_values('log2FoldChange')
        result.to_excel(writer, sheet_name='CDS_'+ pair)
    for pair in res_TE_utr5:
        result = res_TE_utr5[pair].sort_values('log2FoldChange')
        result.to_excel(writer, sheet_name='UTR5_'+ pair)
    for pair in res_TE_utr3:
        result = res_TE_utr3[pair].sort_values('log2FoldChange')
        result.to_excel(writer, sheet_name='UTR3_'+ pair)
    writer.save()

    writer = pd.ExcelWriter('RNA_DE_genes.xlsx', engine='xlsxwriter')

    for pair in res_RNA:
        result = res_RNA[pair].sort_values('log2FoldChange')
        result.to_excel(writer, sheet_name='RNA_'+ pair)

    writer.save()

    writer = pd.ExcelWriter('RPKM_genes.xlsx', engine='xlsxwriter')

    result = rna_exon_rpkm_thres_genename_avg.sort_index()
    result.to_excel(writer, sheet_name='RNA_exon')

    result = rp_cds_rpkm_thres_genename_avg.sort_index()
    result.to_excel(writer, sheet_name='Ribo_CDS')
    result = rp_utr5_rpkm_thres_genename_avg.sort_index()
    result.to_excel(writer, sheet_name='Ribo_utr5')
    result = rp_utr3_rpkm_thres_genename_avg.sort_index()
    result.to_excel(writer, sheet_name='Ribo_utr3')

    writer.save()

    print("all done")

if __name__=='__main__':
    main()


