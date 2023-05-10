# Plakoglobin is a mechanosensitive regulator of naïve pluripotency

This Github repository provides the analysis pipelines used for the publication: 
Timo N. Kohler, Joachim De Jonghe, Anna L. Ellerman *et al.* 2022. “Plakoglobin Is a Mechanoresponsive Regulator of Naïve Pluripotency.” bioRxiv. https://doi.org/10.1101/2022.03.13.484158.

The fastq data is accessible at GEO repository: GSE197643.

The data analysis includes one repository for the bulk RNA-seq analysis and one repository for the single-cell RNA-seq analysis.

Bulk RNA-seq analysis script description
1. Pre-processing the data

Mapping and feature counting scripts can be found in mapping_and_featurecounts.sh

2. downstream analysis using DESeq2

Downstream differential expression analysis between each triplicate can be found in the DESeq2.R script.


Single-cell analysis script description:
1. Pre-processing the data

Note: All the steps can be concatenated in master scripts if automation is needed, the files presented here are intended to be showcased for each step independently.

a) Generating the STAR index

The Mus musculus STAR indexing was ran using the run_STAR_indexer.sh script.

b) bcl2fastq conversion

The BCL files were converted fastq files using the bcl2fastq.sh script.

c) De-multiplexing using the Pheniqs tool from biosails/pheniqs

The run_pheniqs.sh script was used to de-multiplex the fastq files according the run_pheniqs.json json file. The data posted on GSE161947 is already de-multiplexed and corresponds to the data obtained after this step.

d) Running zUMIs to produce count matrices

The zUMIs pipeline was ran using the script run_zumis.sh for each of the samples. To convert the resulting produced dgecounts.rds files to txt matrices counting both intronic and exonic UMI counts using the extract_inex_matrix.R file. The Ensembl names were then converted to gene names using the convert_names.ipynb script.

2. Filtering the data and downstream analysis

The downstream analysis was performed using scanpy and scVelo tools for dimensional reduction, gene expression plotting and latent time analysis. The example script for this analysis is downstream_analysis.ipynb. The differential expression analysis was performed using Seurat and an example script for running the Wilcoxon rank sum test can be found in DGE_Seurat.R. The matrices can be obtained from scanpy by running the follwoing command:
pd.DataFrame(data=adata.raw.X.A, index=adata.obs_names, columns=adata.var_names).T.to_csv("cat_R.csv")
