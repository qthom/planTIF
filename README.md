# planTIF

All TIF-Seq clusters were called using the script TIF-Seq_isoform_clustering_v6.py which also converts them to .bed files directly
The read pairs of identical positions were first clustered together then a threshold of 20bp was applied to merge slightly different clusters within the same 20 bp range from both sides
These raw bed files from after the clustering were filtered for mis primming false positive and kept only if they matched TSS-Seq peaks. These are the bed files used for the downstream analysis and are available at GSE129523. From those.bed files all following analysis were performed in R.

Alternatively: For looking at raw alignments from TIF-Seq;
I produced bed files of all mapped read pair (without clustering) found after filtering and are available here in the Github repository.They have been filtered for polyA misprimming and filtered if no match with TSS-Seq clusters was found. The files are called <*.allreads.bed> and it is possible to load directly into a genome browser

The script "TIF-Seq_pipeline_v2.R" will load the required R packages and download the genomic information for arabidopsis for the rest of the analysis. It will also produce the plots used in the publication. Be carefull to have the files needed for the different scripts in order to repeat the analysis.

If any scripts are not provided here they can be found unchanged at https://github.com/Maxim-Ivanov/Nielsen_et_al_2018. from Nielsen et al, 2019

Before rungin the TIF-Seqp_pipeline_v2.R in R, it is important to read the scripts and to know where every scripts and external files are called from, since each script refers to other scripts containing important functions for the analysis. If the functions  contained in those scripts are not called properly the analysis won't work.
