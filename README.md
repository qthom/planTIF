# planTIF

All TIF-Seq clusters were called using the script TIF-Seq_isoform_clustering_v6.py and then converted to .bed files
The read pairs of identical positions were first clustered together then a threshold of 20bp was applied to merge slightly different clusters with within the same 20 bp range from both sides
The bed files produced from this clustering were used for the downstream analysis and are available at GSE129523. From those.bed files all following analysis were performed in R

Alternatively:
I produced bed files of all reads (without clustering) found after filtering and are available here in the Github repository.They have been filtered for polyA misprimming and filtered if no match with TSS-Seq clusters was found.

The script "TIF-Seq_pipeline_v2.R" will load the required R packages and download the genomic information for arabidopsis for the rest of the scripts

If any scripts are not provided here they can be found unchanged at https://github.com/Maxim-Ivanov/Nielsen_et_al_2018. from Nielsen et al, 2019

Each script refer to other scripts containing important functions for the analysis.


.
