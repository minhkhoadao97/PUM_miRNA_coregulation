# PUM_miRNA_coregulation
Codes for identify potential co-regulation between PUM and several miRNAs involved in mouse stem cell development

Identify a list of genes with the following criteria:
* Contain a miRNA binding and a pumilio motif within 300 nucleotide in the 3'UTR
* Conserved (average phyloP score of region between the miRNA site and PUM site >= 2)
* Change in translation efficieny upon PUM knockdown (log 2 fold change > 1)
* RIP-ChIP demonstrate direct PUM binding

Tracks showing conservation and locations of PUM and miRNA sites were created using pygenometracks.

## Workflow

1. Download a phyloP file for mouse genome (mm39 assembly, BigWig format) from ENSEMBL database
2. Run **identify_PUM_miRNA_motifs.py** script
* Input:
  * *ENSEMBL/mouse_3UTR_annotation.gtf.gz* - annotation file obtained from ENSEMBL and processed to only contain 3'UTRs
  * *ENSEMBL/mm39.phyloP35way.bw* (not provided)- pre-computed 35-way phyloP score from ENSEMBL
  * *mir_prediction/* - folder contains prediction of 8 miRNA binding sites obtained from miRGate tool
* Output: output a text file for each miR ("OUTPUT/*_out.csv") showing whether PUM motifs are found proximal to predicted miR sites

3. Run **validate_PUM_miRNA_prediction.ipynb** jupyter notebook
* Input:
  * Output text file from OUTPUT folder
  * *mouse_PUM_target.csv* - RIP-ChIP data for PUM1 and PUM2
  * *translation_analysis_mouse_dev.xls* - Ribosomal profiling data for PUM1/2
* Output:
  * *PUM_miRNA_prediction_result.xlsx*

