#!/bin/bash

# GO enrichment
#R --slave -f GOstats_analysis.R --args workdir=/Users/shengwei/MEGA/GOstats organism=Staphylococcus_HG001 target_gene_file=HG001_00009__HG001_03232_CopraRNA_result_locustags.tsv ipr2go_file=NZ_CP018205_ipr2go.tsv 

# fetch revigo csv
./fetch_revigo_csv.py HG001_00009__HG001_03232_CopraRNA_result_locustags_GO_enrichment_BP.tsv
./fetch_revigo_csv.py HG001_00009__HG001_03232_CopraRNA_result_locustags_GO_enrichment_CC.tsv
./fetch_revigo_csv.py HG001_00009__HG001_03232_CopraRNA_result_locustags_GO_enrichment_MF.tsv

# revigo visualization
R --slave -f REVIGO_plotter.R --args workdir=/Users/shengwei/MEGA/GOstats REVIGO_BP=HG001_00009__HG001_03232_CopraRNA_result_locustags_GO_enrichment_BP_revigo.csv REVIGO_CC=HG001_00009__HG001_03232_CopraRNA_result_locustags_GO_enrichment_CC_revigo.csv REVIGO_MF=HG001_00009__HG001_03232_CopraRNA_result_locustags_GO_enrichment_MF_revigo.csv out_pdf=HG001_00009__HG001_03232_CopraRNA_result_locustags_GO_enrichment.pdf 

