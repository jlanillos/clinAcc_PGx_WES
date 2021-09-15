# clinAc_PGx_WES

This repository includes all self-made scripts written in Python to analyze the pharmacogenetic (PGx) germline genetic variation of 5,001 individuals referred for diagnostic Whole Exome Sequencing (WES).

Mainly, it explains how to translate genotyping data (VCF-like format) into PGx alleles and phenotypes depicted by all the CPIC guidelines (https://cpicpgx.org/guidelines/)

These scripts are ordered according to different steps in the analysis.

#   IMPORTANT COMMENT  
- Individual-associated data (i.e., variation, gender and location) is not avalable due to data protection matters.
- Genotyping data (VCFs) will be available at a public repository upon publication.
- The folder called "testdata" aims to provide with examples of how the input data would look like, with the aim of helping anyone to examine and understand the code

# HOW TO

1) Prepare your input (per sample) VCF files:

- Use BCFTOOLS to split alleles in each individual VCF, then remove duplicated lines (if needed), then split again (just in case some alleles have been merged again) and merge all individual VCFs together.
- Finally, a last splitting alleles of the merged VCF

 
2) We go to Python scripts to generate MERGED_GVCFS.csv --> Alleles.csv --> (merge Alleles.csv from both panels). Alleles.csv contains all resulting PGx alleles for each individual obtained from the VCF files.

`conda activate py36`

The following .CSV file contains all reference alleles defined in CPIC tables which have been used in this study (actionable alleles)
`HAPREF=path/to/test_data/reference_HAPLOTYPES_20201130_hg38_hg19.csv`



## CREATE Alleles.csv file from MERGED VCF file of each panel (CCP17 and SSV6)

### CCP17 panel

`FILES_DIR=/path/to/test_data/CCP17_panel/`

`cd $FILES_DIR`

`for i in *vcf; do mv $i tmp.vcf; bgzip $i; tabix -p vcf $i.gz; rm tmp.vcf.gz;done`

`for i in *vcf.gz; do mv $i tmp.vcf.gz; bcftools norm -m -both tmp.vcf.gz -Oz -o $i; rm tmp.vcf.gz;done`

`for i in *vcf.gz; do mv $i tmp.vcf.gz; bcftools norm -D tmp.vcf.gz -Oz -o $i; rm tmp.vcf.gz; done`

`for i in *vcf.gz; do mv $i tmp.vcf.gz; bcftools norm -m -both tmp.vcf.gz -Oz -o $i; rm tmp.vcf.gz;done`


`realpath *vcf.gz > list.txt`

`bcftools merge -o tmp.vcf.gz -Oz --file-list list.txt`

`bcftools norm -m -both -o MERGED.vcf -Ov tmp.vcf.gz;rm tmp.vcf.gz`

#Create MERGED with gt info:

`python /path/to/scripts/PGx_0_FromMergedVCF2MergedCSV.py --vcf MERGED.vcf`

#Next, create Alleles.csv for this panel (CCP17)

`python /path/to/scripts/PGx_1_merged_to_alleles.py --fileformat '.vcf' --searchpath $FILES_DIR --hapref $HAPREF --mergedCSV MERGED_GVCFs.csv`

### SSV6 panel

`FILES_DIR=/path/to/test_data/SSV6_panel/`

`cd $FILES_DIR`

`for i in *vcf; do mv $i tmp.vcf; bgzip $i; tabix -p vcf $i.gz; rm tmp.vcf.gz;done`

`for i in *vcf.gz; do mv $i tmp.vcf.gz; bcftools norm -m -both tmp.vcf.gz -Oz -o $i; rm tmp.vcf.gz;done`

`for i in *vcf.gz; do mv $i tmp.vcf.gz; bcftools norm -D tmp.vcf.gz -Oz -o $i; rm tmp.vcf.gz; done`

`for i in *vcf.gz; do mv $i tmp.vcf.gz; bcftools norm -m -both tmp.vcf.gz -Oz -o $i; rm tmp.vcf.gz;done`

`realpath *vcf.gz > list.txt`

`bcftools merge -o tmp.vcf.gz -Oz --file-list list.txt`

`bcftools norm -m -both -o MERGED.vcf -Ov tmp.vcf.gz;rm tmp.vcf.gz`

#Create MERGED with gt info

`python /path/to/scripts/PGx_0_FromMergedVCF2MergedCSV.py --vcf MERGED.vcf`

#Next, create Alleles.csv for this panel (SSV6)

`python /path/to/scripts/PGx_1_merged_to_alleles.py --fileformat '.vcf' --searchpath $FILES_DIR --hapref $HAPREF --mergedCSV MERGED_GVCFs.csv`


## Merge data from both panels (SSv6 aned CCP17) [input files are specified inside the script]

#Merge Alleles.csv from both panels

`python /path/to/scripts/PGx_2_mergeAlleles_CCP17_SSV6_together.py`

#Get per sample diplotype information (haplotypes...csv) using Alleles_2....csv) [input files are specified inside the script]

`python /path/to/scripts/PGx_3_mergeSex_Proc_AllelesperInd.py`

#Convert haplotypes to phenotypes [input files are specified inside the script]

`python /path/to/scripts/PGx4_diplotypes2phenotypes.py`


##Please, find the scripts used to plot figures derived from these output files in `/path/to/scripts/Figures_tables`
