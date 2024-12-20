# clinAc_PGx_WES

This repository includes all self-made scripts written in Python to analyze the pharmacogenetic (PGx) germline genetic variation of 5,001 individuals referred for diagnostic 
Whole Exome Sequencing (WES), published here:

https://www.nature.com/articles/s41525-022-00283-3


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

 
2) We go to Python3 scripts to generate MERGED_GVCFS.csv --> Alleles.csv --> (merge Alleles.csv from both panels). Alleles.csv contains all resulting PGx alleles for each individual obtained from the VCF files.

`conda activate py36`

The following .CSV file contains all reference alleles defined in CPIC tables which have been used in this study (actionable alleles)
`HAPREF=path/to/test_data/Ref_table.csv`



## CREATE Alleles.csv file from MERGED VCF file (both CCP17 and SSV6 panels, separately)
Taking Ref_table.csv file and merged VCF file with all variant calling results of the whole cohort, the aim is to create Alleles.csv. This file contains genotyping data relevant to resolve each pharmacogenetic allele. First, work with either CCP17 or SSV6 data (CES or WES data), and then, merge both panel into a final Alleles.csv file.

Several "normalization steps" are applied using bcftools in order to split multiallelic sites. These steps make easier to parse merged vcf into csv file

NOTE: individual VCF files can be found at https://www.ebi.ac.uk/ena/browser/view/PRJEB48632 

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


## Merge data from both panels (SSv6 and CCP17, a.k.a WES and CES, respectively) [input files are specified inside the script]

#Merge Alleles.csv from both panels

`python /path/to/scripts/PGx_2_mergeAlleles_CCP17_SSV6_together.py`

#Get per sample diplotype information (haplotypes.csv) using Alleles.csv. [input files are specified inside the script]
haplotypes.csv contains per-sample pharmacogenetic allele information

`python /path/to/scripts/PGx_3_mergeSex_Proc_AllelesperInd.py`
(This script is important to correct gene-specific cases which require to remove redundant alleles. For example, if an alleleX is defined by variant A and alleleY is defined by phased variants A and B, the allele with the lowest number of variants is filtered out. All positive cases in our data have been revised and included in this script. This script also adds indels which have been manually revised to avoid missing some variants due to variant caller annotation diversity)


#Convert haplotypes to phenotypes [input files are specified inside the script]

`python /path/to/scripts/PGx4_diplotypes2phenotypes.py`


##Please, find the scripts used to plot figures derived from these output files in `/path/to/scripts/Figures_tables`
