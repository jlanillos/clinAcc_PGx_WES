# These lines of code explain the order of execution to extract all novel variation data for the genes of interest.
SCRIPTSDIR=/path/to/scripts/DISCOVERY_scripts
FILEDIR=/path/to/samples
cd $FILEDIR

GENES=/path/to/file/with/gene/names
OUTDIR=/path/to/output

# Run the scripts in the following order;

python $SCRIPTSDIR/Discovery_PGx_0_getPGxselection.py --searchdir $FILEDIR --fileformat '.tsv' --genes $GENES --outdir $OUTDIR

python $SCRIPTSDIR/Discovery_PGx_1_getTemplate.py --searchpath $FILEDIR

for i in *pgx.csv; do python $SCRIPTSDIR/Discovery_PGx_2_growTemplate.py --file $i --searchpath $FILEDIR; echo $i; done

python $SCRIPTSDIR/Discovery_PGx_3_modifyFinalTemplateCols.py --searchpath $FILEDIR

python $SCRIPTSDIR/Discovery_PGx_4_annotateSamples_updateVarcols.py --searchpath $FILEDIR

# The following lines apply the same script over 5 different genotyping attributes (unique in each individual, called hereinafter as "features")
# The script Discovery_PGx_5_getFeatures_inparallel.py was designed to input a single variable name in order to enable parallelizing its use
# features: genotype,zigosity,VAF,AlleleRatio,AlleleCoverage
python $SCRIPTSDIR/Discovery_PGx_5_getFeatures_inparallel.py --feature genotype --searchpath $FILEDIR &
python $SCRIPTSDIR/Discovery_PGx_5_getFeatures_inparallel.py --feature zigosity --searchpath $FILEDIR &
python $SCRIPTSDIR/Discovery_PGx_5_getFeatures_inparallel.py --feature VAF --searchpath $FILEDIR &
python $SCRIPTSDIR/Discovery_PGx_5_getFeatures_inparallel.py --feature AlleleRatio --searchpath $FILEDIR &
python $SCRIPTSDIR/Discovery_PGx_5_getFeatures_inparallel.py --feature AlleleCoverage --searchpath $FILEDIR &


# Finally, merge all dataframes generated in the previous step.
# Its output file will serve to extract the final results and make Figures/tables related to the novel variation discovery section of clinically actioanble pharmacogenes

python $SCRIPTSDIR/Discovery_PGx_6_mergedAllFeatures.py --searchpath $FILEDIR
