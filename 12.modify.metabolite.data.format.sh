# Code to modify output from metaboprep so it can be read into stata
# run from code locatiom

cd ../../../data/metabolon/step4/alspac/metaboprep_release_2025_04_04/filtered_data/

# need to add in a column header for the sample ID column
# extract header
head -n 1 alspac_only_2025_04_04_Filtered_metabolite_data.txt > ../../data_return/temp_header_line.txt
# extract all but header
tail -n +2 alspac_only_2025_04_04_Filtered_metabolite_data.txt > ../../data_return/temp_metab_data.txt
wc -l ../../data_return/temp_metab_data.txt

# update header
awk 'BEGIN { OFS = "\t" } { $1 = "\"PARENT_SAMPLE_NAME\"\t" $1; print }' ../../data_return/temp_header_line.txt > ../../data_return/temp_header_line2.txt

# add new header back to data
cat ../../data_return/temp_header_line2.txt ../../data_return/temp_metab_data.txt > ../../data_return/temp_metab_data2.txt
wc -l ../../data_return/temp_metab_data2.txt

# remove temp files
rm ../../data_return/temp_header_line.txt
rm ../../data_return/temp_metab_data.txt
rm ../../data_return/temp_header_line2.txt

