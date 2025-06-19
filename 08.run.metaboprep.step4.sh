# copy over feature metadata to each folder
cp ../../../data/metabolon/step3/metaboprep_release_2025_04_04/vol_norm/metaboprep_release_2025_04_04/filtered_data/bbs_plus_alspac_no_outliers_2025_04_04_Filtered_feature_data.txt ../../../data/metabolon/step4/alspac/
cp ../../../data/metabolon/step3/metaboprep_release_2025_04_04/vol_norm/metaboprep_release_2025_04_04/filtered_data/bbs_plus_alspac_no_outliers_2025_04_04_Filtered_feature_data.txt ../../../data/metabolon/step4/bbs/
cp ../../../data/metabolon/step3/metaboprep_release_2025_04_04/vol_norm/metaboprep_release_2025_04_04/filtered_data/bbs_plus_alspac_no_outliers_2025_04_04_Filtered_feature_data.txt ../../../data/metabolon/step4/bbs_alspac/

# run metaboprep on combined dataset
Rscript run_metaboprep_pipeline.R parameter_file_get_sumstats.txt
# move extra files to release
mv figure /REDACTED/working/data/metabolon/step4/bbs_alspac/metaboprep_release_2025_04_04/
mv metaboprep_Report_v0.md /REDACTED/working/data/metabolon/step4/bbs_alspac/metaboprep_release_2025_04_04/

# run metaboprep on alspac dataset
Rscript run_metaboprep_pipeline.R parameter_file_get_sumstats_alspac.txt
# move extra files to release
mv figure /REDACTED/working/data/metabolon/step4/alspac/metaboprep_release_2025_04_04/
mv metaboprep_Report_v0.md /REDACTED/working/data/metabolon/step4/alspac/metaboprep_release_2025_04_04/

# run metaboprep on bbs dataset
Rscript run_metaboprep_pipeline.R parameter_file_get_sumstats_bbs.txt
# move extra files to release
mv figure /REDACTED/working/data/metabolon/step4/bbs/metaboprep_release_2025_04_07/
mv metaboprep_Report_v0.md /REDACTED/working/data/metabolon/step4/bbs/metaboprep_release_2025_04_07/
