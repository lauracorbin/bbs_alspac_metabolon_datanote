# run vol norm through metaboprep
Rscript run_metaboprep_pipeline.R parameter_file_vol_norm_full_dataset.txt
# move output from script dir to release dir
mv figure /REDACTED/working/data/metabolon/step1/metaboprep_release_2025_04_03/vol_norm/metaboprep_release_2025_04_03/
mv metaboprep_Report_v0.md /REDACTED/working/data/metabolon/step1/metaboprep_release_2025_04_03/vol_norm/metaboprep_release_2025_04_03/


# run peak area through metaboprep
Rscript run_metaboprep_pipeline.R parameter_file_peak_area_full_dataset.txt
# move output from script dir to release dir
mv figure /REDACTED/working/data/metabolon/step1/metaboprep_release_2025_04_03/peak_area/metaboprep_release_2025_04_04/
mv metaboprep_Report_v0.md /REDACTED/working/data/metabolon/step1/metaboprep_release_2025_04_03/peak_area/metaboprep_release_2025_04_04/

