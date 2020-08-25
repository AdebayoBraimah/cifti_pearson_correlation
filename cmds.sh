

# Convert CIFTI to text
wb_command -cifti-convert -to-text dr_stage3_ic0000_tfce_tstat_fwep_c1.dscalar.nii dr_stage3_ic0000_tfce_tstat_fwep_c1.dscalar.txt
wb_command -cifti-convert -to-text REST_agg_Atlas_s4.dtseries.nii REST_agg_Atlas_s4.dtseries.txt

# Convert CIFTI to NIFTI-1
wb_command -cifti-convert -to-nifti dr_stage3_ic0000_tfce_tstat_fwep_c1.dscalar.nii dr_stage3_ic0000_tfce_tstat_fwep_c1.dscalar.nii.gz
wb_command -cifti-convert -to-nifti REST_agg_Atlas_s4.dtseries.nii REST_agg_Atlas_s4.dtseries.nii.gz



