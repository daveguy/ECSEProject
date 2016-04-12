#!/bin/bash
FILES=BRATS-2/Image_Data/HG/*
GM=Atlas/mni_icbm152_gm_tal_nlin_sym_09a.nii
WM=Atlas/mni_icbm152_wm_tal_nlin_sym_09a.nii
CSF=Atlas/mni_icbm152_csf_tal_nlin_sym_09a.nii
AP="/home/dave/programs/ant/antsbin/bin/"
for f in $FILES
do
	affine=$f/Transform/VSD_fixed_mni_icbm152_t1_tal_nlin_sym_09a_moving_setting_is_forproduction0GenericAffine.mat
	warp=$f/Transform/VSD_fixed_mni_icbm152_t1_tal_nlin_sym_09a_moving_setting_is_forproduction1Warp.nii.gz
	out=$f/Atlas
	full=$f/VSD.Brain.XX.O.MR_T1/VSD.Brain.XX.O.*
	for fixed in $full
	do
		${AP}antsApplyTransforms -d 3 -i $GM -r $fixed -n linear -t $warp -t $affine -o $out/GM_warped.nii.gz
		${AP}antsApplyTransforms -d 3 -i $WM -r $fixed -n linear -t $warp -t $affine -o $out/WM_warped.nii.gz
		${AP}antsApplyTransforms -d 3 -i $CSF -r $fixed -n linear -t $warp -t $affine -o $out/CSF_warped.nii.gz
	done
done