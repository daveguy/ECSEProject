#!/bin/bash
FILES=BRATS-2/Image_Data/HG/*
GM=Atlas/mni_icbm152_gm_tal_nlin_sym_09a.nii
WM=Atlas/mni_icbm152_wm_tal_nlin_sym_09a.nii
CSF=Atlas/mni_icbm152_csf_tal_nlin_sym_09a.nii
AP="/home/dave/programs/ant/antsbin/bin/"
for f in $FILES
do
	mask=$f/T1_mask.nii.gz
	seg=$$f/T1_seg.nii.gz
	brain=$f/brainS.nii.gz
	full=$f/VSD.Brain.XX.O.MR_T1/VSD.Brain.XX.O.*
	for fixed in $full
	do
		${AP}ThresholdImage 3 $fixed $mask 1 1.e9 #create mask
	
		${AP}ImageMath 3 $mask MD $mask 5 #dilate and shrink to remove holes 
		${AP}ImageMath 3 $mask ME $mask 5

		${AP}ImageMath 3 $mask ClosestSimplifiedHeaderMatrix $mask #needed to simplify header to make origins co-locate
		${AP}ImageMath 3 $brain ClosestSimplifiedHeaderMatrix $fixed

		${AP}Atropos -d 3 -x $mask -c [3,0] -m [0.1,1x1x1] -i kmeans[3] -o $seg -a $brain -v
	done
done

