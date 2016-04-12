#1/bin/bash
FILES=BRATS-2/Image_Data/HG/*
ATLAS=Atlas/mni_icbm152_t1_tal_nlin_sym_09a.nii
for f in $FILES
do
	full=$f/VSD.Brain.XX.O.MR_T1/VSD.Brain.XX.O.*
	out=$f/Transform/
	for fixed in $full
	do
		./registration.sh $fixed $ATLAS forproduction $out
		# printf "fixed: "$fixed"\n""out: "$out"\n"
	done
done