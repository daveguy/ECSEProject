registration/ contains the bash scripts used to run ANTs registration. The average brain atlas was registered to each of the T1 images to produce a warp, uing registerAll.sh which calls registration.sh, registration.sh was taken from the ANT scripts and slightly modified. The warp was then applied to the three probability maps with transformAll.sh

SegmentAll.sh is the bash scritp used to segment all the T1 images.

I have provided .py files as well as .ipynb files for jupyter. Either is fine for inspection but I have not tried the .py files to know if they run. If you try to run it (don't really imagine you will) the filenames will need to be modified.

report/ contains the .tex files etc. for my report