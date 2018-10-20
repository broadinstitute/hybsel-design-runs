#!/bin/bash

# NOTE: The dengue.fasta sequences in here contain duplicates -- it is a copy prior
# to commit 14b613c in the hybsel_design repo. It has 4083 Dengue sequences, but
# only 3785 are unique. After the download_dataset_fastas.py script is run,
# duplicates are removed; the duplicates in here won't affect the end output.

python ../../fasta-processing/split_by_regex.py -i dengue.fasta -r ".*(virus 1[ ,]|type 1[ ,]|Subtype:1|DENV-1).*" -o1 ~/tmp/dengue_1.fasta -o2 ~/tmp/tmp-dengue-non1.fasta

python ../../fasta-processing/split_by_regex.py -i tmp-dengue-non1.fasta -r ".*(virus 2[ ,]|type 2[ ,]|Subtype:2|DENV-2).*" -o1 ~/tmp/dengue_2.fasta -o2 ~/tmp/tmp-dengue-non2.fasta

python ../../fasta-processing/split_by_regex.py -i tmp-dengue-non2.fasta -r ".*(virus 3[ ,]|type 3[ ,]|Subtype:3|DENV-3).*" -o1 ~/tmp/dengue_3.fasta -o2 ~/tmp/tmp-dengue-non3.fasta

python ../../fasta-processing/split_by_regex.py -i tmp-dengue-non3.fasta -r ".*(virus 4[ ,]|type 4[ ,]|Subtype:4|DENV-4).*" -o1 ~/tmp/dengue_4.fasta -o2 ~/tmp/tmp-dengue-non4.fasta

# Now hardcode types for the ones that got missed
python ../../fasta-processing/split_by_regex.py -i tmp-dengue-non4.fasta -r ".*(KR052012|KP772252|KF971871|KF971870|KF971869|JX669466|JX669465|JX669464|JX669463|JX669462|JX669461|HG316482|HG316481|KC759167|CS477306|CS477305|CS477304|CS477265|CS477264|CS477263|A75711|CS805347|CS805343|CS479204|CS479203).*" -o1 dengue_1_missed.fasta -o2 /dev/null

python ../../fasta-processing/split_by_regex.py -i tmp-dengue-non4.fasta -r ".*(KF479233|JX669479|JX669478|JX669477|JX669476|AF489932|CS479202|CS479167|CS479165|CS805348|CS805344|CS477302).*" -o1 dengue_2_missed.fasta -o2 /dev/null

python ../../fasta-processing/split_by_regex.py -i tmp-dengue-non4.fasta -r ".*(KF954949|KF954948|KF954947|KF954946|KF954945|JF808120|HG316484|HG316483|KC261634|CS805342|CS805345|CS479205).*" -o1 dengue_3_missed.fasta -o2 /dev/null

python ../../fasta-processing/split_by_regex.py -i tmp-dengue-non4.fasta -r ".*(KF907503|KC333651|CS805346|CS479206).*" -o1 dengue_4_missed.fasta -o2 /dev/null

mv dengue_1.fasta dengue_1_hit.fasta
mv dengue_2.fasta dengue_2_hit.fasta
mv dengue_3.fasta dengue_3_hit.fasta
mv dengue_4.fasta dengue_4_hit.fasta

cat dengue_1_hit.fasta dengue_1_missed.fasta > dengue_1.fasta
cat dengue_2_hit.fasta dengue_2_missed.fasta > dengue_2.fasta
cat dengue_3_hit.fasta dengue_3_missed.fasta > dengue_3.fasta
cat dengue_4_hit.fasta dengue_4_missed.fasta > dengue_4.fasta

rm tmp-dengue-non1.fasta
rm tmp-dengue-non2.fasta
rm tmp-dengue-non3.fasta
rm tmp-dengue-non4.fasta
