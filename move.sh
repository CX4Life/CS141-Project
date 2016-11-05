#!/bin/bash
cd WT_RigidityAnalysis/data
find -name '*postPG_BBH.xml' -exec cp {} ../../input \;
find -name '*processed.pdb.knr' -exec cp {} ../../input \;
cd ..
cd ..
cd finalPockets
cp * ../input
cd ../input
for file in *.knr
do
mv "$file" "${file%.knr}.txt"
done
cd ..
