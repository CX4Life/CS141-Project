#!/bin/bash
cd WT_RigidityAnalysis/data
find -name '*postPG_BBH.xml' -exec cp {} /home/woodst/Rigidity/input \;
find -name '*processed.pdb.knr' -exec cp {} /home/woodst/Rigidity/input \;
cd ..
cd ..
cd finalPockets
cp * /home/woodst/Rigidity/input
cd ~/Rigidity/input
for file in *.knr
do
mv "$file" "${file%.knr}.txt"
done
cd ..

