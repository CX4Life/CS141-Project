#!/bin/bash
pdb="input/$1"
surf="input/$2"
xml="input/$3"
python3 Pipeable.py $pdb $surf $xml
echo "Python complete!"
now=$(date +%s.%N)
dirName="${pdb:6:4}_$now"
mkdir $dirName
mv *Output.txt ${dirName}
