#!/bin/bash
mkdir -p input
tar -xzvf *.gz
./move.sh
cd input/
allPDBs=`ls *.pdb.txt`
allSurfs=`ls *.Surf*`
allXMLs=`ls *.xml`

arrPDB=($allPDBs)
arrSurf=($allSurfs)
arrXML=($allXMLs)

stopLoop=`expr ${#arrPDB[@]} - 1`
cd ..
mkdir -p output
batch=1
for file in `seq 0 $stopLoop`; do
	# echo ${arrPDB[file]}
	# echo ${arrSurf[file]}
	# echo ${arrXML[file]}
	echo "Sending to Python ${arrPDB[file]}"echo ${arrPDB[file]}
	./pipeToPython.sh ${arrPDB[file]} ${arrSurf[file]} ${arrXML[file]}
	if (($batch > 5)); then
		newDir=$(date +%m-%d-%y.%s)
		mkdir -p tmp
		mv 1* 2* 3* tmp
		cd tmp
		for dir in `ls`; do
			cd $dir
			mv * ~/Rigidity/output/
			cd ..
			rm -r $dir
		done
		cd ..
		rm -r tmp
		cd output
		zip $newDir *
		rm *.txt
		cd ..
		let "batch=0"
	fi
	let "batch++"
done
rm -r WT* final*

# echo "Sending 1HE4 to pipeable"
# echo ${proteinName}
# cat < ${proteinName}
# echo "Those are the two ways..."
# echo ${proteinName} | python3 Pipeable.py
echo "Done!"

