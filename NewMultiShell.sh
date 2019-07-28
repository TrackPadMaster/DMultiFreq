#!/bin/bash

echo "I used to say that I am a test script"
echo "Why have a test branch when everything can be master"
echo "I need a P.txt and Q.txt file to work"
# This whole thing is a hot mess
# First we read the text files into arrays
readarray Pall < ./P.txt
readarray Qall < ./Q.txt
# q is just for simple indexing
echo ${Pall[@]}
echo ${Qall[@]}
q=0
# Cycle through every value of Pall
for p in ${Pall[@]}; do
	# For some infuriating reason, the readarray also saves spaces after the numbers
	# For most functions, this is just fine
	# For bc, the world is ending
	# A quick division by 1 (since they're integers) will fix this
	tempP=$(($p/1))
	tempQ=${Qall[$q]}
	tempQ=$(($tempQ/1))
	q=$(($q+1))
	# "scale" here decides how many digits the ratio will keep
	ratio=$(bc <<< "scale=5;($tempP/$tempQ)")
	# Now we have the ratio that we'll need to save
	awk -v ratio=$ratio '{if($3=="pqratio") {print ratio" !	pqratio"} else {print $0}}' MultiInStart > ./MultiIn0

	# FortranRunner for 20 points, FortranRunner2 for 40 points
	./FortranRunner2.sh
	./DataCollector.sh > ${tempP}over${tempQ}pqratio.dat

done

echo "Hopefully that worked"
