#!/bin/bash

for dir in "CA-GrQc" "CA-CondMat" "Email-EuAll" "Epinions" "slashdot" "dblp" "Amazon0601" "web-Google" "WikiTalk" "as-skitter" "LiveJournal" "uk-2002"
do
	for alg in "NoIt-OE" "EdIt-FE" "EdIt-FV" "EdIt-FVo" "EdIt-FVo2" "EdIt-FC" "EdIt-FCo" "EdIt-OC" "EdIt-OV" "EdIt-rOV"  
	do
		for i in 1 2 3 4 5
		do
			./trilist ${alg} ../data/${dir} >> log.txt
		done
	done
done
