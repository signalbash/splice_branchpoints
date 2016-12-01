#!bin/sh
mkdir data/inputs/altered_inputs

#transforms gtf to columns
cat data/inputs/gencode.v12.annotation.gtf | grep -v "#" | awk '{print $1, $3, $4, $5, $7, $10, $14, $12,$20}' |sed 's/;//g'| sed 's/"//g' > data/inputs/altered_inputs/gencode.v12.condensed.txt
cat data/inputs/gencode.v24.annotation.gtf | grep -v "#" | awk '{print $1, $3, $4, $5, $7, $10, $14, $12,$20}' |sed 's/;//g'| sed 's/"//g' > data/inputs/altered_inputs/gencode.v24.condensed.txt
