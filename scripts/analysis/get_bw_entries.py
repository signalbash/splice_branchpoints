#!/usr/bin/python
## Get wanted entries out of a bigwig file
##	Copyright (c) B. Signal 2016 (b.signal@garvan.org.au)

import argparse
import pandas 
import pyBigWig

#Load in files
parser = argparse.ArgumentParser()
parser.add_argument('-bw','--bw_input',help='Bigwig file name', required=True)
parser.add_argument('-csv','--input_query',help='CSV file with genomic regions', required=True)
parser.add_argument('-o','--output_file',help='Output file name', required=True)

args = parser.parse_args()
 
bw = pyBigWig.open(args.bw_input)

infile = open(args.input_query)

cons=[]

##Goes through each line in the query file
for line in infile:

	chrom=str(line.split(",")[0].strip())
	start=int(line.split(",")[1].strip())-1
	end=int(line.split(",")[2].strip())

	cons.append(bw.values(chrom,start,end))


panda_df=pandas.DataFrame(cons)
panda_df.to_csv(args.output_file, index=False,header=False)
