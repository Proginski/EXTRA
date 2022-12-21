#!/bin/bash

awk '
BEGIN{FS=OFS="\t"}

FNR==NR && $0 !~ /#/ && $3 ~ /^mRNA$/ && $9 ~ /ID=/ && $9 ~ /Parent=/ {

	ID=gensub(/.*ID=([^;]+).*/,"\\1","g",$9)
	mRNA_parent=gensub(/.*Parent=([^;]+).*/,"\\1","g",$9)

	parent[ID]=mRNA_parent
}

FNR!=NR { mRNAs[$1] = 1}

END{

for (mRNA in mRNAs){

	if(mRNA in parent){ print mRNA, parent[mRNA] }
	else { print mRNA, mRNA }

}

}' $1 $2


