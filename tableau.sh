#!/bin/bash

# The tblastn command has been configured with -outfmt "7 std qlen qcov(hsp/s) sframe".

# Get bash named parameters (adapted from https://www.brianchildress.co/named-parameters-in-bash/).
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
   fi
  shift
done

cd $wd

if [ $type == "coding" ]
then
#        file=$( echo "${BLAST_OUT}/${query_species}_TRG_blastp_${subject_species}_CDS")

	file=$( echo "${BLAST_OUT}/${query_species}_TRG_multielongated_blastp_${subject_species}_CDS_elongated")
	suff=""
fi
if [ $type == "non-coding" ]
then
#        file=$( echo "${BLAST_OUT}/${query_species}_TRG_tblastn_${subject_species}")

        file=$( echo "${BLAST_OUT}/${query_species}_TRG_multielongated_tblastn_${subject_species}")
	suff="_processed"
fi

# First remove the output files if they already exist.
rm -f ${file}_hits.txt
rm -f ${file}_best_hits.txt


awk -v query_species="${query_species}" -v subject_species="${subject_species}" -v BLAST_OUT="${BLAST_OUT}" -v type="${type}" -v file="${file}" '
BEGIN{
	FS=OFS="\t"
}
# For non-comment lines with a query coverage by subject  >= 50% (the evalue is assumed to have been filtered within the blast command),
$1 !~ /#/ && $14 >= 50 {

        # Get the query name but without the specified frames.
        $1=gensub(/(.*_elongated)_F[0-9]_[0-9]/,"\\1","g",$1)

	# If this is the first time the programme meet a significant aligment for this query and this subject, add them to the "hits.txt" file.
	if ( !($1$2 in data)){
		# Save the pair.
		data[$1$2]=1
		print query_species,$1,subject_species,$2,type >> file"_hits.txt"
	}

	# If this alignment has the best evalue seen for the current query,
        if ( !($1 in evalue) || $11 < evalue[$1] ){
		# Record the evalue,
		evalue[$1]=$11
		# Record the subject.
		best[$1]=$2
	}
}
# At the end, for each query with has a best-subject, add the query-subject pair to the "best_hits.txt" file.
END {
	for (query in best){
		print query_species,query,subject_species,best[query],type >> file"_best_hits.txt"
	}
}
' ${file}${suff}.out
