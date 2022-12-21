#!/bin/bash

# For a given CDS, this script looks for homologs among the CDS and intergenic regions (i.e. IGR) (BLASTp) of the neighbor species.
# It generates :
# - a FASTA file with the CDS (nucl) seq and the homologs (CDS or IGR)
# - a txt file with the names (species) and types (CDS or IGR) present in the FASTA file
# - a nwk file with evolutionnary relationships bewteen the selected species.

# Get bash named parameters (adapted from https://www.brianchildress.co/named-parameters-in-bash/)
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
   fi
  shift
done


trap "exit 1" TERM
export TOP_PID=$$


cd $main_dir
echo $query


# Since the query's name will be used to name most of the following files,
# some special characters in it can stop the pipeline before the end.
# At this step, the query's name as already been corrected.
# Retrieve the original name.
query_uncorr=$(echo $query | sed "s/__COLON__/:/g")

echo $query_uncorr


# This function check the variables $header and $seq and write them to ${process_dir}/${query}_toalign.fna
add_to_fasta () {

        if [ -z $header ]; then
                echo "*** NO HEADER for the sequence : $seq ***"
                kill -s TERM $TOP_PID
        else 
                if [ -z $seq ]; then
                        echo "*** EMPTY SEQUENCE for the header $header . This EXACT header might be absent of the file you are parsing ***"
                        kill -s TERM $TOP_PID
                else
                        echo $header >> ${process_dir}/${query}_toalign.fna
                        echo $seq >> ${process_dir}/${query}_toalign.fna
                fi
        fi

        header=""
        seq=""
}


# Append the focal (species of interest) short name to the output FASTA file and to the output txt file
header=$(awk -F"," -v focal=$focal_species '$1 == focal {print ">"$2}' ${process_dir}/$phylip_names)
awk -F"," -v focal=$focal_species '$1 == focal {print $2" CDS"}' ${process_dir}/$phylip_names > ${process_dir}/${query}_names_n_types.txt

# Add the focal (query) nucleotide sequence to the FASTA file
# First line : make sur to deal with a linearized (long format) input FASTA
# Second line : get the sequence without the header
seq=$(cat CDS/${focal_species}_CDS.fna | awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}' | \
grep -x ">${query_uncorr}" -m1 -A1  | grep -v \>)

# Safely append header and seq to the output FASTA
add_to_fasta


# Then search for homologs in the neighbor sequences for each neighbor species.
cat GENOMES/$species_list | grep -v -x $focal_species | while read neighbor; do

	echo $neighbor

	# In the BLASTp output of focal vs neighbor :
        # Only keep one row per query ($1) base on minimal evalue ($11)
        # Only keep this row if it satisfies the evalue ($11) threshold and print the subject sequence ($2) AND the match is on a POSITIVE strand (otherwise it cannot be a match with the tested neighbor's CDS) and with a qcov >= 0.5 
	CDS_sbjct=$(awk -F"\t" -v query="${query_uncorr}" '
        
                          BEGIN { pattern = "^"query"_elongated$" }                                                                                 
        
                          $2 ~ pattern {print $4 ; exit }
        
                          ' BLAST_OUT/${focal_species}_TRG_multielongated_blastp_${neighbor}_CDS_elongated_best_hits.txt )

	# If there is a CDS match
	if [ ! -z $CDS_sbjct ]
	then

		echo "CDS match found : $CDS_sbjct"

		# Append the neighbor short name to the output FASTA file and to the output txt file
                header=$(awk -F"," -v neighbor=$neighbor '$1 == neighbor {print ">"$2}' ${process_dir}/$phylip_names)
                awk -F"," -v neighbor=$neighbor '$1 == neighbor {print $2" CDS"}' ${process_dir}/$phylip_names >> ${process_dir}/${query}_names_n_types.txt
                
		# Only add the sequence not the corresponding header.
		seq=$(grep -x -E ">${CDS_sbjct}" CDS/${neighbor}_CDS_elongated.fna -m1 -A1  | grep -v ">")

		# Safely append header and seq to the output FASTA
		add_to_fasta


	else
		echo "NO CDS match, looking for a non-coding match..."
	        # In the tBLASTn output of focal vs neighbor :
	        # Only keep one row per query ($1) base on minimal evalue ($11)
        	# Only keep this row if it satisfies the evalue ($11) threshold and print the subject sequence ($2) and with a qcov >= 0.5 
        	IGR_sbjct=$(awk -F"\t" -v query="${query_uncorr}" '
	
			  BEGIN { pattern = "^"query"_elongated$" }                                            				            
	
			  $2 ~ pattern {print $4 ; exit }
	
			  ' BLAST_OUT/${focal_species}_TRG_multielongated_tblastn_${neighbor}_best_hits.txt  )


	        # If there is an IGR match
	        if [ ! -z $IGR_sbjct ]; then

        	        echo "IGR match found : $IGR_sbjct"
        	        
        	        sbjct=$(awk '{print gensub(/(.*)_([0-9]+)_([0-9]+)/,"\\1","g",$0)}' <(echo "${IGR_sbjct}") )
        	        sstart=$(awk '{print gensub(/(.*)_([0-9]+)_([0-9]+)/,"\\2","g",$0)}' <(echo "${IGR_sbjct}") )
        	        send=$(awk '{print gensub(/(.*)_([0-9]+)_([0-9]+)/,"\\3","g",$0)}' <(echo "${IGR_sbjct}") )
        	        
        	        echo "sbjct $sbjct sstart $sstart send $send"
        	        
        	        IGR_line=$(awk -F"\t" -v query="${query_uncorr}" -v sbjct="${sbjct}" -v sstart="${sstart}" -v send="${send}" ' 
        	        		
        	        		BEGIN { pattern = "^"query"_elongated_F[0-9]_[0-9]$" } 
        	        		
        	        		$1 ~ pattern && $2 ~ sbjct && $9 ~ sstart && $10 ~ send { print ; exit }
        	        		
        	        		' BLAST_OUT/${focal_species}_TRG_multielongated_tblastn_${neighbor}_processed.out)
        	        		
        	        echo "IGR_line : $IGR_line"

			# Append the neighbor short name to the output FASTA file and to the output txt file
                        header=$(awk -F"," -v neighbor=$neighbor '$1 == neighbor {print ">"$2}' ${process_dir}/$phylip_names)
                        awk -F"," -v neighbor=$neighbor '$1 == neighbor {print $2" IGR"}' ${process_dir}/$phylip_names >> ${process_dir}/${query}_names_n_types.txt
			
			# Elongate the hit as suggested by Papadopoulos and print the result as a one-line GFF.
			elongated_hit=$( awk 'BEGIN{FS=OFS="\t"} {
			
				$2=gensub(/(.*)_([0-9]+)_([0-9]+)/,"\\1","g",$2)
			
				qlen=$13
				qs=$7
				qe=$8
				ss=$9
				se=$10

				q_missing_start = qs - 1
				q_missing_end = qlen - qe

				q_ms_nucl = q_missing_start*3
				q_me_nucl = q_missing_end*3

				if(ss < se){ print $2,"ELONGATED_HIT","non_coding", ss-q_ms_nucl , se+q_me_nucl ,".","+",".","X"}
				else{ print $2,"ELONGATED_HIT","non_coding", se-q_me_nucl , ss+q_ms_nucl+1 ,".","-",".","X"}
				
			}' <(echo "$IGR_line")
			)
			
			echo "elongated hit : ${elongated_hit}"
			
			# If the start coordinate of the elongated hit is inferior to one, or if its end coordinate is greter than the whole sequence,
			# bedtools getfasta will be in error : so correct it. 
			elongated_hit_corr=$( awk 'BEGIN{FS=OFS="\t"}{
						
						# For the first file, store the seq as key and its length as value 
						if(NR==FNR){long[$1]=$2}
						
						# For the second file, if the start coordinate is inferior to one, or if the end coordinate is greater 
						# than the whole sequence, bedtools getfasta will be in error : so correct it.  
						else{
							if($4 < 1){$4=1}
							if($5 > long[$1]){$5 = long[$1]}
							print $0
						}
						}' GENOMES/${neighbor}.genome <(echo "$elongated_hit")
					)

			echo "elongated (corrected) hit : ${elongated_hit_corr}"
						
			# Get the corresponding sequence and append it to the alignment FASTA | wihtout the header (already written, see above)
			# (bedtools getfasta does generates long lines FASTA files)
			seq=$(bedtools getfasta -fi GENOMES/${neighbor}.fna -bed <(echo "$elongated_hit_corr") -s | grep -v ">")
			
			# Safely append header and seq to the output FASTA
			add_to_fasta


	        else
	                echo "NO IGR match."
	        fi
	fi
done


# Only retain cases where at least 3 sequences (1 focal + 2 neighbors) have been collected
seq_nb=$(cat ${process_dir}/${query}_names_n_types.txt | wc -l)

if [ $seq_nb -gt 2 ]
then
	echo "OK enough (${seq_nb}) sequences"
	names=$(awk '{print $1}' ${process_dir}/${query}_names_n_types.txt)
	echo "names = $names"
	# Generate the phylogenetic tree with only the species that have a match
	python /SCRIPTS/get_subtree.py -tree ${process_dir}/${tree} -names <(echo "$names") -out ${process_dir}/${query}.nwk
else
	echo "NOT enough (${seq_nb}) sequences"
	mv ${process_dir}/${query}_toalign.fna ${process_dir}/fail
	mv ${process_dir}/${query}_names_n_types.txt ${process_dir}/fail_names_n_types
fi
