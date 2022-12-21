#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script takes as input :
    - two lists of CDS that are upstream and downstream a CDS of interest in 
    genomeA
    - two lists of CDS that are upstream and downstream its homologous 
    intergenic region in genomeB
    - a file containning the pairs of orthologs between genomeA (first column),
    and genomeB (second column)
    - two dictionnaries with CDS as keys and their strand as value (one for   
    each genome)
    - an insertion threshold : only insertions whose length exceeds the below 
    threshold are accounted.
    
It aligns the two genomic regions provided base on the CDS of interest and its
homologous IGR position.

It provided insights to determine whether or not the two regions can be 
considered as syntenic.
"""

import pandas as pd
import numpy as np
import re
from dcjoin import calculate_distance 
import string

empty_char= " "
any_CDS_char = "."


def micro_synteny(query_upstream,
                  query_downstream,
                  IGR_upstream,
                  IGR_downstream,
                  ortho_df,
                  strand_dic_A,
                  strand_dic_B,
                  ins_threshold = 0, 
                  nb_flank = 6) :
    
    # Process the previous upstream and downstream lists for display as a dataframe
    query_upstream_to_display = [ query_upstream[i] if i in range(0,len(query_upstream)) else empty_char for i in range(0,nb_flank) ]
    query_downstream_to_display = [ query_downstream[i] if i in range(0,len(query_downstream)) else empty_char for i in range(0,nb_flank) ]
    query_upstream_to_display.reverse()
    query_to_display = query_upstream_to_display +["CDS"]+ query_downstream_to_display
    
    IGR_upstream_to_display = [ IGR_upstream[i] if i in range(0,len(IGR_upstream)) else empty_char for i in range(0,nb_flank) ]
    IGR_upstream_to_display.reverse()
    IGR_downstream_to_display = [ IGR_downstream[i] if i in range(0,len(IGR_downstream)) else empty_char for i in range(0,nb_flank) ]
    IGR_to_display = IGR_upstream_to_display +["IGR"]+ IGR_downstream_to_display
    
    # Display a first dataframe containing the original CDS names
    df = pd.DataFrame()
    df['genomeA']  = query_to_display
    df['genomeB']  = IGR_to_display
    
    #################################
    
    
    
    
    ### Map orthologs on the local alignment ###
    
    # Each pair of orthologs will be associated with a given letter
    alphabet = list(string.ascii_uppercase)
    i=0
    
    query_alpha=[]
    for CDS in query_to_display:
    
        if CDS == "CDS" : query_alpha.append("CDS")
        
        elif CDS != empty_char : query_alpha.append(alphabet[i]); i = i+1
            
        else : query_alpha.append(empty_char)
    
    
    IGR_alpha=[]
    for CDS in IGR_to_display:
    
        if CDS == "IGR" : IGR_alpha.append("IGR")
        
        elif CDS == empty_char : IGR_alpha.append(empty_char)
            
        # Only CDS that have an ortholog in the query's flanking CDS are 
        # associated to the appropriate letter
        elif CDS in list(ortho_df['B']) :

            ortholog = ortho_df['A'][np.where(ortho_df['B'] == CDS)[0][0]]
            if ortholog in query_to_display :

                letter = query_alpha[query_to_display.index(ortholog)]
                
                # CDS which are not on the same strand than their ortholog are
                # represented upside-down
                if strand_dic_A[ortholog] != strand_dic_B[CDS] :
                    # IGR_alpha.append(upsidedown.transform(letter))
                    IGR_alpha.append(letter+"(r)")
                    
                else : IGR_alpha.append(letter)
                    
            else : IGR_alpha.append(any_CDS_char)
            
        else : 
            IGR_alpha.append(any_CDS_char)
    
    # Display a second dataframe containing where orthologs have the same letter
    df_alpha = pd.DataFrame()
    df_alpha['genomeA']  = query_alpha
    df_alpha['genomeB']  = IGR_alpha
    
    #############################################
    
    
    
    ### Analyse the local alignment ###
    
    # Let the sequence of CDS of genomeA be the reference, and the one of genomeB
    # the analysis vector.
    # The analysis vector will be interpreted to determine the number insertions, 
    # deletions and rearrangments events relatively the the reference (genomeA).
    
    # To obtain the analysis vector, remove the IGR symbol and non-orthologs that
    # are before the first ortholog or after the last one
    analysis_vec = IGR_alpha[:]
    analysis_vec.remove("IGR")
    # letter_indexes = [i for i in range(0,len(analysis_vec)) if analysis_vec[i].isalpha()]
    letter_indexes = [i for i in range(0,len(analysis_vec)) if analysis_vec[i].isalpha() or analysis_vec[i][0].isalpha()]
    
    if len(letter_indexes) > 0 :
        analysis_vec = analysis_vec[min(letter_indexes):max(letter_indexes)+1]
        
        # Tranform the analysis vec to a string with all letter normally oriented.
        # analysis_string = [letter if letter != any_CDS_char and upsidedown.transform(letter) not in alphabet[0:nb_flank] else upsidedown.transform(letter) for letter in analysis_vec ]
        analysis_string = [letter if letter != any_CDS_char and "(r)" not in letter not in alphabet[0:nb_flank] else letter[0] for letter in analysis_vec ]
        analysis_string = ''.join(analysis_string)
        
        # Remaining symbols that are not letters are inserted CDS.
        insertions = list(filter(None, re.split("[A-Z]", analysis_string)))
        detected_insetions = [ insertion for insertion in insertions if len(insertion) > ins_threshold ]
        
        # Letters (CDS) from genomeA local sequence that are not in genomeB local 
        # sequence are deletions.
        # deletions = [ letter for letter in query_alpha if letter.isalpha() and all(element not in IGR_alpha for element in [letter,upsidedown.transform(letter)]) and letter != "CDS" ]
        deletions = [ letter for letter in query_alpha if letter.isalpha() and all(letter not in element for element in IGR_alpha) and letter != "CDS" ]
        
        
        # The number of rearrangment events is computed using a customized version of 
        # https://github.com/mlliou112/py-dcj, which implements the algorithm proposed
        # in Braga et al. (2010). For our purpose, it is mainly based in the reference 
        # paper Yancopoulos et al. (2005).
        # Consider only the order and the orientation (straight or upside-down) of
        # the letters in the analysis vector.
        only_letters = [ letter for letter in analysis_vec if letter.isalpha() or letter[0].isalpha() ]
        # only_letters_stranded = [ letter if upsidedown.transform(letter) not in alphabet[0:nb_flank] else upsidedown.transform(letter).lower() for letter in only_letters]
        only_letters_stranded = [ letter if letter in alphabet[0:nb_flank] else letter[0].lower() for letter in only_letters]
        
        
        only_letters_stranded_uniq = []
        previous = ""
        for letter in only_letters_stranded :
            if letter != previous :
                only_letters_stranded_uniq.append(letter)
            else :
                print(letter + " is present more than once in a row in the second genomic region. Maybe another sequence is interspersed between its exons ?")
                
            previous = letter
            
        # Flanking "." represent the ends of a linear sequence.
        # "y" and "z" fix the strand orientation
        genome_a = ['.y'+''.join(alphabet[0:nb_flank])+'z.']
        # genome_b = ['.'+''.join(only_letters_stranded)+'.']
        genome_b = ['.y'+''.join(only_letters_stranded_uniq)+'z.']
        
        print("Number of insertions that exceed the provided threshold : "+str(len(detected_insetions)))
        print("Number of deletions : "+str(len(deletions)))
        
        try :
            nb_rearrangement = calculate_distance(genome_a, genome_b, method="dcj")
            print("Number of rearrangments : "+str(int(nb_rearrangement)))
        except ValueError:
            print("ValueError: Duplicated markers are not allowed.")
            
        ## NEW
        # If "IGR" has at least one letter at each flank :
        IGR_alpha_ol = [ letter for letter in IGR_alpha if letter.isalpha() or letter[0].isalpha() ]
        if "IGR" not in [IGR_alpha_ol[0],IGR_alpha_ol[-1]] : print("Two-sides anchors : TRUE")
        else : print("Two-side anchors : FALSE")
        
    else : print("No orthologs.")
    
    
    
    
    print("\n - ORIGINAL SYNTENY - ")
    print(df)
    print("\n - MAPPED ORTHOLOGS - ")
    print(df_alpha)


