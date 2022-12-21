#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 12:19:21 2022

@author: paul.roginski

The idea behind this script is that all possible most recent common ancestors 
of : a given species on the first hand, and its neighbors on the other hand,
are simply all the (node-represented) ancestors of the given species.
AND, all the ancestors-nodes of a given species, can be order by counting the 
number of edges between them and the root.
"""


import argparse
import dendropy
import pandas as pd
import math

# Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("-tree", required=True, help="newick file for the input phylogeny tree")
parser.add_argument("-focal", required=True, help="name of the focal species")
parser.add_argument("-names", required=True, help="tsv file with focal,cds,species,seq,type")
parser.add_argument("-out", required=False, help="output file name")
args = parser.parse_args()

  
# newick file for the phylogeny tree
tree_file = args.tree

#tree_file = "/home/paul.roginski/Bureau/Saccharomyces_species.nwk"

tree = dendropy.Tree.get(path=tree_file,
                         schema='newick',
                         preserve_underscores=True,
                         rooting='force-rooted')



#print(tree.as_ascii_plot())

#test="/home/paul.roginski/Bureau/rna-NM_001178926.3_homologs.tsv"
#homologs_df = pd.read_csv(test, sep='\t', header = None)
homologs_df = pd.read_csv(args.names, sep='\t', header = None)
homologs_df.columns =['focal', 'cds', 'species', 'seq', 'type']
taxon_labels = homologs_df.species.tolist()
coding_labels = homologs_df.loc[homologs_df.type == "coding"].species.tolist()
noncoding_labels = homologs_df.loc[homologs_df.type == "non-coding"].species.tolist()

# Retrieve in the taxon_namespace attribute of the tree object, the element 
# corresponding with the focal species.
#focal_name="Scer_NCBI"
focal_name = args.focal
print("focal name : {}".format(focal_name))
print("names : {}".format(taxon_labels))
index = [ i for i in range(0,len(tree.taxon_namespace)) if tree.taxon_namespace[i].label == focal_name ][0]
focal = tree.taxon_namespace[index]


# Set every edge length to 1. By doing so, we are sure the calc_node_root_distances can be used
# Besides, we do not care about the actual value of each edge lentgh in this case.
for node in tree.preorder_node_iter():
    node.edge.length = 1

# Adds attribute “root_distance” to each node, with value set to the sum of 
# edge lengths from the node to the root. Returns list of distances. 
tree.calc_node_root_distances(return_leaf_distances_only=False) 

# For each taxon provided in names, get its most recent common ancestor with 
# the focal species (= a node), and get its distance to the tree root.
root_distance = { name:tree.mrca(taxon_labels=[focal_name, name]).root_distance for name in taxon_labels if name != focal_name }
print("root_distance : {}".format(root_distance))

print("coding_labels : {}".format(coding_labels))
coding_root_distance = { label:root_distance[label] for label in coding_labels if label != focal_name }
print("coding root distance : {}".format(coding_root_distance))

if len(coding_root_distance) > 0 :
    min_coding_root_distance = min(coding_root_distance.values())
else :
    min_coding_root_distance = math.inf

print("noncoding_labels : {}".format(noncoding_labels))
# There MUST be a non-coding homolog farer than the farest coding homolog.
noncoding_root_distance = { label:root_distance[label] for label in noncoding_labels if label != focal_name if root_distance[label] < min_coding_root_distance}
print("noncoding root distance : {}".format(noncoding_root_distance))

if len(noncoding_root_distance) > 0 :
    max_coding_root_distance = max(noncoding_root_distance.values())

    out_group = [ name for name in root_distance.keys() if root_distance[name] == max_coding_root_distance ] 

    out_group_rows = homologs_df.loc[(homologs_df.type == "non-coding") & (homologs_df.species.isin(out_group))][["cds", "species", "seq"]]
    detailed_seq = out_group_rows.seq.str.rsplit("_", n = 2, expand = True)
    out_group_rows.seq = detailed_seq[0]
    out_group_rows['start'] = detailed_seq[1]
    out_group_rows['stop'] = detailed_seq[2]
            
    if args.out : 
        out_group_rows.to_csv(args.out, mode='w', index=False, header=False, sep="\t")
        
else :
    print("There is no outgroup.")
    
    
