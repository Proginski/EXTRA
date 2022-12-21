#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 12:19:21 2022

@author: paul.roginski
"""

import argparse
import dendropy


# Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("-names", required=True, help="list of names to retain", default=False)
parser.add_argument("-tree", required=True, help="newick file for the input phylogeny tree")
parser.add_argument("-out", required=False, default="./pruned_tree.nwk",help="output newick file")
parser.add_argument("-prune_extra_taxa", required=False, default="False", help="If some taxa in the tree are not in the names file, whther or not to adjust (prune) the tree to the names file and continue.")
args = parser.parse_args()

# argsnames = "/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/SCER_NCBI/test"
with open(args.names) as f:
    names = [name.strip() for name in f]
    
# newick file for the phylogeny tree
tree_file = args.tree


# names = ["Saccharomyces cerevisiae","Saccharomyces paradoxus","Saccharomyces bayanus"]
# tree_file = "/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/ORFdate/liste_carvunis.nwk"
# names = ["Scer_NCBI","Spar","Sbay"]
# tree_file = "/home/me/STORE/Eukaryotes/PHYLO/SCER_NCBI/GENOMES/Saccharomyces_species_corr_doubleoutput.nwk"


tree = dendropy.Tree.get(path=tree_file, schema='newick', preserve_underscores=True)

tree_taxa = [ tree.taxon_namespace[i].label for i in range(0,len(tree.taxon_namespace)) ]
names_not_in_tree_taxa = [ name for name in names if name not in tree_taxa]
tree_taxa_not_in_names = [ label for label in tree_taxa if label not in names]

if len(names_not_in_tree_taxa) > 0 :
    print("ERROR : the following names are in the {} file, but not in {} : {}".format(args.names, args.tree, names_not_in_tree_taxa))
    exit(1)

if len(tree_taxa_not_in_names) > 0 :
    
    if args.prune_extra_taxa == "True" :
        print("The following names are in {}, but not in the {} file : {}.\nThe tree will be pruned and the execution wiil continue as the '-prune-extra-taxa parameter' is set to 'True'".format(args.tree, args.names, tree_taxa_not_in_names))
        
    else :
        print("ERROR : The following names are in {}, but not in the {} file : {}.\nThe tree execution wiil be stopped as the '-prune-extra-taxa parameter' is not set to 'True'".format(args.tree, args.names, tree_taxa_not_in_names))
        exit(1)
        

        
new_tree = tree.extract_tree_with_taxa_labels(labels=names)

print("Final tree")
print(new_tree.as_ascii_plot())

new_tree.write(path=args.out, schema="newick")