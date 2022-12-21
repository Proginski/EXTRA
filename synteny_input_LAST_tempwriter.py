# -*- coding: utf-8 -*-
"""
This script performs micro-synteny for a reference genome (A) and its 
neighbour (B).
It takes as input : 
    - genomeA's gff
    - genomeB's gff
    - a file containning the pairs of orthologs between genomeA (first column),
    and genomeB (second column)
    - a file containning CDS from genomeA (first column), and their homologous
    intergenic regions in genomeB (columns 2 (sequence), 3 (start), and 4 (stop))
"""

### TODO : parallelization (not necessary 'til now)

from pybedtools.bedtool import BedTool
import pandas as pd
import re
from synteny import micro_synteny
import multiprocessing
import argparse
import time
import operator
from os import makedirs
from os.path import exists
from random import sample
import sys
import tempfile

# Output uniq values from a list while preserving their order
# adapted from https://www.peterbe.com/plog/uniqifiers-benchmark
def ordered_set(seq, exclude=[]):

    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x) or x in exclude)]


# With a gff line as bedtool interval object : get the ID
def getID(interval, attr, attrs_index = -1):
    
    attributes = re.split(pattern='=|;', string=interval[attrs_index])

    try : index = attributes.index(attr)+1
    except ValueError : print("{} is not found in {}".format(attr, interval.fields))
    
    try : id_value = attributes[index]
    except IndexError : print("There is no such index in attributes of {}".format(interval.fields)) 
    return(id_value)  


# For a given CDS-IGR pair, retrieve their flanking genes along with their
# strand and launch the micro_synteny function.
def launch_synteny(query):
    
    # temp_file
    tmp = tempfile.NamedTemporaryFile(delete=False, dir="temp")
    # print(tmp.name)
    
    # fd, path = tempfile.mkstemp(suffix=".txt",prefix="abc") #can use anything 
    # tmpo = os.fdopen(fd, 'w')
    
    # Redirect all prints to the output file
    sys.stdout = open(tmp.name, 'w')
    # sys.stdout = open(tmpo, 'w')


    print("Query {}/{}".format(queries.index(query)+1,len(queries)))
    

    query_cds,IGR_seq,IGR_start,IGR_stop = query
    
    print("\nCDS : {}\nIGR : {}:{}-{}".format(query_cds,IGR_seq,IGR_start,IGR_stop))
    
    # Convert 0-base BED start to 1-based GFF start
    IGR_start = IGR_start +1 
    
    
    # Based on an "closestX" object, retrieve "nb" flanking upstream or 
    # downstream genes of a sequence along with their strand.
    def get_flankers(closest_obj, side, nb, genome):

        # Define sign for comparison
        if side == "upstream" :
            comp = operator.lt
            query_bound = 3
            flanker_bound = 12
        elif side == "downstream":
            comp = operator.gt
            query_bound = 4
            flanker_bound = 13
        else : print("No valid side.")
        
        flankers_dic = {}
        
        includers_dic = {}
        
        if genome == "A" :
            
            # Restrict the iteration to precomputed indexes in "indsX" objects
            for interval in closest_obj[indsA[query_cds]["start"]:indsA[query_cds]["end"]] :
                
                # Stop after suffiscient number of flankers
                if len(flankers_dic) >= nb:
                    break
                
                # Ignore empty flankers
                elif interval.fields[-1] == "." : 
                    print("empty flanker")
                    continue
                
                else :
                    # For now only consider neighbour objects with distance 
                    # strictly greater or lower than zero
                    # if "ID={}".format(query_cds) in interval.fields[8] and comp(int(interval.fields[-1]), 0) :

                    #     flankers_dic[getID(interval, attr = ortho_attrA, attrs_index = -2)] = interval[-4]
                    # If the interval contains indeed the query,
                    if "{}={}".format(queries_attr, query_cds) in interval.fields[8] :
                        
                        # print("the interval contains the query")
                        
                        # If the flanker is at the given side,
                        if comp(int(interval.fields[-1]), 0) :
                            # Add it as key to dict with its strand as value.
                            flankers_dic[getID(interval, attr = ortho_attrA, attrs_index = -2)] = interval[-4]                        
                        # If the flankers overlapps the query :
                        elif interval.fields[-1] == 0 :
                            
                            query_start = int(interval.fields[3])
                            query_stop = int(interval.fields[4])
                            flanker_start = int(interval.fields[12])
                            flanker_stop = int(interval.fields[13])
                            
                            # If the flankers includes the query :
                            if flanker_start <= query_start and flanker_stop >= query_stop :
                                
                                print("inclusion : {}".format(interval.fields))
                                flankers_dic[getID(interval, attr = ortho_attrA, attrs_index = -2)] = interval[-4]
                                includers_dic[getID(interval, attr = ortho_attrA, attrs_index = -2)] = flanker_stop - flanker_start +1
                            
                            # If the flankers overlapps only the given side's end of the query,
                            elif comp(int(interval.fields[flanker_bound]),int(interval.fields[query_bound])) :
                                flankers_dic[getID(interval, attr = ortho_attrA, attrs_index = -2)] = interval[-4]
                        
                            else : continue
                    # else : print("the interval DOES NOT contain the query")

            
            if includers_dic != {} :
                # Remove the shortest includer gene as it is supposed to be the gene that contains the query CDS.
                shortest_includer = [includer for includer in includers_dic if includers_dic[includer] == min(includers_dic.values())][0]
                del flankers_dic[shortest_includer]
            
        elif genome == "B" :
            # "start-1" to get 0-based start 
            subset = indsB["{}:{}-{}".format(IGR_seq,IGR_start-1,IGR_stop)]
            for interval in closest_obj[subset["start"]:subset["end"]] :
                if len(flankers_dic) >= nb:
                    break
                elif interval.fields[-2] == "." :
                    continue
                else :
                    # "int(interval.fields[1])+1" to get 1-base start
                    if interval.fields[0] == IGR_seq and int(interval.fields[1])+1 == IGR_start and int(interval.fields[2]) == IGR_stop and comp(int(interval.fields[-1]), 0) :
        
                        flankers_dic[getID(interval, attr = ortho_attrB, attrs_index = -2)] = interval[-4]
                        
        else : print("No valid genome.")
                    
        return(flankers_dic)
    
    
    # The use of ordered_set function allows to remove redondant elements
    # including the query sequence if present. 
    query_upstream = get_flankers(closestA, "upstream", flankA, genome = "A")
    query_upstream_cds = query_upstream.keys()
    query_upstream_cds = ordered_set(query_upstream_cds, exclude=[query_cds])
    
    query_downstream = get_flankers(closestA, "downstream", flankA, genome = "A")
    query_downstream_cds = query_downstream.keys()
    query_downstream_cds = ordered_set(query_downstream_cds, exclude=[query_cds])
    
    IGR_upstream = get_flankers(closestB, "upstream", flankB, genome = "B")
    IGR_upstream_cds = IGR_upstream.keys()
    IGR_upstream_cds = ordered_set(IGR_upstream_cds)
    
    IGR_downstream = get_flankers(closestB, "downstream", flankB, genome = "B")
    IGR_downstream_cds = IGR_downstream.keys()
    IGR_downstream_cds = ordered_set(IGR_downstream_cds)
    
    
    micro_synteny(query_upstream_cds,
                  query_downstream_cds,
                  IGR_upstream_cds,
                  IGR_downstream_cds,
                  ortho_df,
                  strand_dic_A = {**query_upstream, **query_downstream},
                  strand_dic_B = {**IGR_upstream, **IGR_downstream},
                  ins_threshold = ins_threshold,
                  nb_flank = flankB)

    print("\n___________________________________________\n")
    
    
    return(tmp.name)
    
    
    
    
########## MAIN #############################################    

very_start_time = time.time()


print('Microsynteny is running...')
# Save a reference to the original standard output
original_stdout = sys.stdout 


# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("-A", "--gffA", help="gff for genomeA")
parser.add_argument("-B", "--gffB", help="gff for genomeB")
parser.add_argument("-X", "--ortho", help="a file containning the pairs of orthologs between genomeA (first column), and genomeB (second column)")
parser.add_argument("-L", "--list", help="a file containning four columns : 1 = CDS, 2 = IGR sequence 3 = IGR start 4 = IGR stop")
parser.add_argument("-N","--ncpus", type=int, help="number of cpus to use", default=1)
parser.add_argument("-C","--flankA", type=int, help="number of flanking CDS for genomeA", default=2)
parser.add_argument("-F","--flankB", type=int, help="number of flanking CDS for genomeB", default=4)
parser.add_argument("-T","--insertion", type=int, help="insertion threshold", default=0)
parser.add_argument("-O","--out", help="output_file")


args = parser.parse_args()


# Set the number of CDS to retrieve on each side of genome A or B
flankA = args.flankA
flankB = args.flankB
# Set the number of contiguous insertions to ignore before reporting one
ins_threshold = args.insertion


# gffA = "GENOMES/sorted_Scer_NCBI.gff"
# gffA = "GENOMES/sorted_Scer_SGD.gff"
gffA = BedTool(args.gffA)
gffB = BedTool(args.gffB)
# gffB="GENOMES/sorted_Spar.gff"
# gffB="GENOMES/sorted_Smik.gff"


# Redirect all prints to the output file
sys.stdout = open(args.out, 'w')


# Filter a gff object with respect to its third column
def filter_gff(gff, gff_types = ["gene", "pseudogene"], title = "the current gff"):
    
    filtered_gff = BedTool([inter for inter in gff if inter.fields[2] in gff_types])

    if gff_types == ["gene", "pseudogene"] :
        # Because in some cases gff containts "pseudogene" objects but no "gene" objects.
        gene_gff = BedTool([inter for inter in filtered_gff if inter.fields[2] in ["gene"]])
        if len(gene_gff) == 0 : 
            filtered_gff = BedTool([inter for inter in gff if inter.fields[2] in gff_types])
            print("As no 'gene' features where found in {}, filtering with 'CDS' instead".format(title))
            filtered_gff = filter_gff(gff, gff_types = ["CDS"])
        
    
    return(filtered_gff.sort())

# The filtered GFFs will be used to retrieve the flanking genes.
print("--- %s seconds for initialization---" % (time.time() - very_start_time))
start_time = time.time()
filtered_gffA = filter_gff(gffA, title = "gffA")
filtered_gffB = filter_gff(gffB, title = "gffB")
print("--- %s seconds to filter gffA and gffB---" % (time.time() - start_time))
start_time = time.time()   
# Build a pandas DataFrame with columns A and B for managing matching orthologs.
ortho_file = args.ortho
# ortho_file = "SYNTENY/Scer_NCBI_Spar_orthologues.csv"
ortho_df = pd.read_csv(ortho_file, sep="\t", names = ["A","B"], dtype=str)


# Build a pandas DataFrame with four columns : 1 = CDS, 2 = IGR sequence 3 = IGR start 4 = IGR stop
candidates = args.list
# candidates = "SYNTENY/Spar_synteny_list_with_IGR.txt"
candidates_df = pd.read_csv(candidates, sep="\t", names= ["query_cds","IGR_seq","IGR_start","IGR_stop"])
print("--- %s seconds for candidates---" % (time.time() - start_time))
start_time = time.time() 

# Determine which attribute to use as ID
def get_attr(gff, test_list, title, verbose = False, ) :
    
    # gff_read = BedTool(gff)
    gff_read = gff
    
    # Test on a sample of 10 sequences of "ortho_df"
    testers_ind = sample(range(0,len(test_list)),min(10,len(test_list)))    
    testers = [ test_list[ind] for ind in testers_ind ] 
          
    def find_attr(tester) :
        for interval in gff_read:
            attributes = re.split(pattern='=|;', string=interval[-1])
            match_attrs = []
            if tester in attributes : 
                match_attrs = [attributes[index-1] for index in range(0,len(attributes)) if attributes[index] == tester]
                if verbose : print(tester+str(match_attrs))
                break
        if match_attrs == [] : 
            print("{} is found in no attributes of the gff".format(tester))
            sys.stderr.write("{} is found in no attributes of the gff\n".format(tester))
            exit
        return(match_attrs)
    
    match_attrs_list = [ find_attr(tester) for tester in testers ]
    
    try : 
        attr_to_use = list(set.intersection(*map(set,match_attrs_list)))[0]
        # print("Using '{}' as attribute for '{}' based on the following names : {}".format(attr_to_use, gff,testers))
        # print("Using '{}' as attribute for {} based on the following names : {}".format(attr_to_use, title, testers))
        print("Using '{}' as attribute for {}".format(attr_to_use, title))

    except : 
        if verbose : return()
        print("No attribute in the given gff seems to contain all features for {} based on the following names : {}".format(title, testers))
        get_attr(gff, title=title, test_list=test_list, verbose = True)
        

    return(attr_to_use)

ortho_attrA = get_attr(filtered_gffA, title ="flanking sequences for genome A", test_list=ortho_df['A'])
ortho_attrB = get_attr(filtered_gffB, title ="flanking sequences for genome B", test_list=ortho_df['B'])
print("--- %s seconds to get ortho attributes---" % (time.time() - start_time))
start_time = time.time() 

# The "queries" object contains a list of every single row of the "candidates"
# DataFrame and will be iterated by the worker. 
queries = candidates_df.to_numpy().tolist()


# Build two objects : a GFF with only tested CDS of genome A
# and a BED with their matching IGR in genome B 
queries_cds = candidates_df["query_cds"].tolist()
print("--- %s seconds to get list from candidates---" % (time.time() - start_time))
start_time = time.time()
query_gffA = filter_gff(gffA, gff_types = ["CDS"])
print("--- %s seconds to filter gffA with CDS---" % (time.time() - start_time))
start_time = time.time()
queries_attr = get_attr(query_gffA, title ="queries", test_list=queries_cds)
print("--- %s seconds to get queries attributes---" % (time.time() - start_time))
start_time = time.time()
# queriesgff = BedTool([ line for line in BedTool(query_gffA).sort() if getID(line, attr=queries_attr) in queries_cds ])
queriesgff = BedTool([ line for line in query_gffA if getID(line, attr=queries_attr) in queries_cds ])
print("--- %s seconds to build querygff---" % (time.time() - start_time))
start_time = time.time() 
igr_bed = BedTool.from_dataframe(candidates_df[["IGR_seq","IGR_start","IGR_stop"]]).sort()
print("--- %s seconds to build igr_bed---" % (time.time() - start_time))
start_time = time.time()

# Retrieve the k closest objects (genes) to :
# the tested CDS in genome A (searching in genome A),
closestA = BedTool(queriesgff).closest(filtered_gffA, k= 1000, D ="ref") 
# their corresponding IGR in genome B (searchig in genome B)
closestB = BedTool(igr_bed).closest(filtered_gffB, k= 1000, D ="ref")
print("--- %s seconds for closest objects---" % (time.time() - start_time))
start_time = time.time()

# For both "closestX" objects : get the start and stop index of each single query
# Will be used to efficiently get the flanking objects for a given query.
indsA = {}
for i, interval in enumerate(closestA) :
    ID = getID(interval, attr=queries_attr, attrs_index= 8)
    if ID not in indsA.keys() :
        indsA[str(ID)] = {}
        indsA[str(ID)]["start"] = i
    indsA[str(ID)]["end"] = i+1

indsB = {}
for i, interval in enumerate(closestB) :
    ID = "{}:{}-{}".format(interval.fields[0],interval.fields[1],interval.fields[2])
    if ID not in indsB.keys() :
        indsB[str(ID)] = {}
        indsB[str(ID)]["start"] = i
    indsB[str(ID)]["end"] = i+1
print("--- %s seconds for indexation of closest objects---" % (time.time() - start_time))


print('\n{} queries to process.\n'.format(len(queries)))


# Check whether the temp dir exists or not. It will be used to store temp files
temp_dir = "temp"
if not exists(temp_dir):
    print("creating a temp dir")
    makedirs(temp_dir)

# Number of cpus to use for parallelization            
ncpus = args.ncpus
if ncpus == 0 : ncpus = 1
# Parallel work
pool = multiprocessing.Pool(ncpus)
tmp_names = pool.map(launch_synteny, queries)
pool.close()
pool.join()

# print("temporary files to append :")
# print(tmp_names)

# Append all temp files to the output file
# with open(args.out, 'a') as outfile:
#     for tmp_name in tmp_names:
#         with open(tmp_name) as infile:
#             outfile.write(infile.read())
for tmp_name in tmp_names:
      with open(tmp_name, 'r') as f:           
          print(f.read())

# # Try to solve the non-run pb
# from time import sleep
# sleep(10)

print("JOB FULLY COMPLETED.")
print("--- %s seconds in total---" % (time.time() - very_start_time))

# Reset the standard output to its original value
sys.stdout = original_stdout 
