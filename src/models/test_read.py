#This code will first test/experiment with finding the fitness of a given 
# orientation based on a tally file.

import deap
import random
import numpy as np
import networkx as nx
import pandas as pd


#This should be replaced by determining the number of unique contigs in a 
# tally file.
filename='tallies'

fh=open(filename,'rb')

#Read in tally file, assuming a 6-column format, no duplicated edges
G=nx.read_edgelist(fh,nodetype=int, data=(('w1',int),('w2',int),('w3',int),('w4',int)))

#Set all edges to 'flipable'
nx.set_node_attributes(G, 'flippable', True)


#G.add_weighted_edges_from(list), list has weights in list.

