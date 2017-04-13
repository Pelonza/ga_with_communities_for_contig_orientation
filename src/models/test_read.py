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
IND_SIZE= 10  #Hardcodes an individual of size 10




creator.create("FitnessMax", base.Fitness, weights=(1.0,))

#G.add_weighted_edges_from(list), list has weights in list.

