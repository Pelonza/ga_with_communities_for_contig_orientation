#This code will first test/experiment with finding the fitness of a given 
# orientation based on a tally file.

import deap
import random
import numpy as np
import networkx as nx
import pandas as pd


input_filename='basic.tally'

#Have networkX use an ordereddict to track the nodes, to guarantee a consistent
# mapping of orientations from individuals to contigs.
from collections import OrderedDict
class OrderedNodeGraph(nx.Graph):
    node_dict_factory=OrderedDict

#Read in tally file, assuming a 6-column format, no duplicated edges
G=nx.read_edgelist(input_filename,nodetype=int, data=(('w1',int),('w2',int),('w3',int),('w4',int)), create_using=OrderedNodeGraph())

#Set all edges to 'flipable'
nx.set_node_attributes(G, 'flippable', True)
i=0
for n in nx.nodes_iter(G):
    G.node[n]['idx']=i
    i=i+1

Initial_matepairs=(0,0,0) #(Total, Good, Bad) Total_matepairs=0
#Initial_Good=0 #Only adding first column
#Initial_Bad=0 #Only adding second column
for (u,v,d) in G.edges(data=True):
    Initial_matepairs[0]=Initial_matepairs[0]+d['w1']+d['w2']+d['w3']+d['w4']
    Initial_matepairs[1]=Initial_matepairs[1]+d['w1']
    Initial_matepairs[2]=Initial_matepairs[2]+d['w2']
    
IND_SIZE=number_of_nodes(G)

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", np.ndarray, fitness=creator.FitnessMax)

toolbox=base.Toolbox()
toolbox.register("individual",tools.initRepeat,creator.Individual,np.random.randint(2), n=IND_SIZE)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.population(n=100)

def evaluate(individual, G, Total_matepairs):
    #Inputs are one individual (with array of flip/no-flip orientations), and 
    # a Graph structure (with nodes/edges from tally file). It returns the 
    # fitness of the individual.
    satisfied
    for (u,v,d) in G.edges(data=True):
        






