#This code will first test/experiment with finding the fitness of a given 
# orientation based on a tally file.

from deap import base
from deap import creator
from deap import tools
import random
import numpy as np
import networkx as nx
import pandas as pd
import os

os.chdir('C:\\Users\\kschmit1\\Documents\\GitHub\\ga_with_communities_for_contig_orientation\\data\\test')


input_filename='basic.tally'

## NetworkX Initializations and data reading.

#Have networkX use an ordereddict to track the nodes, to guarantee a consistent
# mapping of orientations from individuals to contigs.
from collections import OrderedDict
class OrderedNodeGraph(nx.Graph):
    node_dict_factory=OrderedDict

#Read in tally file, assuming a 6-column format, no duplicated edges
G=nx.read_edgelist(input_filename,nodetype=int, data=(('w1',int),('w2',int),('w3',int),('w4',int)), create_using=OrderedNodeGraph())

#Set all edges to 'flipable' and give each node a fixed index for referencing.
nx.set_node_attributes(G, 'flippable', True)
i=0
for n in nx.nodes_iter(G):
    G.node[n]['idx']=i
    i=i+1

#Compute basic data about the tally file/input data.
Initial_matepairs=list([0,0,0]) #(Total, Good, Bad) Total_matepairs=0
for (u,v,d) in G.edges(data=True):
    Initial_matepairs[0]=Initial_matepairs[0]+d['w1']+d['w2']+d['w3']+d['w4']
    Initial_matepairs[1]=Initial_matepairs[1]+d['w1']
    Initial_matepairs[2]=Initial_matepairs[2]+d['w2']

#Set some GA parameters:    
IND_SIZE = nx.number_of_nodes(G)  #Size of each individual in the population
MUT_IDPB = 0.1 #Independent probability of an attribute (orientation) mutating
POP_SIZE = 100  #Size of the overall population
NGEN     = 100  #Number of generations to evolve for
MUT_PB   = 0.10 #Probability that an offspring will mutate
CX_PB    = 0.20 #Probability that two children will perform a crossover 

## End data reading with networkX

#Define an fitness-evaluation function.
def evaluate(individual, G=G, Init_pairs=Initial_matepairs):
    #Inputs are one individual (with array of flip/no-flip orientations), and 
    # a Graph structure (with nodes/edges from tally file). It returns the 
    # fitness of the individual.       
    Orient_pairs=list(Init_pairs)
    print(Orient_pairs)
    for (u,v,d) in G.edges(data=True):
        if individual[G.node[u]['idx']]!=individual[G.node[v]['idx']]:
            Orient_pairs[1]=Orient_pairs[1]-d['w1']+d['w2']
            Orient_pairs[2]=Orient_pairs[2]-d['w2']+d['w1']
    
    print(Orient_pairs)
    return (Orient_pairs[1]/Orient_pairs[0], )

#Define objects for/from DEAP's toolbox
#Currently modeled largely off of their one_max example.
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox=base.Toolbox()
toolbox.register("attr_bool", random.randint, 0, 1)
toolbox.register("individual",tools.initRepeat,creator.Individual,toolbox.attr_bool, n=IND_SIZE)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate",evaluate, G=G, Init_pairs=Initial_matepairs)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=MUT_IDPB)
toolbox.register("select", tools.selTournament, tournsize=3)

stats= tools.Statistics()
stats.register("mean", np.mean)
stats.register("max", max)
stats.register("min", min)
stats.register("std", np.std)


## Main function that actually uses all the toolbox...
def trial(POP=POP_SIZE, CX=CX_PB, MUT=MUT_PB, GENS=NGEN):
    

    pop= toolbox.population(n=POP)
    fitnesses=list(map(toolbox.evaluate,pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values=fit
        
    #Begin Evolution
    for g in range(NGEN):
        print("-- Generation %i --" % g)
        
        evo_pop_offspring = toolbox.select(pop, len(pop))
        offspring = list((map(toolbox.clone,evo_pop_offspring)))
   
        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < CX:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:
            if random.random() < MUT:
                toolbox.mutate(mutant)
                del mutant.fitness.values   
                
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        
        pop[:] = offspring
        
        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]
        stats.compile(fits)
        print(stats.compile(fits))
        
    return pop, stats

def main():
    (outpop, stats)=trial()
    print(outpop)
    print(stats)

if __name__ == "__main__":
    main()
    
        
        
        
        
        









