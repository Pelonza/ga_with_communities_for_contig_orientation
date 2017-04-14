# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 08:57:38 2017

@author: Karl R.B. Schmitt

This script contains definitions of all the functions used within the DEAP 
package for GA's. This includes unique implementations of crossovers, evaluate
etc.

This module contains a basic GA maximization optimization using the DEAP package.
It uses a uniform crossover and 'flipbit' mutation.

Running as an independent function requires the "load_data.py" script to be in
the same folder as this code.

Tests against 'basic.tally' in the ..\\..\\data\\test folder. 
"""

#Can try commenting this out and see if function still works!
from deap import base
from deap import creator
from deap import tools
from deap import algorithms

#multiprocessing
#import multiprocessing
#from pathos.multiprocessing import ProcessingPool as Pool
from scoop import futures

import random
import numpy as np
import networkx as nx  #Defines G.edges and G.node 

#Define an fitness-evaluation function.
def evaluate(individual, G, Init_pairs):
    #Inputs are one individual (with array of flip/no-flip orientations), and 
    # a Graph structure (with nodes/edges from tally file). It returns the 
    # fitness of the individual.       
    Orient_pairs=list(Init_pairs)
    for (u,v,d) in G.edges(data=True):
        if individual[G.node[u]['idx']]!=individual[G.node[v]['idx']]:
            Orient_pairs[1]=Orient_pairs[1]-d['w1']+d['w2']
            Orient_pairs[2]=Orient_pairs[2]-d['w2']+d['w1']
    
    return (Orient_pairs[1]/Orient_pairs[0], )

def gaopt_Uni(G, Initial_matepairs, pool, params=list([100,100,0.10,0.20,0.1,0.2])):
    #Note this is included for readability/reference.
    #Set some GA parameters:
    IND_SIZE = nx.number_of_nodes(G)    #Size of each individual in the population
    POP_SIZE = params[0]                #Size of the overall population
    NGEN     = params[1]                #Number of generations to evolve for
    MUT_PB   = params[2]                #Probability that an offspring will mutate
    CX_PB    = params[3]                #Probability that two children will perform a crossover 
    MUT_IDPB = params[4]                #Independent probability of an attribute (orientation) mutating
    CX_IDPB  = params[5]                #Independent probability of an attribute (orientation) crossingover
    
    #Define objects for/from DEAP's toolbox
    #Currently modeled largely off of their one_max example.
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)
    
    toolbox=base.Toolbox()
    toolbox.register("attr_bool", random.randint, 0, 1)
    toolbox.register("individual",tools.initRepeat,creator.Individual,toolbox.attr_bool, n=IND_SIZE)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    
    toolbox.register("evaluate", evaluate, G=G, Init_pairs=Initial_matepairs)
    toolbox.register("mate", tools.cxUniform, indpb=CX_IDPB)
    toolbox.register("mutate", tools.mutFlipBit, indpb=MUT_IDPB)
    toolbox.register("select", tools.selTournament, tournsize=3)
    
    #toolbox.register("map", pool.map)
    toolbox.register("map", futures.map)
    
    hof= tools.HallOfFame(1)
    
    stats= tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("mean", np.mean)
    stats.register("max", max)
    stats.register("min", min)
    stats.register("std", np.std)
    
    logbook=tools.Logbook()
    
    random.seed(1)
    pop= toolbox.population(n=POP_SIZE)
    
    pop, logbook = algorithms.eaSimple(pop,toolbox,cxpb=CX_PB, mutpb=MUT_PB, ngen=NGEN, stats=stats, halloffame=hof, verbose=False)
    
    return pop, stats, hof, logbook


if __name__=='__main__':
    """
    Test Evaluate & combined optimization.
    
    Depends on module: load_data being in same folder as ga_functions.py
    """
    
    import os
    import load_data as ld
    
    os.chdir("..\\..\\data\\test")
    
    G, Init = ld.load_data('basic.tally')
    
    #Check to make sure the test file loaded looks like we expect.
    assert nx.number_of_nodes(G)==5
    assert nx.number_of_edges(G)==6
    assert Init[0]==32
    assert Init[1]==24
    assert Init[2]==8
    
    #Basic orientations possible in an individual
    no_flip=list([0,0,0,0,0])
    one_flip=list([0,0,0,0,1])
    mny_flip=list([1,1,1,1,0])
    
    a = evaluate(no_flip, G, Init)
    
    assert a[0]==0.75
    #assert a[1]==None
    
    a = evaluate(one_flip, G, Init)
    
    assert a[0]==0.84375
    #assert b==None
    
    a = evaluate(mny_flip, G, Init)
    
    assert a[0]==0.84375
    #assert b==None
    
    #Check for correct clearing of repeated evaluates
    a = evaluate(one_flip, G, Init)
    
    assert a[0]==0.84375
    #assert b==None
    
    a = evaluate(mny_flip, G, Init)
    
    assert a[0]==0.84375
    #assert b==None
    
    #Confirm our initial data hasn't gotten modified
    assert nx.number_of_nodes(G)==5
    assert nx.number_of_edges(G)==6
    assert Init[0]==32
    assert Init[1]==24
    assert Init[2]==8
    
    print("Testing of evaluate in ga_functions has succeeded")
    
    #----------------------
    print("Start testing the full optimization scheme")
    
    G, Initial_matepairs = ld.load_data("basic.tally")
    outpop, outstats, outhof, outlog = gaopt_Uni(G, Initial_matepairs, list([100,100,0.10,0.20,0.1,0.2]))
    
    try:
        assert outhof[0]==list([0,0,0,0,1])
    except:
        assert outhof[0]==list([1,1,1,1,0])
        
    print("Passed optimization of basic.tally")