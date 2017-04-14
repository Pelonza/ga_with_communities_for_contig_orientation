# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 08:50:14 2017

@author:Karl R. B. Schmitt


"""


#import ga_functions as ga_func #Contains self-defined functions for GA
#import load_data as ld #Contains the load_data function.
import os
import json
#from pathos.multiprocessing import ProcessingPool as Pool
import networkx as nx
import numpy as np
import random

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

from scoop import futures

from collections import OrderedDict
class OrderedNodeGraph(nx.Graph):
    node_dict_factory=OrderedDict

def load_data(input_filename):    
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
    
    return G, Initial_matepairs

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

#def trial(args, pool):
#    
#    G, Init_mps = load_data(args.ifilename)
#    
#    mypop, mystats, myhof, mylog = ga_func.gaopt_Uni(G, Init_mps, pool)
#    
#    with open(args.oorient, 'w+') as f:
#        json.dump(myhof[0], f)
#    with open(args.output_stats, 'w+') as f:
#        json.dump(mylog,f)
        
def is_valid_file(parser, arg):
    """
    Check if arg is a valid file that already exists on the file system.

    Parameters
    ----------
    parser : argparse object
    arg : str

    Returns
    -------
    arg
    """
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


def get_parser():
    """Get parser object for script xy.py."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--ifile",
                        dest="ifilename",
                        type=lambda x: is_valid_file(parser, x),
                        help="Read tallies from FILE",
                        metavar="FILE")
    parser.add_argument("-o", "--ofile",
                        dest="oorient",
                        help="Write best orientation to FILE",
                        metavar="FILE")
    parser.add_argument("-s", "--sfile",
                        dest="output_stats",
                        help="Write statistics to FILE",
                        metavar="FILE")    
    parser.add_argument("-n",
                        dest="n",
                        default=10,
                        type=int,
                        help="how many lines get printed")
    parser.add_argument("-q", "--quiet",
                        action="store_false",
                        dest="verbose",
                        default=True,
                        help="don't print status messages to stdout")
    return parser


if __name__ == "__main__":
    """
    Main that calls a single trial with currently default parameters.
    """
    
    args = get_parser().parse_args()
    
    #pool = Pool(15)
    
    G, Init_mps = load_data(args.ifilename)
    
    params=list([100,100,0.10,0.20,0.1,0.2])
    
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
    
    toolbox.register("evaluate", evaluate, G=G, Init_pairs=Init_mps)
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

    with open(args.oorient, 'w+') as f:
        json.dump(hof[0], f)
    with open(args.output_stats, 'w+') as f:
        json.dump(logbook,f)
    
    #trial(args, pool)


    
        
    
        
        
        
        
        









