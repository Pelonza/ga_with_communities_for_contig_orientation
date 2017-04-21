# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 08:50:14 2017

@author:Karl R. B. Schmitt


"""
from deap import base
from deap import creator
from deap import tools
from deap import algorithms

import random
import os
import multiprocessing

import json
#import dill as pickle
import pickle
import igraph as ig

import numpy as np
import networkx as nx  # Defines G.edges and G.node

from scoop import futures #imports the scoop distributed computing package.


# %%
# These have to be global for multiprocessing and scoop (?).
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

# Register all the data we want to track, in a global box for multiprocessing
stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("mean", np.mean)
stats.register("max", max)
stats.register("min", min)
stats.register("std", np.std)

hof_local = tools.HallOfFame(1)
hof_trials = tools.HallOfFame(1)

# Register elements of the toolbox that apply to both full and community GA.
# This does not necessarily need to be global (I think?)
toolbox = base.Toolbox()
toolbox.register("attr_bool", random.randint, 0, 1)
toolbox.register("select", tools.selTournament, tournsize=3)

# Making a global i-graphs to hold data for all workers to reduce overhead
# Only using one to preserve reusuability of trials function.
#global G_global
#G_global = ig.Graph()


# %%
def load_data(input_filename):
    # Read in tally file, assuming a 6-column format, no duplicated edges
    G = nx.read_edgelist(input_filename, nodetype=int,
                         data=(('w1', int), ('w2', int), ('w3', int),
                               ('w4', int)))

    # Set all edges to 'flipable' -- To be modified later if non-flippable.
    nx.set_node_attributes(G, 'flippable', True)
    
    # Write to GML then load as an igraph graph for faster computations on
    # other elementals.
    nx.write_gml(G, "tally.gml")
    ig_G = ig.Graph()
    ig_G = ig_G.Read_GML("tally.gml")

    
    for edge in range(ig_G.ecount()):
        # Define an additional edge attribute of total pairs
        ig_G.es[edge]['mates'] = ig_G.es[edge]['w1'] + ig_G.es[edge]['w2'] 
        # If col3/4 included    # + ig_G.es[edge]['w3'] + ig_G.es[edge]['w4'])

    Initial_matepairs = compute_initial_pairs(ig_G)
        
    return ig_G, Initial_matepairs


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
                        dest="ntrials",
                        default=50,
                        type=int,
                        help="how many trials to do per parameter set")
    parser.add_argument("-q", "--quiet",
                        action="store_false",
                        dest="verbose",
                        default=True,
                        help="don't print status messages to stdout")
    return parser


# %%
# Define an fitness-evaluation function to be used repeatedly in DEAP's toolbox
def evaluate(individual, G, Init_pairs):
    # Inputs are one individual (with array of flip/no-flip orientations), and
    # a Graph structure (with nodes/edges from tally file). It returns the
    # fitness of the individual.
    Orient_pairs = list(Init_pairs)
    for edge in range(G.ecount()):
        if individual[G.es[edge].tuple[0]] != individual[G.es[edge].tuple[1]]:
            Orient_pairs[1] = Orient_pairs[1]-G.es[edge]['w1'] + G.es[edge]['w2']
            Orient_pairs[2] = Orient_pairs[2]-G.es[edge]['w2'] + G.es[edge]['w1']

    return (Orient_pairs[1]/Orient_pairs[0], )


def Run_GA(G, params):
    # Sets up and runs the GA using DEAP on the input graph.
    
    # Setup GA
    # Note this is included for readability/reference.
    # Set some GA parameters:
    IND_SIZE = G.vcount()   # Size of individual in population
    POP_SIZE = params[0]    # Size of the overall population
    NGEN     = params[1]    # Number of generations to evolve for
    MUT_PB   = params[2]    # Probability that an offspring will mutate
    CX_PB    = params[3]    # Probability that two children will crossover 
    MUT_IDPB = params[4]    # Independent probability of an attribute mutating
    #CX_IDPB  = params[5]    # Independent probability of an attribute crossover
    
    # Define objects for/from DEAP's toolbox
    Init_mps = compute_initial_pairs(G)
    
    try:
        assert Init_mps[0] != 0
    except:
        print("No matepairs in graph for Run_GA")

    toolbox.register("individual", tools.initRepeat, creator.Individual,
                     toolbox.attr_bool, n=IND_SIZE)
    toolbox.register("population", tools.initRepeat, list,
                     toolbox.individual)
    toolbox.register("evaluate", evaluate, G=G,
                     Init_pairs=Init_mps)
    #toolbox.register("mate", tools.cxUniform, indpb=CX_IDPB)
    toolbox.register("mate", tools.cxOnePoint)
    toolbox.register("mutate", tools.mutFlipBit, indpb=MUT_IDPB)
    
    hof_local.clear()
    
    pop = toolbox.population(n=POP_SIZE)
    pop, tlogbook = algorithms.eaSimple(pop, toolbox, cxpb=CX_PB, mutpb=MUT_PB,
                                        ngen=NGEN, stats=stats,
                                        halloffame=hof_local,
                                        verbose=False)
        
    best_ort = hof_local[0]
    
    # Clean-up the toolbox from this call to run a GA
    toolbox.unregister("individual")
    toolbox.unregister("population")
    toolbox.unregister("evaluate")
    toolbox.unregister("mate")
    toolbox.unregister("mutate")
    
    return pop, best_ort, tlogbook


def merge_ort(orient1, orient2, outort):
    for i in range(len(orient1)):
        outort[i]=orient1[i]^orient2[i]        
    return


def compute_initial_pairs(G):
    # Compute basic data about the tally file/input data.
    # Uses an igraph -style graph.
    # Assumes edges in G have "mates", "w1" and "w2" properties
    Initial_matepairs = list([0, 0, 0])  # (Total, Good, Bad) Total_matepairs=0
    for edge in range(G.ecount()):
        Initial_matepairs[0] = Initial_matepairs[0] + G.es[edge]['mates']        
        Initial_matepairs[1] = Initial_matepairs[1] + G.es[edge]['w1']
        Initial_matepairs[2] = Initial_matepairs[2] + G.es[edge]['w2']
    
    return Initial_matepairs


def update_graph(G, orient):
    for edge in range(G.ecount()):
        if orient[G.es[edge].tuple[0]] != orient[G.es[edge].tuple[1]]:
            #Swap 'w1' and 'w2' (also 'w3' and 'w4' to reduce fixing later)
            tmp = G.es[edge]['w1']
            G.es[edge]['w1'] = G.es[edge]['w2']
            G.es[edge]['w2'] = tmp
            tmp = G.es[edge]['w3']
            G.es[edge]['w3'] = G.es[edge]['w4']
            G.es[edge]['w4'] = tmp
    
    return


def trials(allparam):
    # Attempt to move trial loop to a function for scoop mapping.
    
    param_logbook=tools.Logbook() #Create a logbook to hold the trials.
    # Loop here
    for i in range(allparam[0]):
        
        # ----------------
        # This chunk runs a GA then records it as a trial.
        pop, tbestort, tlogbook = Run_GA(allparam[1], allparam[2])
        param_logbook.record(trial=i, tmax=tlogbook.select('max'),
                            tbort=tbestort, tgen=tlogbook.select('gen'),
                            tmin=tlogbook.select('min'),
                            tstd=tlogbook.select('std'),
                            tmean=tlogbook.select('mean'),
                            tparam=allparam[2])
        hof_trials.update(pop)
    
    return param_logbook


# %%
if __name__ == "__main__":
    """
    Main that calls a single parameter trial, with multiple runs (args=n)
    with currently default parameters.

    Some naming conventions:
        *_full  -- Variables related to the full, original tally file.
        *_comm  -- Variables related to the reduce, community based data
    """
    
    #global G_global  # Allow G_global to be modified in this script.

    args = get_parser().parse_args()
    random.seed(1)  # For testing/reproducibility.
    
    #pool = multiprocessing.Pool()
    #toolbox.register("map", pool.map)
    
    # Located near top to ease changing them. May be overwritten in testing.
    #params_full = list([200, 10, 0.10, 0.20, 0.1, 0.2])
    params_full = list([200, 10, 0.10, 0.20, 0.1]) # Not Uniform CX
    #params_comm = list([200, 10, 0.10, 0.20, 0.1, 0.2])
    
    #Replicated from Run_GA for reference.
    #POP_SIZE = params[0]    # Size of the overall population
    #NGEN     = params[1]    # Number of generations to evolve for
    #MUT_PB   = params[2]    # Probability that an offspring will mutate
    #CX_PB    = params[3]    # Probability that two children will crossover 
    #MUT_IDPB = params[4]    # Independent probability of an attribute mutating
    #CX_IDPB  = params[5]    # Independent probability of an attribute crossover
    
    # %%
    # ===========
    # Load data and (if required) determine community structure and remake
    # tally into community-tally.
    G_full, Init_mps_full = load_data(args.ifilename)
    
    toolbox.register("orient", tools.initRepeat, creator.Individual,
                 int, n=G_full.vcount())
    
    base_orient = toolbox.orient()  # For storing the original orientation(s)
    tmp_orient = toolbox.orient()   # For use in merging/swapping orientations.
    cls_orient = toolbox.orient()   # To unwrap clustered orientation into.

    try:
        assert Init_mps_full[0] != 0
    except:
        print("No matepairs in graph")

    # Generate a reduced by community graph.
    # Cluster membership of a given node can be dereferenced by:
    # *_clusters.membership[node] thus the orientation of a node from cluster
    # is: clus_ort[*_clusters.membership[node]] 
    G_full_dendrogram = G_full.community_fastgreedy(weights="mates")
    G_full_clusters = G_full_dendrogram.as_clustering()
    G_comm = G_full_clusters.cluster_graph(combine_edges=sum)

    
    # %%
    
    logbook=tools.Logbook()
    full_logbook=tools.Logbook()  # Logbook for running just the base GA
    
    #Parameter Sweep
    #Base Parameters:
    #param = list([50, 2000, 0.1, 0.1, 0.05, 0.1])
    param = list([50, 1000, 0.1, 0.1, 0.05])
    
    # Set which graph we are testing on.
    G_global = G_full
    
    # ============
    # When using "update_graph" MUST MAKE COPY OF GRAPH FOR EACH DISTRIBUTED
    # COMPUT NODE! -- Python auto-passes by REFERENCE, so otherwise solutions
    # will get all super jumbled up!
    # ============
    
    
    # Pre-declare mapdata as a list of lists.
    mapdata = list(list())
    for mut_pb in range(0,10,5):
        #Set sweep parameter really by 0.05
        param[2] = mut_pb/100

        mapdata.append(list([4, G_full, list(param)]))
        
    full_logbook=list(futures.map(trials, mapdata))

        
    print("Finished trials of parameter") #: ", param[2])
        
    with open(args.output_stats, 'w') as f:
        json.dump(full_logbook,f)





# Mouse line:
    #-i ..\\..\\data\\raw\\mouse_tallies -o ..\\..\\data\\interim\\mouse.orient -s ..\\..\\data\\interim\\mouse.stat
    
# Test cluster lines:
# -i ..\\..\\data\\test\\basic2_cluster2.tally -s ..\\..\\data\\interim\\cluster.tally        
        
        
        
        









