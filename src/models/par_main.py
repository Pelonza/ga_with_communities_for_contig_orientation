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
import dill as pickle
import igraph as ig

import numpy as np
import networkx as nx  # Defines G.edges and G.node

# This may not be needed actually after moving to igraph.
#from collections import OrderedDict  # For modifying networkX graph-dic


# Define an fitness-evaluation function.
def evaluate(individual, G, Init_pairs):
    # Inputs are one individual (with array of flip/no-flip orientations), and
    # a Graph structure (with nodes/edges from tally file). It returns the
    # fitness of the individual.
    Orient_pairs = list(Init_pairs)
    for edge in range(G.ecount()):
        if individual[G.es[edge].tuple[0] != G.es[edge].tuple[1]]:
            Orient_pairs[1] = Orient_pairs[1]-G.es[edge]['w1'] + G.es[edge]['w2']
            Orient_pairs[2] = Orient_pairs[2]-G.es[edge]['w2'] + G.es[edge]['w1']

    return (Orient_pairs[1]/Orient_pairs[0], )


# Needs to be outside the function definition so it can be picklable.

# class OrderedNodeGraph(nx.Graph):
#    node_dict_factory = OrderedDict


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

    # Compute basic data about the tally file/input data.
    Initial_matepairs = list([0, 0, 0])  # (Total, Good, Bad) Total_matepairs=0
    for edge in range(ig_G.ecount()):
        # Define an additional edge attribute of total pairs
        ig_G.es[edge]['mates'] = ig_G.es[edge]['w1'] + ig_G.es[edge]['w2'] 
        # If col3/4 included    # + ig_G.es[edge]['w3'] + ig_G.es[edge]['w4'])

        Initial_matepairs[0] = Initial_matepairs[0] + ig_G.es[edge]['mates']        
        Initial_matepairs[1] = Initial_matepairs[1] + ig_G.es[edge]['w1']
        Initial_matepairs[2] = Initial_matepairs[2] + ig_G.es[edge]['w2']
        

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
                        default=2,
                        type=int,
                        help="how many trials to do per parameter set")
    parser.add_argument("-q", "--quiet",
                        action="store_false",
                        dest="verbose",
                        default=True,
                        help="don't print status messages to stdout")
    return parser

# %%
# These have to be global for multiprocessing.
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

# Register all the data we want to track, in a global box for multiprocessing
stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("mean", np.mean)
stats.register("max", max)
stats.register("min", min)
stats.register("std", np.std)

# Register elements of the toolbox that apply to both full and community GA.
# This does not necessarily need to be global (I think?)
toolbox = base.Toolbox()
toolbox.register("attr_bool", random.randint, 0, 1)
toolbox.register("select", tools.selTournament, tournsize=3)

# %%

if __name__ == "__main__":
    """
    Main that calls a single parameter trial, with multiple runs (args=n)
    with currently default parameters.

    Some naming conventions:
        *_full  -- Variables related to the full, original tally file.
        *_comm  -- Variables related to the reduce, community based data
    """

    args = get_parser().parse_args()
    random.seed(1)  # For testing/reproducibility.
    
    # Located near top to ease changing them. May be overwritten in testing.
    params_full = list([200, 10, 0.10, 0.20, 0.1, 0.2])
    params_comm = list([200, 10, 0.10, 0.20, 0.1, 0.2])
    # %%
    # ===========
    # Load data and (if required) determine community structure and remake
    # tally into community-tally.
    G_full, Init_mps_full = load_data(args.ifilename)
    
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

    # Compute basic data about the tally file/input data.
    Init_mps_comm = list([0, 0, 0])  # (Total, Good, Bad) Total_matepairs=0
    for edge in range(G_comm.ecount()):
        Init_mps_comm[0] = Init_mps_comm[0] + G_comm.es[edge]['mates']        
        Init_mps_comm[1] = Init_mps_comm[1] + G_comm.es[edge]['w1']
        Init_mps_comm[2] = Init_mps_comm[2] + G_comm.es[edge]['w2']
    
    try:
        assert Init_mps_comm[0] != 0
    except:
        print("No matepairs (remained) in collapsed graph")

    # %%
    # =======================
    # Build toolbox for GA on the full original tally file given as input.
    # =======================
    
    # Note this is included for readability/reference.
    # Set some GA parameters:
    IND_SIZE = G_full.vcount()   # Size of individual in population
    POP_SIZE = params_full[0]    # Size of the overall population
    NGEN     = params_full[1]    # Number of generations to evolve for
    MUT_PB   = params_full[2]    # Probability that an offspring will mutate
    CX_PB    = params_full[3]    # Probability that two children will crossover 
    MUT_IDPB = params_full[4]    # Independent probability of an attribute mutating
    CX_IDPB  = params_full[5]    # Independent probability of an attribute crossover
    
    # Define objects for/from DEAP's toolbox
    # Currently modeled largely off of their one_max example.

    toolbox.register("individual_full", tools.initRepeat, creator.Individual,
                     toolbox.attr_bool, n=IND_SIZE)
    toolbox.register("population_full", tools.initRepeat, list,
                     toolbox.individual_full)
    toolbox.register("evaluate_full", evaluate, G=G_full,
                     Init_pairs=Init_mps_full)

    # Above could be done outside of "__main__" ??
    # ====================================

    # These have to be modified/looped over to modify the
    #  INDIVIDUAL crossover or mutation rate of orientations.
    toolbox.register("mate", tools.cxUniform, indpb=CX_IDPB)
    toolbox.register("mutate", tools.mutFlipBit, indpb=MUT_IDPB)

    pops_full = [toolbox.population_full(n=POP_SIZE) for i in range(args.ntrials)]
    # ===============================================
    # %%
    # =======================
    # Build toolbox for GA on the community reduced graph..
    # =======================
    
    # Note this is included for readability/reference.
    # Set some GA parameters:
    IND_SIZE = G_comm.vcount()   # Size of individual in population
    POP_SIZE = params_comm[0]    # Size of the overall population
    NGEN     = params_comm[1]    # Number of generations to evolve for
    MUT_PB   = params_comm[2]    # Probability that an offspring will mutate
    CX_PB    = params_comm[3]    # Probability that two children will crossover 
    MUT_IDPB = params_comm[4]    # Independent probability of an attribute mutating
    CX_IDPB  = params_comm[5]    # Independent probability of an attribute crossover
    
    # Define objects for/from DEAP's toolbox
    # Currently modeled largely off of their one_max example.

    toolbox.register("individual_comm", tools.initRepeat, creator.Individual,
                     toolbox.attr_bool, n=IND_SIZE)
    toolbox.register("population_comm", tools.initRepeat, list,
                     toolbox.individual_comm)
    toolbox.register("evaluate_comm", evaluate, G=G_comm,
                     Init_pairs=Init_mps_comm)

    # Above could be done outside of "__main__" ??
    # ====================================

    # These have to be modified/looped over to modify the
    #  INDIVIDUAL crossover or mutation rate of orientations.
    toolbox.register("mate", tools.cxUniform, indpb=CX_IDPB)
    toolbox.register("mutate", tools.mutFlipBit, indpb=MUT_IDPB)

    pops_comm = [toolbox.population_comm(n=POP_SIZE) for i in range(args.ntrials)]
    # ===============================================
    # %%

    # Below copies DEAP's multiprocessing one-max example inside __main__
    # Includes a loop over multiple trials.
    #
    # Tried to use "scoop" for distributed computing (I think that'd be better)
    #  but it didn't seem to work right. Would need to determine what actually
    #  can be pickled/dill'd or not from DEAP, networkX and igraph.
    
    pool = multiprocessing.Pool()
    toolbox.register("map", pool.map)

    # I think these need to be unregistered to get the multiprocessing to work.
    toolbox.unregister("attr_bool")
    toolbox.unregister("individual_full")
    toolbox.unregister("population_full")
    toolbox.unregister("individual_comm")
    toolbox.unregister("population_comm")

    hof = tools.HallOfFame(1)

    logbook = tools.Logbook()  # Logbook for all trials.

    # Loop here to change crossover or mutation rates (CX_PB and MUT_PB).

    # Have to explicitly register an 'evaluate' function to use simple alg.
    toolbox.register("evaluate", toolbox.evaluate_comm)
    
    for jj in range(args.ntrials):
        pop, tlogbook = algorithms.eaSimple(pops_comm[jj], toolbox, cxpb=CX_PB,
                                            mutpb=MUT_PB, ngen=NGEN,
                                            stats=stats, halloffame=hof,
                                            verbose=False)
        logbook.record(trial=jj, tmax=tlogbook.select('max'),
                       tbestort=hof[0], tgen=tlogbook.select('gen'),
                       tmin=tlogbook.select('min'),
                       tstd=tlogbook.select('std'),
                       tmean=tlogbook.select('mean'))

    # End run trials of GA on full tally
    # ====================================
    
    # %%



    print(logbook.select('tbestort', 'trial'))

    pool.close()
    
    
#    for i in range(args.n):
#        trialpop, trialstats, trialhof, triallog = ga_func.gaopt_Uni(G, Init_mps, params)
#        
#        #print(triallog)
#        
#        #This un-packs the logbook generated from the optimization performed 
#        # above, by keeping each track of information connected to a trial.
#        logbook.record(trial=i, best=trialhof[0], genmin=triallog.select("min"), 
#                       genmax=triallog.select("max"), genstd=triallog.select("std"), 
#                       genavg=triallog.select("mean"), gengen=triallog.select("gen"))
#        
##    with open(args.oorient, 'w+') as f:
##        json.dump(myhof[0], f)
#    with open(args.output_stats, 'w') as f:
#        json.dump(logbook,f)




# Mouse line:
    #-i ..\\..\\data\\raw\\mouse_tallies -o ..\\..\\data\\interim\\mouse.orient -s ..\\..\\data\\interim\\mouse.stat
    
        
        
        
        
        









