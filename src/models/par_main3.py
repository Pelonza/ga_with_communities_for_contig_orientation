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
from copy import deepcopy

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
    parser.add_argument("-t", "--type",
                        dest="type",
                        default=1,
                        help="Optimization Scheme: 1 - Node-Centric, 2 - GA, 3 - GA-Comm, 4 - Comm-GA")
    parser.add_argument("-l", "--loops",
                        dest = "cycle",
                        default = 1,
                        help="How many times to repeat optimization scheme")
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

def Internal_External(myG, myClusters, CurrOrt):
    #  This function will return the internal/external mate-pairs for
    #  each cluster, based on the input orientation.
    intercluster = myClusters.crossing()
    
    num_cluster = np.max(myClusters.membership)+1
    cluster_score = dict(igood=[False]*num_cluster, ibad=[False]*num_cluster,
                         egood=[False]*num_cluster, ebad=[False]*num_cluster,
                         ifit =[False]*num_cluster, efit=[False]*num_cluster)
    for edge in myG.es:
        #  Split inner/outer based on crossing true/false.
        if intercluster[edge.index]: # External Edges
            #  Check orientations to correctly use edge weights.
            if CurrOrt[edge.source] == CurrOrt[edge.target]:
                cluster_score["egood"][myClusters.membership[edge.source]]+=edge["w1"]
                cluster_score["egood"][myClusters.membership[edge.target]]+=edge["w1"]
                cluster_score["ebad"][myClusters.membership[edge.source]]+=edge["w2"]
                cluster_score["ebad"][myClusters.membership[edge.target]]+=edge["w2"]
            else:
                cluster_score["egood"][myClusters.membership[edge.source]]+=edge["w2"]
                cluster_score["egood"][myClusters.membership[edge.target]]+=edge["w2"]
                cluster_score["ebad"][myClusters.membership[edge.source]]+=edge["w1"]
                cluster_score["ebad"][myClusters.membership[edge.target]]+=edge["w1"]                
        else: # Internal Edges
            #  Check orientations to correctly use edge weights.
            if CurrOrt[edge.source] == CurrOrt[edge.target]:
                cluster_score["igood"][myClusters.membership[edge.source]]+=edge["w1"]
                cluster_score["igood"][myClusters.membership[edge.target]]+=edge["w1"]
                cluster_score["ibad"][myClusters.membership[edge.source]]+=edge["w2"]
                cluster_score["ibad"][myClusters.membership[edge.target]]+=edge["w2"]
            else:
                cluster_score["igood"][myClusters.membership[edge.source]]+=edge["w2"]
                cluster_score["igood"][myClusters.membership[edge.target]]+=edge["w2"]
                cluster_score["ibad"][myClusters.membership[edge.source]]+=edge["w1"]
                cluster_score["ibad"][myClusters.membership[edge.target]]+=edge["w1"]                

                
        
    for i in range(num_cluster):
        # Denominators split for line-length.
        tmp_denom = (cluster_score["igood"][i]+cluster_score["ibad"][i])
        if tmp_denom != 0 :
            cluster_score["ifit"][i] = cluster_score["igood"][i]/ tmp_denom
        else:
            cluster_score["ifit"][i] = np.NaN
        
        tmp_denom = (cluster_score["egood"][i]+cluster_score["ebad"][i])
        if tmp_denom != 0 :
            cluster_score["efit"][i] = cluster_score["egood"][i]/ tmp_denom
        else:
            cluster_score["efit"][i] = np.NaN
        
        
            
    return cluster_score
            
def node_centric(G):
    node_net=[False]*G.vcount()
    node_good = [False]*G.vcount()
    node_bad = [False]*G.vcount()
    
    for edge in G.es:
        # Source
        node_good[edge.source] += edge['w1']
        node_bad[edge.source] += edge['w2']
        
        # Target
        node_good[edge.target] +=edge['w1']
        node_bad[edge.target] += edge['w2']
    
    for i in range(G.vcount()):
        node_net[i] = node_good[i] - node_bad[i]

    node_ort = [False]*G.vcount()
    
    # Find the max -> min ordering of the node-happiness
    sort_ord = np.argsort(node_net)
    while (node_net[sort_ord[-1]] < 0) :  #  Keep looping till all statified or ort
    # Note: To help this loop condition, nodes with flipped ort will be inf
    
        # Set ort to 'true' (flipped) and unhappy value to 'inf'
        node_net[sort_ord[-1]] = np.inf
        node_ort[sort_ord[-1]] = True
        
        # Loop through neighbors and update their good-bad values.
        for i in G.vs[sort_ord[-1]].neighbors():
            edge = G.get_eid(sort_ord[-1], i)
            node_net[i] = node_net[i] + G.es[edge]['w2'] - G.es[edge]['w1']
            
        sort_ord = np.argsort(node_net)


    return node_ort
        

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

def one_series_trial_GA_CM(allparam):
    # One trial with GA then Group GA for scoop mapping
    myG = allparam[1]  # Just for easier referencing -- myG should be a node local copy of the graph.
    param_logbook=tools.Logbook()
    
    pop, tbestort, tlogbook = Run_GA(myG, allparam[2])

    param_logbook.record(trial='GA-Full', tmax=tlogbook.select('max'),
                            tbort=tbestort, tgen=tlogbook.select('gen'),
                            tmin=tlogbook.select('min'),
                            tstd=tlogbook.select('std'),
                            tmean=tlogbook.select('mean'),
                            tparam=allparam[2])
    hof_trials.update(pop)

    print("Finished base GA")
    
    # Reorient base graph using best orientation.
    update_graph(myG, tbestort)
    
    full_best_ort = deepcopy(tbestort)  # Preserve best orientation for merging later.
    #tmp_ort = deepcopy(tbestort) # This will be overwritten by unmapping community orientation
    #final_ort = deepcopy(tbestort)  # This will hold the final merged orientation.

    # Generate a reduced by community graph.
    # Cluster membership of a given node can be dereferenced by:
    # *_clusters.membership[node] thus the orientation of a node from cluster
    # is: clus_ort[*_clusters.membership[node]] 
    G_full_dendrogram = myG.community_fastgreedy(weights="mates")
    G_full_clusters = G_full_dendrogram.as_clustering()
    myG_comm = G_full_clusters.cluster_graph(combine_edges=sum)        
    
    pop, tbestort, tlogbook = Run_GA(myG_comm, allparam[3])

    param_logbook.record(trial='GA-Comm', tmax=tlogbook.select('max'),
                            tbort=tbestort, tgen=tlogbook.select('gen'),
                            tmin=tlogbook.select('min'),
                            tstd=tlogbook.select('std'),
                            tmean=tlogbook.select('mean'),
                            tparam=allparam[3])
    hof_trials.update(pop)
    
    # Unmap community orientation
    final_ort=[None]*len(full_best_ort)
    for i in range(len(full_best_ort)):
        final_ort[i] = tbestort[G_full_clusters.membership[i]]^full_best_ort[i]
    
    # Get and store cluster-fitness scores.
    tmp_clusterscr = Internal_External(G_full, G_full_clusters, full_best_ort)
    tmp_gaclsscr = Internal_External(G_full, G_full_clusters, final_ort)
    param_logbook.record(merged_ort = final_ort, 
                         postcommclsscr = tmp_clusterscr,
                         postgaclsscr = tmp_gaclsscr,
                         order = "Comm - GA")
    
    param_logbook.record(merged_ort = final_ort)
    
    print("Finished trial")
    return param_logbook

def one_series_trial_CM_GA(allparam):
    # One trial with GA then Group GA for scoop mapping
    myG = allparam[1]  # Just for easier referencing -- myG should be a node local copy of the graph.
    param_logbook=tools.Logbook()
 
    #  Run optimization on community-grouped graph.
    G_full_dendrogram = myG.community_fastgreedy(weights="mates")
    G_full_clusters = G_full_dendrogram.as_clustering()
    myG_comm = G_full_clusters.cluster_graph(combine_edges=sum)        
    
    pop, tbestort, tlogbook = Run_GA(myG_comm, allparam[3])

    param_logbook.record(trial='GA-Comm', tmax=tlogbook.select('max'),
                            tbort=tbestort, tgen=tlogbook.select('gen'),
                            tmin=tlogbook.select('min'),
                            tstd=tlogbook.select('std'),
                            tmean=tlogbook.select('mean'),
                            tparam=allparam[3])
    hof_trials.update(pop)
    
    print("Finished CM")
    
    #  Get individual contig orientations from group orientations to set graph
    unmap_ort = [None]*myG.vcount()
    for i in range(unmap_ort):
        unmap_ort[i] = tbestort[G_full_clusters.membership[i]]
    
    update_graph(myG, unmap_ort)
    
    
    #  Run optimization on full graph
    pop, tbestort, tlogbook = Run_GA(myG, allparam[2])

    param_logbook.record(trial='GA-Full', tmax=tlogbook.select('max'),
                            tbort=tbestort, tgen=tlogbook.select('gen'),
                            tmin=tlogbook.select('min'),
                            tstd=tlogbook.select('std'),
                            tmean=tlogbook.select('mean'),
                            tparam=allparam[2])
    hof_trials.update(pop)

    print("Finished GA")
    
    
    full_best_ort = deepcopy(tbestort)  # Preserve best orientation for merging later.
    #tmp_ort = deepcopy(tbestort) # This will be overwritten by unmapping community orientation
    #final_ort = deepcopy(tbestort)  # This will hold the final merged orientation.

    # Generate a reduced by community graph.
    # Cluster membership of a given node can be dereferenced by:
    # *_clusters.membership[node] thus the orientation of a node from cluster
    # is: clus_ort[*_clusters.membership[node]] 

    
    # Unmap community orientation
    final_ort=[None]*len(full_best_ort)
    for i in range(len(full_best_ort)):
        final_ort[i] = unmap_ort[i]^full_best_ort[i]
    
    # Merge the unmapped orientation.... except. not.
#    final_ort=[None]*len(full_best_ort)
#    for i in range(len(full_best_ort)):
#        final_ort[i]=full_best_ort[i]^tmp_ort[i]  
    
    tmp_clusterscr = Internal_External(G_full, G_full_clusters, unmap_ort)
    tmp_gaclsscr = Internal_External(G_full, G_full_clusters, final_ort)
    param_logbook.record(merged_ort = final_ort, 
                         postcommclsscr = tmp_clusterscr,
                         postgaclsscr = tmp_gaclsscr,
                         order = "Comm - GA")
    
    
    
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
    params_full = list([50, 1500, 0.30, 0.90, 0.05]) # Not Uniform CX
    params_comm = list([50, 1500, 0.30, 0.70, 0.025])
    
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
    #G_full_dendrogram = G_full.community_fastgreedy(weights="mates")
    #G_full_clusters = G_full_dendrogram.as_clustering()
    #G_comm = G_full_clusters.cluster_graph(combine_edges=sum)

    
    # %%
    
    logbook=tools.Logbook()
    full_logbook=tools.Logbook()  # Logbook for running just the base GA
    
    #Parameter Sweep
    #Base Parameters:
    #param = list([50, 2000, 0.1, 0.1, 0.05, 0.1])
    param = list([50, 3000, 0.1, 0.1, 0.05])

    
    # --------
    # This section runs the various optimization schemes.
    #  Note that when passing the graph we actually pass a copy of it to avoid
    #  referencing/overwrite issues when doing distributed computing.
    # --------
    
    # Pre-declare mapdata as a list of lists, then fill it with data to dist
    mapdata = list(list())
    for optfull in range(args.ntrials):
        if args.type == 1:
            break
        elif args.type == 2:
            mapdata.append(list(G_full.copy(), list(param)))
        elif args.type == 3:
            mapdata.append(list([1, G_full.copy(), list(params_full), list(params_comm)]))
        elif args.type == 4:
            mapdata.append(list([1, G_full.copy(), list(params_full), list(params_comm)]))
        else:
            print("Invalid scheme")
            exit()
    
    #  Use scoop to actually map the data out to the compute nodes.
    if args.type == 1:
        full_logbook = node_centric(G_full)
    elif args.type == 2:
        full_logbook = list(futures.map(Run_GA, mapdata))
    elif args.type == 3:
        full_logbook = list(futures.map(one_series_trial_GA_CM, mapdata))
    elif args.type == 4:
        full_logbook = list(futures.map(one_series_trial_CM_GA, mapdata))


        
    print("Finished trials of parameter") #: ", param[2])
        
    with open(args.output_stats, 'w') as f:
        json.dump(full_logbook,f)

    exit()



# Mouse line:
    #-i ..\\..\\data\\raw\\mouse_tallies -o ..\\..\\data\\interim\\mouse.orient -s ..\\..\\data\\interim\\mouse.stat
    
# Test cluster lines:
# -i ..\\..\\data\\test\\basic2_cluster2.tally -s ..\\..\\data\\interim\\cluster.tally        
        
        
        
        









