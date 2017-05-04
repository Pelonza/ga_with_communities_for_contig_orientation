# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 08:50:14 2017

@author:Karl R. B. Schmitt


"""
import random

import json
#import dill as pickle
import pickle
import igraph as ig

import numpy as np
import networkx as nx  # Defines G.edges and G.node

from scoop import futures #imports the scoop distributed computing package.


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

def Internal_External(myG, myClusters, CurrOrt):
    #  This function will return the internal/external mate-pairs for
    #  each cluster, based on the input orientation.
    intercluster = myClusters.crossing()
    
    num_cluster = np.max(myClusters.membership)+1
    cluster_score = dict(igood=[False]*num_cluster, ibad=[False]*num_cluster,
                         egood=[False]*num_cluster, ebad=[False]*num_cluster,
                         ifit =[False]*num_cluster, efit=[False]*num_cluster, 
                         icount=[False], ecount=[False])
    for edge in myG.es:
        #  Split inner/outer based on crossing true/false.
        if intercluster[edge.index]: # External Edges
            #  Check orientations to correctly use edge weights.
            cluster_score['ecount'] += edge["w1"] + edge["w2"]
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
            cluster_score['icount'] += edge["w1"]+edge["w2"]
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
    
    cluster_score["emean"] = np.nanmean(cluster_score["efit"])
    cluster_score["imean"] = np.nanmean(cluster_score["ifit"])
        
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
    while (node_net[sort_ord[0]] < 0) :  #  Keep looping till all statified or ort
    # Note: To help this loop condition, nodes with flipped ort will be inf
    
        # Set ort to 'true' (flipped) and unhappy value to 'inf'
        node_net[sort_ord[0]] = np.inf
        node_ort[sort_ord[0]] = True
        
        # Loop through neighbors and update their good-bad values.
        for vertex in G.vs[sort_ord[0]].neighbors():
            indx = vertex.index
            edge = G.get_eid(sort_ord[0], indx)
            node_net[indx] = node_net[indx] + G.es[edge]['w2'] - G.es[edge]['w1']
            
        sort_ord = np.argsort(node_net)


    return node_ort

# %%
if __name__ == "__main__":
    """
    This main tests a few of the functions included in par_main3
    It's also used to test and pickle some of the cluster scores for various 
    turkey runs.
    
    """
    
    ifilename = "../../data/raw/Turkey/orientations_tallies"
    
    G_full, Init_mps_full = load_data(ifilename)
        
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

        
    #perform node-centric orientation
    node_ort = node_centric(G_full)
    nodescr = evaluate(node_ort, G_full, Init_mps_full)
    print("Node Centric Fitness:", nodescr)
    
    #Score base communities.
    null_ort = [False]*G_full.vcount()
    cls_scr = Internal_External(G_full, G_full_clusters, null_ort)
    
    for i in range(len(cls_scr["ifit"])):
        print("Cluster ", i, " I-Fit: ", cls_scr["ifit"][i], " E-Fit: ", cls_scr["efit"][i])
    print("Mean Cluster I-Fit: ", cls_scr['imean']," Mean Cluster E-fit: ", cls_scr['emean'])
    
    
    




# Mouse line:
    #-i ..\\..\\data\\raw\\mouse_tallies -o ..\\..\\data\\interim\\mouse.orient -s ..\\..\\data\\interim\\mouse.stat
    
# Test cluster lines:
# -i ..\\..\\data\\test\\basic2_cluster2.tally -s ..\\..\\data\\interim\\cluster.tally        
        
        
        
        









