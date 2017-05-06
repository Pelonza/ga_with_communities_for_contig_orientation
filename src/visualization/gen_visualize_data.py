"""
By Karl Schmitt, 4/19/17

This script performs basic visualizations of the logbooks from using GA's to 
solve the contig orientation problem.

Dependencies:
    Bokeh
    DEAP
    json

Second visualization file to reduce total number of graphs being plotted per
run to speed up plotting/loading of data.


Files being visualized in script 1:
    mouse_cx.stat    
    mouse_mut.stat
    mouse_cx_comm_B.stat
    mouse_mut_comm.stat
    mouse_mutidpb.stat
    mouse_mutID_comm2
    
Files visualized in script 2:
    mouse_mutidpb_2d5_15_by_2d5_comm.stat
    mouse_mutidpb_2d5_15_by_2d5.stat
    mouse_2xswp_15_30mut_20_30cx.stat
    mouse_2xswp_15_40mut_15_35cx.stat
    
Files visualized in script 3:
    mouse_twostage_full.stat
    mouse_twostage_full_cmga.stat

Files visualized in script 4:
    turkey_longGA_(A/B/C).stat
    turkey_twostage_cmga_(A/B/C).stat
    turkey_twostage_full.stat
    
"""

from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column, gridplot
from bokeh.palettes import d3
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn
import json
import os
import numpy as np
import networkx as nx
import igraph as ig
import pickle
from scoop import futures
    
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
    
def Internal_External(param_tuple):
    #  This function will return the internal/external mate-pairs for
    #  each cluster, based on the input orientation.
    myG = param_tuple[0]
    myClusters = param_tuple[1]
    CurrOrt = param_tuple[2]
        
    intercluster = myClusters.crossing()
    
    num_cluster = np.max(myClusters.membership)+1
    cluster_score = dict(igood=[False]*num_cluster, ibad=[False]*num_cluster,
                         egood=[False]*num_cluster, ebad=[False]*num_cluster,
                         ifit =[False]*num_cluster, efit=[False]*num_cluster, 
                         icount=0, ecount=0)
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

def Internal_External_naive(myG, myClusters):
    #  This function will return the internal/external mate-pairs for
    #  each cluster, based on picking each edges "ideal" orientation
    intercluster = myClusters.crossing()
    
    num_cluster = np.max(myClusters.membership)+1
    cluster_score = dict(igood=[False]*num_cluster, ibad=[False]*num_cluster,
                         egood=[False]*num_cluster, ebad=[False]*num_cluster,
                         ifit =[False]*num_cluster, efit=[False]*num_cluster, 
                         icount=0, ecount=0)
    for edge in myG.es:
        #  Split inner/outer based on crossing true/false.
        if intercluster[edge.index]: # External Edges
            #  Check orientations to correctly use edge weights.
            cluster_score['ecount'] += edge["w1"] + edge["w2"]
            cluster_score["egood"][myClusters.membership[edge.source]]+=np.max([edge["w1"], edge["w2"]])
            cluster_score["egood"][myClusters.membership[edge.target]]+=np.max([edge["w1"], edge["w2"]])
            cluster_score["ebad"][myClusters.membership[edge.source]]+=np.min([edge["w1"], edge["w2"]])
            cluster_score["ebad"][myClusters.membership[edge.target]]+=np.min([edge["w1"], edge["w2"]])
            
        else: # Internal Edges
            #  Check orientations to correctly use edge weights.
            cluster_score['icount'] += edge["w1"]+edge["w2"]
            cluster_score["igood"][myClusters.membership[edge.source]]+=np.max([edge["w1"], edge["w2"]])
            cluster_score["igood"][myClusters.membership[edge.target]]+=np.max([edge["w1"], edge["w2"]])
            cluster_score["ibad"][myClusters.membership[edge.source]]+=np.min([edge["w1"], edge["w2"]])
            cluster_score["ibad"][myClusters.membership[edge.target]]+=np.min([edge["w1"], edge["w2"]])
        
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

    #Initial_matepairs = compute_initial_pairs(ig_G)
        
    return ig_G

def pickle_clsdata(ifile, ofile, ifile2 = None):
    # Little helper function to hide repeated load, open, and pickle lines
    # This also parallizes the computations of the cluster scores!
    
    global G_full
    global G_full_clusters
    with open(ifile,'r') as f:
        df = json.load(f)
    
    if ifile2:
        with open(ifile2,'r') as f:
            df += json.load(f)
            
        ortsdata = [df[k][1] for k in range(len(df))]
    else:
        ortsdata = [df[k][2]['merged_ort'] for k in range(len(df))]

    par_data = zip([G_full]*len(df), [G_full_clusters]*len(df), ortsdata)
    cls_scr = list(futures.map(Internal_External, par_data))
    
    with open(ofile, 'wb') as f:
        pickle.dump(cls_scr, f)
    
    return

# %%
if __name__ == "__main__":

    global G_full
    global G_full_clusters
    
    # Load the turkey orientations and compute clusters for them.
    G_full = load_data('../../data/raw/Turkey/orientations_tallies')
    G_full_dendrogram = G_full.community_fastgreedy(weights="mates")
    G_full_clusters = G_full_dendrogram.as_clustering()
    


    #  Add a preoriented comm-ga point. 
    ifile = '../../data/turkey_prcmga.stat'
    ofile = '../../data/interim/t-prcmga-cls'
    pickle_clsdata(ifile, ofile)
    
#   Some quick code for testing if things got pickled and unpickled correctly
#    with open('../../data/interim/t-prcmga-cls', 'rb') as f:
#        df = pickle.load(f)
#    
#    print(len(df))
#    print(df[0].keys())

    ifile = '../../data/turkey_prgacm.stat'
    ofile = '../../data/interim/t-prgacm-cls'
    pickle_clsdata(ifile, ofile)
    
    ifile = '../../data/turkey_prga_A.stat'
    ifile2 = '../../data/turkey_prga_B.stat'
    ofile = '../../data/interim/t-prga-cls'
    pickle_clsdata(ifile, ofile, ifile2)
    
    ifile = '../../data/turkey_cmga.stat'
    ofile = '../../data/interim/t-cmga-cls'
    pickle_clsdata(ifile, ofile)
    
    ifile = '../../data/turkey_gacm.stat'
    ofile = '../../data/interim/t-gacm-cls'
    pickle_clsdata(ifile, ofile)
    
    ifile = '../../data/turkey_longGA_A.stat'
    ifile2 = '../../data/turkey_longGA_B.stat'
    ofile = '../../data/interim/t-longGA-cls'
    pickle_clsdata(ifile, ofile, ifile2)
    
    # We don't need to parallize the next two since there's only ort each!
    node_ort = node_centric(G_full)
    cls_scr = [Internal_External(list([G_full, G_full_clusters, node_ort]))]
    with open('../../data/interim/t-node-cls','wb') as f:
        pickle.dump(cls_scr,f)
    
    cls_scr = [Internal_External_naive(G_full, G_full_clusters)]
    with open('../../data/interim/t-naive-cls','wb') as f:
        pickle.dump(cls_scr,f)