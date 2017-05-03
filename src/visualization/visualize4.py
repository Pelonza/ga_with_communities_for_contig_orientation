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
import json
import os
import numpy as np
import networkx as nx
import igraph as ig

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
#    parser.add_argument("-o", "--ofile",
#                        dest="oorient",
#                        help="Write best orientation to FILE",
#                        metavar="FILE")
#    parser.add_argument("-s", "--sfile",
#                        dest="output_stats",
#                        help="Write statistics to FILE",
#                        metavar="FILE")
#    parser.add_argument("-n",
#                        dest="ntrials",
#                        default=50,
#                        type=int,
#                        help="how many trials to do per parameter set")
#    parser.add_argument("-q", "--quiet",
#                        action="store_false",
#                        dest="verbose",
#                        default=True,
#                        help="don't print status messages to stdout")
    return parser

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
    
def Internal_External(myG, myClusters, InOrt):
    #  This function will return the internal/external mate-pairs for
    #  each cluster, based on the input orientation.
    intercluster = myClusters.crossing()
    
    num_cluster = np.max(myClusters.membership)+1
 
    imean = [None]*len(InOrt)
    emean = [None]*len(InOrt)
    
    for j in range(len(InOrt)):
        CurrOrt = InOrt[j]
    
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
        
        imean[j] = np.mean(cluster_score["ifit"])
        emean[j] = np.mean(cluster_score["efit"])     
        
        print("Finished scoring orientation ", str(j))
        
    return list([imean, emean])

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
# %%
if __name__ == "__main__":
    
    
    args = get_parser().parse_args()
    output_file("test4.html")

    # Load the turkey orientations and compute clusters for them.
    G_full = load_data('../../data/raw/Turkey/orientations_tallies')
    G_full_dendrogram = G_full.community_fastgreedy(weights="mates")
    G_full_clusters = G_full_dendrogram.as_clustering()
    #G_comm = G_full_clusters.cluster_graph(combine_edges=sum)

    # Plot the two-stage version, GA then Comm
    f = open('../../data/interim/turkey_twostage_full.stat')
    twostage = json.load(f)
    f.close()
    
    # Get list of GA orientations
    ga_orts = [ twostage[k][0]['tbort'] for k in range(len(twostage))]
    ga_cls_scr = map(Internal_External, ga_orts)  #  This should generate cluster scores for each best orientation
    
    mrg_orts = [ twostage[k][2]['merged_ort'] for k in range(len(twostage))]
    mrg_cls_scr = map(Internal_External, mrg_orts)  # This should generate cluster scores for each merged orientation
    
    
    
    # Note: If I make the ranges dynamic based on the runs, I can reuse these
    # no matter how many trials I do.... 
    gagrp_fig = figure(title = "Turkey GA with GA-Comm")
    df_2stage = [ twostage[k][0]['tmax'] for k in range(50)]
    df_2stage_grp = [twostage[k][1]['tmax'] for k in range(50)]
    avg_tmax_2stg = (np.mean(df_2stage, axis=0).tolist())
    avg_tmax_2stg_grp = (np.mean(df_2stage_grp, axis=0).tolist())
    xdata = twostage[0][0]['tgen']
    xdata_grp = [ (twostage[0][1]['tgen'][k]+len(xdata)) for k in range(len(twostage[0][1]['tgen'])) ]
    #xdata = xdata + tmpx
    gagrp_fig.line(xdata, avg_tmax_2stg,
                   legend = 
                   "GA - Mut. = "+str(twostage[0][0]['tparam'][2])+
                   " , Cx = "+str(twostage[0][0]['tparam'][3]),
                   line_color = d3['Category20'][3][0])
    gagrp_fig.line(xdata_grp, avg_tmax_2stg_grp,
                   legend =
                   "GA-Comm - Mut. = "+str(twostage[0][1]['tparam'][2])+
                   " , Cx = "+str(twostage[0][1]['tparam'][3]),
                   line_color = d3['Category20'][3][2])
    
    gagrp_fig.legend.location = "bottom_right"
    gagrp_fig.legend.click_policy = "mute"
    
    # Plot the two-stage version, Comm then GA
    grpga_fig = figure(title = "Turkey Comm-GA with GA")
        
    f = open('../../data/interim/turkey_twostage_cmga_A.stat')
    twostage = json.load(f)
    f.close()
    
    df_2stage = [ twostage[k][0]['tmax'] for k in range(16)]
    df_2stage_grp = [twostage[k][1]['tmax'] for k in range(16)]
    
    f = open('../../data/interim/turkey_twostage_cmga_B.stat')
    twostage = json.load(f)
    f.close()
    
    df_2stage = df_2stage + [ twostage[k][0]['tmax'] for k in range(16)]
    df_2stage_grp = df_2stage_grp + [twostage[k][1]['tmax'] for k in range(16)]
    
    f = open('../../data/interim/turkey_twostage_cmga_C.stat')
    twostage = json.load(f)
    f.close()
    
    df_2stage = df_2stage + [ twostage[k][0]['tmax'] for k in range(16)]
    df_2stage_grp = df_2stage_grp + [twostage[k][1]['tmax'] for k in range(16)]        

    avg_tmax_2stg = (np.mean(df_2stage, axis=0).tolist())
    avg_tmax_2stg_grp = (np.mean(df_2stage_grp, axis=0).tolist())
    xdata = twostage[0][0]['tgen']
    xdata_grp = [ (twostage[0][1]['tgen'][k]+len(xdata)) for k in range(len(twostage[0][1]['tgen'])) ]
    #xdata = xdata + tmpx
    grpga_fig.line(xdata, avg_tmax_2stg,
                   legend = 
                   "GA - Mut. = "+str(twostage[0][0]['tparam'][2])+
                   " , Cx = "+str(twostage[0][0]['tparam'][3]),
                   line_color = d3['Category20'][3][0])
    grpga_fig.line(xdata_grp, avg_tmax_2stg_grp,
                   legend =
                   "GA-Comm - Mut. = "+str(twostage[0][1]['tparam'][2])+
                   " , Cx = "+str(twostage[0][1]['tparam'][3]),
                   line_color = d3['Category20'][3][2])

    grpga_fig.legend.location = "bottom_right"
    grpga_fig.legend.click_policy = "mute"
    
    f=open('../../data/interim/turkey_longGA_A.stat','r')
    tlonga = json.load(f)
    f.close()
    tmp = [ [tlonga[k][2][j]['max'] for j in range(10001)] for k in range(16)]
    
    f=open('../../data/interim/turkey_longGA_B.stat','r')
    tlongb = json.load(f)
    f.close()
    tmp = tmp + [ [tlongb[k][2][j]['max'] for j in range(10001)] for k in range(16)]
    
    f=open('../../data/interim/turkey_longGA_C.stat','r')
    tlongc = json.load(f)
    f.close()
    tlong = tmp + [ [tlongc[k][2][j]['max'] for j in range(10001)] for k in range(16)]

    tlong_fig = figure(title = "High Generation GA on Turkey")
    avg_max_long = (np.mean(tlong, axis=0).tolist())
    xdata = [ tlonga[0][2][j]['gen'] for j in range(10001)]
    tlong_fig.line(xdata, avg_max_long,
               legend = "Mut. = 0.3" + " , Cx. = 0.9")
    
    tlong_fig.legend.location = "bottom_right"


    show(gridplot([[tlong_fig, None],
                   [gagrp_fig, grpga_fig]]))
#    show(column(cxfig, mut_fig, mut_comfig, mutidfig, mutid_comfig) )  