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
import numpy as np
import networkx as nx
import igraph as ig
import pickle

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

def load(ifile):
    # Little helper function to hide repeated lines for opening files.
    f = open(ifile,'r')
    df = json.load(f)
    f.close()
    
    return df

def append_data(cls_scr, xdata, ydata, xmstd, ymstd):
    # Little helper function to hide repeated data appends.
    xdata += [np.mean([cls_scr[k]['emean'] for k in range(len(cls_scr))])]
    ydata += [np.mean([cls_scr[k]['imean'] for k in range(len(cls_scr))])]
    xmstd += [np.std([cls_scr[k]['emean'] for k in range(len(cls_scr))])/np.sqrt(len(cls_scr))]
    ymstd += [np.std([cls_scr[k]['imean'] for k in range(len(cls_scr))])/np.sqrt(len(cls_scr))]
    return

def append_data_file(ifile, xdata, ydata, xmstd, ymstd, tdata, tmstd):
    # Little helper function to hide repeated data appends.
    with open(ifile, 'rb') as f:
        cls_scr = pickle.load(f)
        
    xdata += [np.mean([cls_scr[k]['emean'] for k in range(len(cls_scr))])]
    ydata += [np.mean([cls_scr[k]['imean'] for k in range(len(cls_scr))])]
    xmstd += [np.std([cls_scr[k]['emean'] for k in range(len(cls_scr))])/np.sqrt(len(cls_scr))]
    ymstd += [np.std([cls_scr[k]['imean'] for k in range(len(cls_scr))])/np.sqrt(len(cls_scr))]
    
    #  Compute a bunch of intermediate values
    tigs = [np.sum(x) for x in [cls_scr[k]['igood'] for k in range(len(cls_scr))]]
    tegs = [np.sum(x) for x in [cls_scr[k]['egood'] for k in range(len(cls_scr))]]
    tebs = [np.sum(x) for x in [cls_scr[k]['ebad'] for k in range(len(cls_scr))]]
    tibs = [np.sum(x) for x in [cls_scr[k]['ibad'] for k in range(len(cls_scr))]]
    tgood = np.sum([tigs, tegs], axis = 0)/2  #  total good
    tmps = np.sum([tigs, tegs, tibs, tebs], axis = 0 )/2  #  total matepairs
    aryfits = tgood/tmps
    tdata += [np.mean(aryfits)]
    tmstd += [np.std(aryfits)/np.sqrt(len(aryfits))]

    return

def errorbar(fig, x, y, xerr=None, yerr=None, color='red', 
             point_kwargs={}, error_kwargs={}):

  fig.circle(x, y, color=color, **point_kwargs)

  if xerr:
      x_err_x = []
      x_err_y = []
      for px, py, err in zip(x, y, xerr):
          x_err_x.append((px - err, px + err))
          x_err_y.append((py, py))
      fig.multi_line(x_err_x, x_err_y, color=color, **error_kwargs)

  if yerr:
      y_err_x = []
      y_err_y = []
      for px, py, err in zip(x, y, yerr):
          y_err_x.append((px, px))
          y_err_y.append((py - err, py + err))
      fig.multi_line(y_err_x, y_err_y, color=color, **error_kwargs)
      
# %%
if __name__ == "__main__":
    
    output_file("test-MultiAlgs_Inn_Out_fit_plot_table.html")

    # Load the turkey orientations and compute clusters for them.
    G_full = load_data('../../data/raw/Turkey/orientations_tallies')
    G_full_dendrogram = G_full.community_fastgreedy(weights="mates")
    G_full_clusters = G_full_dendrogram.as_clustering()
    #G_comm = G_full_clusters.cluster_graph(combine_edges=sum)

    xdata = []
    ydata = []
    tdata = []
    tmstd = []
    xmstd = []
    ymstd = []
    labels = []
    color = []
    
    p = figure(x_axis_label = 'Average External Fitness', y_axis_label = 'Average Internal Fitness')

    #  Add a preoriented comm-ga point. 
    ifile = '../../data/interim/t-prcmga-cls'
    append_data_file(ifile, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    labels += ['Node-Centric with GRP - GA']
    color += [d3['Category20'][20][0]]    

    #  Add a preoriented ga-comm point. 
    ifile = '../../data/interim/t-prgacm-cls'
    append_data_file(ifile, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    labels += ['Node-Centric with GA-GRP']
    color += [d3['Category20'][20][1]]

    #  Add the preoriented ga point.
    ifile = '../../data/interim/t-prga-cls'
    append_data_file(ifile, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    labels += ['Node-Centric with GA']
    color += [d3['Category20'][20][2]]
    
    # Add the Comm-GA point
    ifile = '../../data/interim/t-cmga-cls'
    append_data_file(ifile, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    labels += ['Comm - GA']
    color += [d3['Category20'][20][3]]

    #  Add the ga-comm data point.
    ifile = '../../data/interim/t-gacm-cls'
    append_data_file(ifile, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    labels += ['GA-Comm']
    color += [d3['Category20'][20][4]]

    #  Add the preoriented ga point.
    ifile = '../../data/interim/t-longGA-cls'
    append_data_file(ifile, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    labels += ['GA']
    color += [d3['Category20'][20][5]]


    # Add the node-centric point
    ifile = '../../data/interim/t-node-cls'
    append_data_file(ifile, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    xmstd[-1:] = ' '  #  Since only 1 orientation, 'error/std' is meaningless
    ymstd[-1:] = ' '  #  Since only 1 orientation, 'error/std' is meaningless
    labels += ['Node-Centric']
    color += [d3['Category20'][20][6]]
    
#    # Add the naive point.    
#    ifile = '../../data/interim/t-naive-cls'
#    append_data_file(ifile, xdata, ydata, xmstd, ymstd, tdata, tmstd)
#    xmstd[-1:] = ' '  #  Since only 1 orientation, 'error/std' is meaningless
#    ymstd[-1:] = ' '  #  Since only 1 orientation, 'error/std' is meaningless
#    labels += ['Naive Ideal']
#    color += [d3['Category20'][20][7]]

    # Turns out the error bars are TINY! ... so don't display them. Included
    # line and function though incase other runs/trials have larger errorbars.
#    errorbar(p, xdata, ydata, xmstd, ymstd)
    source = ColumnDataSource(dict(x=xdata, y=ydata, colors= color, 
                                   label = labels, xerr=xmstd, yerr=ymstd,
                                   t = tdata, terr = tmstd))
    
    p.circle(x = 'x', y = 'y', color = 'colors', legend = 'label', source = source, size = 15, muted_alpha = 0.2)

    p.legend.click_policy = "mute"
    p.legend.location = "bottom_right"
    p.title.text = 'Average Interior vs. Exterior Fitness Per Optimization'
    
    columns = [
            TableColumn(field = "label", title = "Algorithm(s)"),
            TableColumn(field = "x", title = "Ext. Fitness"),
            TableColumn(field = "xerr", title = "E. Fit. Std. M. Err"),
            TableColumn(field = "y", title = "Int. Fitness"),
            TableColumn(field = "yerr", title = "Int. Fit. Std. M. Err"),
            TableColumn(field = "t", title = "Avg. Fitness"),
            TableColumn(field = "terr", title = "Fit. Std. M. Err")]
    
    #source2 = ColumnDataSource(data=dict())
    #source2.data = { 'algs':labels, 'xdata':xdata, 'xmstd':xmstd, 'ydata':ydata, 'ymstd':ymstd}
    
    data_table = DataTable(source=source, columns=columns, row_headers= False, width=800)
    
    
    
    show(column(p, data_table))