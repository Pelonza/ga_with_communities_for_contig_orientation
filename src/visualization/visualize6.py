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

def append_data_file(ifile, xdata, ydata, xmstd, ymstd):
    # Little helper function to hide repeated data appends.
    with open(ifile, 'rb') as f:
        cls_scr = pickle.load(f)
        
    xdata += [np.mean([cls_scr[k]['emean'] for k in range(len(cls_scr))])]
    ydata += [np.mean([cls_scr[k]['imean'] for k in range(len(cls_scr))])]
    xmstd += [np.std([cls_scr[k]['emean'] for k in range(len(cls_scr))])/np.sqrt(len(cls_scr))]
    ymstd += [np.std([cls_scr[k]['imean'] for k in range(len(cls_scr))])/np.sqrt(len(cls_scr))]
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
    
    output_file("MultiAlgs_Inn_Out_fit_plot_table.html")

    # Load the turkey orientations and compute clusters for them.
    G_full = load_data('../../data/raw/Turkey/orientations_tallies')
    G_full_dendrogram = G_full.community_fastgreedy(weights="mates")
    G_full_clusters = G_full_dendrogram.as_clustering()
    #G_comm = G_full_clusters.cluster_graph(combine_edges=sum)

    xdata = []
    ydata = []
    xmstd = []
    ymstd = []
    labels = []
    color = []
    
    p = figure(x_axis_label = 'Average External Fitness', y_axis_label = 'Average Internal Fitness')

    #  Add a preoriented comm-ga point. 
    ifile = '../../data/interim/t-prcmga-cls'
    append_data_file(ifile, xdata, ydata, xmstd, ymstd)
    labels += ['Node-Centric with GRP - GA']
    color += [d3['Category20'][20][0]]    

    #  Add a preoriented ga-comm point. 
    ifile = '../../data/interim/t-prgacm-cls'
    append_data_file(ifile, xdata, ydata, xmstd, ymstd)
    labels += ['Node-Centric with GA-GRP']
    color += [d3['Category20'][20][1]]

    #  Add the preoriented ga point.
    ifile = '../../data/interim/t-prga-cls'
    append_data_file(ifile, xdata, ydata, xmstd, ymstd)
    labels += ['Node-Centric with GA']
    color += [d3['Category20'][20][2]]
    
    ifile = '../../data/turkey_cmga.stat'
    ofile = '../../data/t-cmga-cls'
    pickle_clsdata(ifile, ofile)
    
    ifile = '../../data/turkey_gacm.stat'
    ofile = '../../data/t-gacm-cls'
    pickle_clsdata(ifile, ofile)
    
    ifile = '../../data/turkey_longGA_A.stat'
    ifile2 = '../../data/turkey_longGA_B.stat'
    ofile = '../../data/t-longGA-cls'
    pickle_clsdata(ifile, ofile, ifile2)
    
    # We don't need to parallize the next two since there's only ort each!
    node_ort = node_centric(G_full)
    cls_scr = [Internal_External(G_full, G_full_clusters, node_ort)]
    with open('../../data/interim/t-node-cls','wb') as f:
        pickle.dump(cls_scr,f)
    
    cls_scr = [Internal_External_naive(G_full, G_full_clusters)]
    with open('../../data/interim/t-naive-cls','wb') as f:
        pickle.dump(cls_scr,f)


    
    #  Add the Preoriented ga-comm data point.
    df = load('../../data/turkey_prgacm.stat')
        
    cls_scr = [Internal_External(G_full,G_full_clusters, df[k][2]['merged_ort']) for k in range(len(df))]
    append_data(cls_scr, xdata, ydata, xmstd, ymstd)
    labels += ['Node-Centric with GA-Comm']
    color += [d3['Category20'][20][1]]
    

#    
    # Add the Comm-GA point
    df = load('../../data/turkey_cmga.stat')
        
    cls_scr = [Internal_External(G_full,G_full_clusters, df[k][2]['merged_ort']) for k in range(len(df))]
    append_data(cls_scr, xdata, ydata, xmstd, ymstd)
    labels += ['Comm - GA']
    color += [d3['Category20'][20][3]]
    
    #  Add the ga-comm data point.
    df = load('../../data/turkey_gacm.stat')
    
    cls_scr = [Internal_External(G_full,G_full_clusters, df[k][2]['merged_ort']) for k in range(len(df))]
    append_data(cls_scr, xdata, ydata, xmstd, ymstd)
    labels += ['GA-Comm']
    color += [d3['Category20'][20][4]]

    #  Add the preoriented ga point.
    df = load('../../data/interim/turkey_longGA_A.stat')
    df += load('../../data/interim/turkey_longGA_B.stat')
    
    cls_scr = [Internal_External(G_full,G_full_clusters, df[k][1]) for k in range(len(df))]
    append_data(cls_scr, xdata, ydata, xmstd, ymstd)
    labels += ['GA']
    color += [d3['Category20'][20][5]]

    node_ort = node_centric(G_full)
    cls_scr = [Internal_External(G_full, G_full_clusters, node_ort)]
    append_data(cls_scr, xdata, ydata, xmstd, ymstd)
    labels += ['Node-Centric']
    color += [d3['Category20'][20][6]]
    
    cls_scr = [Internal_External_naive(G_full, G_full_clusters)]
    append_data(cls_scr, xdata, ydata, xmstd, ymstd)
    labels += ['Naive Ideal']
    color += [d3['Category20'][20][7]]
    
    source = ColumnDataSource(dict(x=xdata, y=ydata, colors= color, label = labels, xerr = xmstd, yerr = ymstd))
    
    p.circle(x = 'x', y = 'y', color = 'colors', legend = 'label', source = source)
    
    errorbar(p, xdata, ydata, xmstd, ymstd)
    
    p.legend.click_policy = "hide"
    p.legend.location = "bottom_right"
    
    columns = [
            TableColumn(field = "label", title = "Algorithm(s)"),
            TableColumn(field = "x", title = "Ext. Fitness"),
            TableColumn(field = "xerr", title = "E. Fit. Std. Mean Error"),
            TableColumn(field = "y", title = "Int. Fitness"),
            TableColumn(field = "yerr", title = "Int. Fit. Std. Mean Error")]
    
    #source2 = ColumnDataSource(data=dict())
    #source2.data = { 'algs':labels, 'xdata':xdata, 'xmstd':xmstd, 'ydata':ydata, 'ymstd':ymstd}
    
    data_table = DataTable(source=source, columns=columns, row_headers= False, width=800)
    
    
    
    show(column(p, data_table))