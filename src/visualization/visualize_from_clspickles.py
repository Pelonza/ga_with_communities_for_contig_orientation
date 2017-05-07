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
from bokeh.models.widgets.tables import NumberFormatter
import json
import numpy as np
import networkx as nx
import igraph as ig
import pickle
import scipy.stats as scs

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
    aryfits = compute_totalfit(cls_scr)
    tdata += [np.mean(aryfits)]
    tmstd += [np.std(aryfits)/np.sqrt(len(aryfits))]

    return

def compute_totalfit(cls_scr):
    tigs = [np.sum(x) for x in [cls_scr[k]['igood'] for k in range(len(cls_scr))]]
    tegs = [np.sum(x) for x in [cls_scr[k]['egood'] for k in range(len(cls_scr))]]
    tebs = [np.sum(x) for x in [cls_scr[k]['ebad'] for k in range(len(cls_scr))]]
    tibs = [np.sum(x) for x in [cls_scr[k]['ibad'] for k in range(len(cls_scr))]]
    tgood = np.sum([tigs, tegs], axis = 0)/2  #  total good
    tmps = np.sum([tigs, tegs, tibs, tebs], axis = 0 )/2  #  total matepairs
    aryfits = tgood/tmps
    
    return aryfits
    
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

def mannw(files, labels):
    # this produces 3 matrices of mannwhitneyu tests.
    cls_scrs = []
    for fname in files:
        with open(fname, 'rb') as f:
            cls_scr += [pickle.load(f)]
    
    mwu_e = np.zeros((len(files), len(files)))
    mwu_i = np.zeros((len(files), len(files)))
    mwu_t = np.zeros((len(files), len(files)))

    for i in range(len(files)):
        for j in range(len(files)):
            x = [cls_scrs[i][k]['emean'] for k in range(len(cls_scrs[i]))]
            y = [cls_scrs[j][k]['emean'] for k in range(len(cls_scrs[j]))]
            mwu_e[i,j] = scs.mannwhitneyu(x,y, alternative = 'two-sided')
            x = [cls_scrs[i][k]['imean'] for k in range(len(cls_scrs[i]))]
            y = [cls_scrs[j][k]['imean'] for k in range(len(cls_scrs[j]))]
            mwu_i[i,j] = scs.mannwhitneyu(x,y, alternative = 'two-sided')
            c_i_fits = compute_totalfit(cls_scrs[i])
            c_j_fits = compute_totalfit(cls_scrs[j])
            mwu_t[i,j] = scs.mannwhitneyu(c_i_fits, c_j_fits, alternative = 'two-sided')
            
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
    ifile2 = '../../data/interim/t-prgacm-cls'
    append_data_file(ifile2, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    labels += ['Node-Centric with GA-GRP']
    color += [d3['Category20'][20][1]]

    #  Add the preoriented ga point.
    ifile3 = '../../data/interim/t-prga-cls'
    append_data_file(ifile3, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    labels += ['Node-Centric with GA']
    color += [d3['Category20'][20][2]]
    
    # Add the Comm-GA point
    ifile4 = '../../data/interim/t-cmga-cls'
    append_data_file(ifile4, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    labels += ['Comm - GA']
    color += [d3['Category20'][20][3]]

    #  Add the ga-comm data point.
    ifile5 = '../../data/interim/t-gacm-cls'
    append_data_file(ifile5, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    labels += ['GA-Comm']
    color += [d3['Category20'][20][4]]

    #  Add the preoriented ga point.
    ifile6 = '../../data/interim/t-longGA-cls'
    append_data_file(ifile6, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    labels += ['GA']
    color += [d3['Category20'][20][5]]


    # Add the node-centric point
    ifile7 = '../../data/interim/t-node-cls'
    append_data_file(ifile7, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    xmstd[-1:] = ' '  #  Since only 1 orientation, 'error/std' is meaningless
    ymstd[-1:] = ' '  #  Since only 1 orientation, 'error/std' is meaningless
    labels += ['Node-Centric']
    color += [d3['Category20'][20][6]]
    
    # Add the naive point.    
    ifile8 = '../../data/interim/t-naive-cls'
    append_data_file(ifile8, xdata, ydata, xmstd, ymstd, tdata, tmstd)
    xmstd[-1:] = ' '  #  Since only 1 orientation, 'error/std' is meaningless
    ymstd[-1:] = ' '  #  Since only 1 orientation, 'error/std' is meaningless
    labels += ['Naive Ideal']
    color += [d3['Category20'][20][7]]
    
    mannw_results = mannw(list([ifile, ifile2, ifile3, ifile4, ifile5, ifile6]), labels[:-3])

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
            TableColumn(field = "x", title = "Ext. Fitness", formatter = NumberFormatter(format = '0.000[00]')),
            TableColumn(field = "xerr", title = "E. Fit. Std. M. Err", formatter = NumberFormatter(format = '0.000[00]')),
            TableColumn(field = "y", title = "Int. Fitness", formatter = NumberFormatter(format = '0.000[00]')),
            TableColumn(field = "yerr", title = "Int. Fit. Std. M. Err", formatter = NumberFormatter(format = '0.000[00]')),
            TableColumn(field = "t", title = "Avg. Fitness", formatter = NumberFormatter(format = '0.000[00]')),
            TableColumn(field = "terr", title = "Fit. Std. M. Err", formatter = NumberFormatter(format = '0.000[00]'))]
    
    #source2 = ColumnDataSource(data=dict())
    #source2.data = { 'algs':labels, 'xdata':xdata, 'xmstd':xmstd, 'ydata':ydata, 'ymstd':ymstd}
    
    data_table = DataTable(source=source, columns=columns, row_headers= False, width=800)
    
    
    
    show(column(p, data_table))