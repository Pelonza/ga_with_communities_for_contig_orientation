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
    mouse_2xswp_30_50mut_35_50cx.stat
    
Files visualized in script 3:
    mouse_twostage_full.stat
    mouse_twostage_full_cmga.stat
    mouse_longGA.stphy
    
Files to be visualized yet:    
    turkey_longGA_(A/B/C).stat
    turkey_twostage_cmga_(A/B/C).stat
    turkey_twostage_full.stat
    
"""

from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column, gridplot
import bokeh.palettes as pal
from bokeh.palettes import d3
import json
import os
import numpy as np


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
# %%
if __name__ == "__main__":
    
    
    args = get_parser().parse_args()
    output_file("test2.html")
    #Open additional mutid plots and add them.
    f = open('../../data/interim/mouse_mutidpb_2d5_15_by_2d5_comm.stat','r')
    More_mutidcomm = json.load(f)
    f.close()
    
    f = open('../../data/interim/mouse_mutidpb_2d5_15_by_2d5.stat','r')
    More_mutid = json.load(f)
    f.close()

    mutid_comm_extr = figure(title="Finer Sweep of Independent Mutation on Communities")
    mutid_extr = figure(title="Finer Sweep of Independent Mutation")
    for j in range(0,12,2):
        df_extraC = [More_mutidcomm[j][k]['tmax'] for k in range(25)]
        df_extraC = df_extraC + [More_mutidcomm[j+1][k]['tmax'] for k in range(25)]
        avg_tmax_extraC = (np.mean(df_extraC, axis=0)).tolist()
        mutid_comm_extr.line(More_mutidcomm[0][0]['tgen'], avg_tmax_extraC,
                          legend = "Ind. Mut. = "+str(More_mutidcomm[j][0]['tparam'][4]),
                          line_color=d3['Category20'][12][j], alpha = 1)
        
        df_extra = [More_mutid[j][k]['tmax'] for k in range(25)]
        df_extra = df_extra + [More_mutid[j+1][k]['tmax'] for k in range(25)]
        avg_tmax_extra = (np.mean(df_extra, axis=0)).tolist()
        mutid_extr.line(More_mutid[0][0]['tgen'], avg_tmax_extra,
                          legend = "Ind. Mut. = "+str(More_mutid[j][0]['tparam'][4]),
                          line_color=d3['Category20'][12][j], alpha = 1)
    
    # Load 2x Parameter sweeps. These are funny indexing!
    f = open('../../data/interim/mouse_2xswp_15_30mut_20_30cx.stat','r')
    swp2x_A_df = json.load(f)
    f.close()

    f = open('../../data/interim/mouse_2xswp_15_40mut_15_35cx.stat','r')
    swp2x_B_df = json.load(f)
    f.close()

    swp2x_fig = figure(title = "Sweeping Two Parameters - Mouse")
    for j in range(3):
        df = [swp2x_A_df[j][k]['tmax'] for k in range(50)]
        avg_tmax = (np.mean(df, axis=0)).tolist()
        swp2x_fig.line(swp2x_A_df[0][0]['tgen'], avg_tmax, 
                       legend="Cx. = "+str(swp2x_A_df[j][0]['tparam'][3])
                       +" Mut. = "+str(swp2x_A_df[j][0]['tparam'][2])+
                       " Id. Mut. = "+str(swp2x_A_df[j][0]['tparam'][4]),
                       line_color= d3['Category20'][20][j],
                       muted_alpha=0.2, alpha=1)

    for j in range(15):
        df = [swp2x_B_df[j][k]['tmax'] for k in range(50)]
        avg_tmax = (np.mean(df, axis=0)).tolist()
        swp2x_fig.line(swp2x_B_df[0][0]['tgen'], avg_tmax, 
                       legend="Cx. = "+str(swp2x_B_df[j][0]['tparam'][3])
                       +" Mut. = "+str(swp2x_B_df[j][0]['tparam'][2])+
                       " Id. Mut. = "+str(swp2x_B_df[j][0]['tparam'][4]),
                       line_color= d3['Category20'][20][j+3],
                       muted_alpha=0.2, alpha=1)
        
    f = open('../../data/interim/mouse_2xswp_30_50mut_35_50cx.stat','r')
    swp2x_C_df = json.load(f)
    f.close()

    swp2xB_fig = figure(title = "Sweeping Two Parameters B- Mouse")
    for j in range(20):
        df = [swp2x_C_df[j][k]['tmax'] for k in range(len(swp2x_C_df[0][0]))]
        avg_tmax = (np.mean(df, axis=0)).tolist()
        swp2xB_fig.line(swp2x_C_df[0][0]['tgen'], avg_tmax, 
                       legend="Cx. = "+str(swp2x_C_df[j][0]['tparam'][3])
                       +" Mut. = "+str(swp2x_C_df[j][0]['tparam'][2])+
                       " Id. Mut. = "+str(swp2x_C_df[j][0]['tparam'][4]),
                       line_color= d3['Category20'][20][j],
                       muted_alpha=0.2, alpha=1)
        
        
    #  Move all the legends and set the interactions to hide unwanted lines.    
    mutid_comm_extr.legend.location = "bottom_right"
    mutid_comm_extr.legend.click_policy = "hide"    
    
    mutid_extr.legend.location = "bottom_right"
    mutid_extr.legend.click_policy = "hide"    

    swp2x_fig.legend.location = "bottom_right"
    swp2x_fig.legend.click_policy = "hide"    

    swp2xB_fig.legend.location = "bottom_right"
    swp2xB_fig.legend.click_policy = "hide"
    
    show(gridplot([[mutid_extr, mutid_comm_extr],
                   [swp2x_fig, swp2xB_fig]]))
#    show(column(cxfig, mut_fig, mut_comfig, mutidfig, mutid_comfig) )  