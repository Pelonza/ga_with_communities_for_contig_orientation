"""
By Karl Schmitt, 4/19/17

This script performs basic visualizations of the logbooks from using GA's to 
solve the contig orientation problem.

Dependencies:
    Bokeh
    DEAP
    json
    
"""

import bokeh as bk
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column, gridplot
import bokeh.palettes as pal
from bokeh.palettes import d3
from deap import tools
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
    
    mypal = pal.inferno(20)
    
    f = open('../../data/interim/mouse_cx.stat','r')
    M_cx = json.load(f)
    f.close()
    
    f = open('../../data/interim/mouse_mut.stat','r')
    M_mut = json.load(f)
    f.close()
    
    f = open('../../data/interim/mouse_cx_comm_B.stat','r')
    M_cx_comm = json.load(f)
    f.close()
    
    f = open('../../data/interim/mouse_mut_comm.stat','r')
    M_mut_comm = json.load(f)
    f.close()
    
    f = open('../../data/interim/mouse_mutidpb.stat','r')
    M_mutid = json.load(f)
    f.close()
    
    f = open('../../data/interim/mouse_mutID_comm2.stat','r')
    M_mutid_comm = json.load(f)
    f.close()
    
    
    output_file("test.html")
    cxfig=figure(title="Crossover Parameter Sweep")
    cx_comfig=figure(title="Crossover Parameter Sweep on Communities")
    
    mut_comfig=figure(title="Mutation Parameter Sweep on Communities")
    mut_fig=figure(title="Mutation Parameter Sweep")
    
    mutidfig = figure(title="Independent Mutation Parameter Sweep")
    mutid_comfig = figure(title="Independent Mutation Parameter Sweep on Communities")

    #  Loops over the 20-value parameter sweeps (cx/mut for full/comm)
    for j in range(2,19):
        df = [M_cx[j][k]['tmax'] for k in range(50)]
        avg_tmax = (np.mean(df, axis=0)).tolist()
        cxfig.line(M_cx[0][0]['tgen'], avg_tmax, legend="Crossover = "+
                   str(M_cx[j][0]['tparam'][3]),
                   line_color= d3['Category20'][20][j], muted_alpha=0.2, alpha=1)
        
        dfc = [M_cx_comm[j][k]['tmax'] for k in range(50)]
        avg_tmaxc = (np.mean(dfc, axis=0)).tolist()
        cx_comfig.line(M_cx_comm[0][0]['tgen'], avg_tmaxc, 
                       legend="Crossover = "+str(M_cx_comm[j][0]['tparam'][3]),
                       line_color=d3['Category20'][20][j], muted_alpha=0.2, alpha=1)
        
        df2 = [M_mut[j][k]['tmax'] for k in range(50)]
        avg_tmax2 = (np.mean(df2, axis=0)).tolist()
        mut_fig.line(M_mut[0][0]['tgen'], avg_tmax2, legend="Mutation = "+
                   str(M_mut[j][0]['tparam'][2]),
                   line_color= d3['Category20'][20][j], muted_alpha=0.2, alpha=1)        
        
        dfc2 = [M_mut_comm[j][k]['tmax'] for k in range(50)]
        avg_tmaxc2 = (np.mean(dfc2, axis=0)).tolist()
        mut_comfig.line(M_mut_comm[0][0]['tgen'], avg_tmaxc2, 
                       legend="Mutation = "+str(M_mut_comm[j][0]['tparam'][2]),
                       line_color=d3['Category20'][20][j], muted_alpha=0.2, alpha=1)    

    #  Loops over the 10-value parameter sweeps (independent mutation rates)        
    for j in range(10):
        df = [M_mutid[j][k]['tmax'] for k in range(50)]
        avg_tmax = (np.mean(df, axis=0)).tolist()
        mutidfig.line(M_mutid[0][0]['tgen'], avg_tmax, 
                      legend="Ind. Mut. = "+str(M_mutid[j][0]['tparam'][4]),
                      line_color=d3['Category20'][15][j], alpha=1)
        df2 = [M_mutid_comm[j][k]['tmax'] for k in range(50)]
        
        avg_tmax2 = (np.mean(df2, axis=0)).tolist()
        mutid_comfig.line(M_mutid_comm[0][0]['tgen'], avg_tmax2, 
                      legend="Ind. Mut. = "+str(M_mutid_comm[j][0]['tparam'][4]),
                      line_color=d3['Category20'][15][j], alpha=1)
    
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
    

    #  Move all the legends and set the interactions to hide unwanted lines.    
    cxfig.legend.location = "bottom_right"
    cxfig.legend.click_policy = "hide"
    
    cx_comfig.legend.location = "bottom_right"
    cx_comfig.legend.click_policy = "hide"
    
    mut_comfig.legend.location = "bottom_right"
    mut_comfig.legend.click_policy = "hide"
    
    mut_fig.legend.location = "bottom_right"
    mut_fig.legend.click_policy = "hide"    

    mutidfig.legend.location = "bottom_right"
    mutidfig.legend.click_policy = "hide"

    mutid_comfig.legend.location = "bottom_right"
    mutid_comfig.legend.click_policy = "hide"

    mutid_comm_extr.legend.location = "bottom_right"
    mutid_comm_extr.legend.click_policy = "hide"    
    
    # Plot the two-stage version, GA then Comm
    f = open('../../data/interim/mouse_twostage_test.stat')
    twostage = json.load(f)
    f.close()
    
    gagrp_fig = figure(title = "GA with GA-Comm")
    df_2stage = [ twostage[k][0]['tmax'] for k in range(2)]
    df_2stage_grp = [twostage[k][1]['tmax'] for k in range(2)]
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
    
    show(gridplot([[cxfig, cx_comfig], [mut_fig, mut_comfig], [mutidfig, mutid_comfig], [mutid_extr, mutid_comm_extr], [None, gagrp_fig]]) )  
#    show(column(cxfig, mut_fig, mut_comfig, mutidfig, mutid_comfig) )  