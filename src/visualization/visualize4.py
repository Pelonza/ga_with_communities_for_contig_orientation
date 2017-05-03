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
    output_file("test4.html")

    # Plot the two-stage version, GA then Comm
    f = open('../../data/interim/mouse_twostage_full.stat')
    twostage = json.load(f)
    f.close()
    
    gagrp_fig = figure(title = "GA with GA-Comm")
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
    f = open('../../data/interim/mouse_twostage_full_cmga.stat')
    twostage = json.load(f)
    f.close()
    
    grpga_fig = figure(title = "Comm-GA with GA")
    df_2stage = [ twostage[k][0]['tmax'] for k in range(50)]
    df_2stage_grp = [twostage[k][1]['tmax'] for k in range(50)]
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
    
    f=open('../../data/interim/mouse_longGA.stphy','r')
    mlong = json.load(f)
    f.close()
    
    mlong_fig = figure(title = "High Generation GA on Mouse")
    dflong = [mlong[k][0]['tmax'] for k in range(50)]
    avg_max_long = (np.mean(dflong, axis=0).tolist())
    mlong_fig.line(mlong[0][0]['tgen'], avg_max_long,
               legend = "Mut. = "+str(mlong[0][0]['tparam'][2]) +
               " , Cx. = "+str(mlong[0][0]['tparam'][3]))
    
    mlong_fig.legend.location = "bottom_right"
    
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