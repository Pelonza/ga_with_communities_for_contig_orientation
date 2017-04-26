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
from bokeh.layouts import row, column
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
    
    f = open('../../data/interim/mouse_cx_comm.stat','r')
    M_cx_comm = json.load(f)
    f.close()
    
    f = open('../../data/interim/mouse_mutpb.stat','r')
    M_mut = json.load(f)
    f.close()
    
    f = open('../../data/interim/mouse_mutidpb.stat','r')
    M_mutid = json.load(f)
    f.close()
    
    f = open('../../data/interim/mouse_mutID_comm2.stat','r')
    M_mutid_comm = json.load(f)
    f.close()
    
    
    output_file("test.html")
    cxfig=figure()
    cx_comfig=figure()
    for j in range(2,19):
        df = [M_cx[j][k]['tmax'] for k in range(50)]
        avg_tmax = (np.mean(df, axis=0)).tolist()
        cxfig.line(M_cx[0][0]['tgen'], avg_tmax, legend="Crossover = "+
                   str(M_cx[j][0]['tparam'][3]),
                   line_color= d3['Category20'][20][j], muted_alpha=0.2, alpha=1)
        
        dfc = [M_cx_comm[j][k]['tmax'] for k in range(50)]
        avg_tmaxc = (np.mean(dfc, axis=0)).tolist()
        cx_comfig.line(M_cx_comm[0][0]['tgen'], avg_tmaxc, 
                       legend="Crossover = "+str(M_cx[j][0]['tparam'][3]),
                       line_color=d3['Category20'][20][j], muted_alpha=0.2, alpha=1)
    
    cxfig.legend.location = "bottom_right"
    cxfig.legend.click_policy = "hide"
    
    cx_comfig.legend.location = "bottom_right"
    cx_comfig.legend.click_policy = "hide"
    
    mutidfig = figure()
    mutid_comfig = figure()
    
    for j in range(10):
        df = [M_mutid[j][k]['tmax'] for k in range(50)]
        avg_tmax = (np.mean(df, axis=0)).tolist()
        mutidfig.line(M_mutid[0][0]['tgen'], avg_tmax, 
                      legend="Ind. Mut. = "+str(M_mutid[j][0]['tparam'][4]),
                      line_color=d3['Category10'][10][j], alpha=1)
        df2 = [M_mutid_comm[j][k]['tmax'] for k in range(50)]
        
        avg_tmax2 = (np.mean(df2, axis=0)).tolist()
        mutid_comfig.line(M_mutid_comm[0][0]['tgen'], avg_tmax, 
                      legend="Ind. Mut. = "+str(M_mutid_comm[j][0]['tparam'][4]),
                      line_color=d3['Category10'][10][j], alpha=1)
        
    
    mutidfig.legend.location = "bottom_right"
    mutidfig.legend.click_policy = "hide"
    mutid_comfig.legend.location = "bottom_right"
    mutid_comfig.legend.click_policy = "hide"
    
    show(column(cxfig,cx_comfig, mutidfig, mutid_comfig) )  

    
    
    """
    Process:
        
        load logbook
        df = [ logbook[9][k]['tmax'] for k in range(50)] #Generate sub-data frame
        avg_tmax=np.mean(df, axis=0)  # Average Data
        avg_tmax_lst=avg_tmax.tolist()  #  Un-array the averages
        plot! : p.line(logbook[x][0]['tgen'], avg_tmax_lst)
        
        
        average tmax: avg_max=numpy.mean()
    """