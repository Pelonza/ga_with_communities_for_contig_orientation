"""
By Karl Schmitt, 4/19/17

This script performs basic visualizations of the logbooks from using GA's to 
solve the contig orientation problem.

Dependencies:
    Bokeh
    DEAP
    json
    
Files being visualized in this script:
    mouse_cx.stat    
    mouse_mut.stat
    mouse_cx_comm_B.stat
    mouse_mut_comm.stat
    mouse_mutidpb.stat
    mouse_mutID_comm2
    
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

# %%
def load(ifile):
    # Little helper function to hide repeated lines for opening files.
    f = open(ifile,'r')
    df = json.load(f)
    f.close()
    
    return df

if __name__ == "__main__":
    
    # =================
    # Load all the data
    # =================
    
    # Load all the 20-value data files
    M_cx = load('../../data/interim/mouse_cx.stat')
    M_mut = load('../../data/interim/mouse_mut.stat')
    M_cx_comm = load('../../data/interim/mouse_cx_comm_B.stat')
    M_mut_comm = load('../../data/interim/mouse_mut_comm.stat')
    swp2x_C_df = load('../../data/interim/mouse_2xswp_30_50mut_35_50cx.stat')
    
    # Put sources that have 20 values to plot in list.
    sources_20v = [M_cx, M_mut, M_cx_comm, M_mut_comm, swp2x_C_df]
    
    # Load all the 10-value date files
    M_mutid = load('../../data/interim/mouse_mutidpb.stat')
    M_mutid_comm = load('../../data/interim/mouse_mutID_comm2.stat')
    
    # Put sources with 10-values together
    sources_10v = [M_mutid, M_mutid_comm]
    
    # Load 6 - value, 2-paired trial data
    More_mutidcomm = load('../../data/interim/mouse_mutidpb_2d5_15_by_2d5_comm.stat')
    More_mutid = load('../../data/interim/mouse_mutidpb_2d5_15_by_2d5.stat')
    
    # Put sources together
    sources_6v = [More_mutidcomm, More_mutid]
    
    # Load 2x Parameter sweeps. These are funny indexing!
    swp2x_A_df = load('../../data/interim/mouse_2xswp_15_30mut_20_30cx.stat')
    swp2x_B_df = load('../../data/interim/mouse_2xswp_15_40mut_15_35cx.stat')
    
    
    # Plan is to create lists of the figures, zip them with list of sources
    # then every loop that has multiple data entries can be condensed. 
    # See sample below:
#    sources = ['mx1', 'mx2']
#
#    figs = ['fig1', 'fig2']
#    
#    for source, fig in zip(sources, figs):
#        print(source)
#        print(fig)

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
    for j in range(1,10):
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
    
 

#    mutid_comm_extr = figure(title="Finer Sweep of Independent Mutation on Communities")
#    mutid_extr = figure(title="Finer Sweep of Independent Mutation")
    for j in range(0,12,2):
        df_extraC = [More_mutidcomm[j][k]['tmax'] for k in range(25)]
        df_extraC = df_extraC + [More_mutidcomm[j+1][k]['tmax'] for k in range(25)]
        avg_tmax_extraC = (np.mean(df_extraC, axis=0)).tolist()
        mutid_comfig.line(More_mutidcomm[0][0]['tgen'], avg_tmax_extraC,
                          legend = "Ind. Mut. = "+str(More_mutidcomm[j][0]['tparam'][4]),
                          line_color=d3['Category20'][12][j], alpha = 1)
        
        df_extra = [More_mutid[j][k]['tmax'] for k in range(25)]
        df_extra = df_extra + [More_mutid[j+1][k]['tmax'] for k in range(25)]
        avg_tmax_extra = (np.mean(df_extra, axis=0)).tolist()
        mutidfig.line(More_mutid[0][0]['tgen'], avg_tmax_extra,
                          legend = "Ind. Mut. = "+str(More_mutid[j][0]['tparam'][4]),
                          line_color=d3['Category20'][12][j], alpha = 1)
 
    
    
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
        
        
#    #  Move all the legends and set the interactions to hide unwanted lines.    
#    mutid_comm_extr.legend.location = "bottom_right"
#    mutid_comm_extr.legend.click_policy = "hide"    
#    
#    mutid_extr.legend.location = "bottom_right"
#    mutid_extr.legend.click_policy = "hide"    

    swp2x_fig.legend.location = "bottom_right"
    swp2x_fig.legend.click_policy = "hide"    

    swp2xB_fig.legend.location = "bottom_right"
    swp2xB_fig.legend.click_policy = "hide"
    
#    #Open additional mutid plots and add them.
#    f = open('../../data/interim/mouse_mutidpb_2d5_15_by_2d5_comm.stat','r')
#    More_mutidcomm = json.load(f)
#    f.close()
#    
#    f = open('../../data/interim/mouse_mutidpb_2d5_15_by_2d5.stat','r')
#    More_mutid = json.load(f)
#    f.close()
#
#
#    mutid_comm_extr = figure(title="Finer Sweep of Independent Mutation on Communities")
#    mutid_extr = figure(title="Finer Sweep of Independent Mutation")
#    for j in range(0,12,2):
#        df_extraC = [More_mutidcomm[j][k]['tmax'] for k in range(25)]
#        df_extraC = df_extraC + [More_mutidcomm[j+1][k]['tmax'] for k in range(25)]
#        avg_tmax_extraC = (np.mean(df_extraC, axis=0)).tolist()
#        mutid_comm_extr.line(More_mutidcomm[0][0]['tgen'], avg_tmax_extraC,
#                          legend = "Ind. Mut. = "+str(More_mutidcomm[j][0]['tparam'][4]),
#                          line_color=d3['Category20'][12][j], alpha = 1)
#        
#        df_extra = [More_mutid[j][k]['tmax'] for k in range(25)]
#        df_extra = df_extra + [More_mutid[j+1][k]['tmax'] for k in range(25)]
#        avg_tmax_extra = (np.mean(df_extra, axis=0)).tolist()
#        mutid_extr.line(More_mutid[0][0]['tgen'], avg_tmax_extra,
#                          legend = "Ind. Mut. = "+str(More_mutid[j][0]['tparam'][4]),
#                          line_color=d3['Category20'][12][j], alpha = 1)
    

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

#    mutid_comm_extr.legend.location = "bottom_right"
#    mutid_comm_extr.legend.click_policy = "hide"    
#    
#    # Plot the two-stage version, GA then Comm
#    f = open('../../data/interim/mouse_twostage_test.stat')
#    twostage = json.load(f)
#    f.close()
#    
#    gagrp_fig = figure(title = "GA with GA-Comm")
#    df_2stage = [ twostage[k][0]['tmax'] for k in range(2)]
#    df_2stage_grp = [twostage[k][1]['tmax'] for k in range(2)]
#    avg_tmax_2stg = (np.mean(df_2stage, axis=0).tolist())
#    avg_tmax_2stg_grp = (np.mean(df_2stage_grp, axis=0).tolist())
#    xdata = twostage[0][0]['tgen']
#    xdata_grp = [ (twostage[0][1]['tgen'][k]+len(xdata)) for k in range(len(twostage[0][1]['tgen'])) ]
#    #xdata = xdata + tmpx
#    gagrp_fig.line(xdata, avg_tmax_2stg,
#                   legend = 
#                   "GA - Mut. = "+str(twostage[0][0]['tparam'][2])+
#                   " , Cx = "+str(twostage[0][0]['tparam'][3]),
#                   line_color = d3['Category20'][3][0])
#    gagrp_fig.line(xdata_grp, avg_tmax_2stg_grp,
#                   legend =
#                   "GA-Comm - Mut. = "+str(twostage[0][1]['tparam'][2])+
#                   " , Cx = "+str(twostage[0][1]['tparam'][3]),
#                   line_color = d3['Category20'][3][2])
    
    show(gridplot([[cxfig, cx_comfig], [mut_fig, mut_comfig], 
                   [mutidfig, mutid_comfig]]) )  
#    show(column(cxfig, mut_fig, mut_comfig, mutidfig, mutid_comfig) )  