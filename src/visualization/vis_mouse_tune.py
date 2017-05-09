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

from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column, gridplot
from bokeh.palettes import d3
from bokeh.models import Title, Range1d
from bokeh.models.widgets import Panel, Tabs
import json
import numpy as np

# %%
def load(ifile):
    # Little helper function to hide repeated lines for opening files.
    f = open(ifile,'r')
    df = json.load(f)
    f.close()
    
    return df

if __name__ == "__main__":
    
    output_file("mouse_tuning.html")

    # =================
    # Load all the data
    # =================
    
    # Load all the 20-value data files
    M_cx = load('../../data/interim/mouse_cx.stat')
    M_mut = load('../../data/interim/mouse_mut.stat')
    M_cx_comm = load('../../data/interim/mouse_cx_comm_B.stat')
    M_mut_comm = load('../../data/interim/mouse_mut_comm.stat')
       
    # Put sources that have 20 values to plot in list.
    sources_20v = [M_cx, M_mut, M_cx_comm, M_mut_comm]
    
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
    swp2x_C_df = load('../../data/interim/mouse_2xswp_30_50mut_35_50cx.stat')
    swp2x_D_df = load('../../data/interim/mouse_2xswp_high.stat')
    
    sources_other = [swp2x_A_df, swp2x_B_df, swp2x_C_df, swp2x_D_df]
    # ==========================
    # Create some of the figures
    # Also, add sub-titles with fixed parameters.
    # ==========================
    
    cxfig=figure(title="Crossover Parameter Sweep: Mouse")
    cx_comfig=figure(title="Crossover Parameter Sweep on Communities: Mouse")
    
    for src, fig in zip([M_cx, M_cx_comm],[cxfig, cx_comfig]):
        fig.add_layout( Title(text = "Fixed Mutation: " + 
                              str(src[0][0]['tparam'][2]) +
                              " Fixed Ind. Mutation: " + 
                              str(src[0][0]['tparam'][4]), align="center" ), "below")   
    
    mut_comfig=figure(title="Mutation Parameter Sweep on Communities: Mouse")
    mut_fig=figure(title="Mutation Parameter Sweep: Mouse")
    
    for src, fig in zip([M_mut, M_mut_comm],[mut_fig, mut_comfig]):
        fig.add_layout( Title(text = "Fixed Crossover: " + 
                              str(src[0][0]['tparam'][3]) +
                              " Fixed Ind. Mutation: " + 
                              str(src[0][0]['tparam'][4]), align="center" ), "below")
    
    mutidfig = figure(title="Independent Mutation Parameter Sweep: Mouse")
    mutid_comfig = figure(title="Independent Mutation Parameter Sweep on Communities: Mouse")
 
    for src, fig in zip([M_mutid, M_mutid_comm],[mutidfig, mutid_comfig]):
        fig.add_layout( Title(text = "Fixed Crossover: " + 
                              str(src[0][0]['tparam'][3]) +
                              " Fixed Mutation: " + 
                              str(src[0][0]['tparam'][2]), align="center" ), "below")

    swp2x_fig = figure(title = "Sweeping Two Parameters: Mouse")
    swp2xB_fig = figure(title = "Sweeping Two Parameters: Mouse")
    swp2xC_fig = figure(title = "Sweeping Two Parameters: Mouse")

    # Figure lists for iterating over.
    figs_20v = [cxfig, mut_fig, cx_comfig, mut_comfig]    
    figs_10v = [mutid_comfig, mutidfig]
    figs_comm = [cx_comfig, mut_comfig, mutid_comfig]
    figs_full = [cxfig, mut_fig, mutidfig]
    figs_2x = [swp2x_fig, swp2xB_fig, swp2xC_fig]
      
    # =====================
    # Plot the actual data!
    # =====================
    
   
    var_20v = ['Crossover', 'Mutation', 'Crossover', 'Mutation']
    indx_20v = [3, 2, 3, 2]

    for source, fig, var, indx in zip(sources_20v, figs_20v, var_20v, indx_20v):
       for j in range(1,19):
           df = [source[j][k]['tmax'] for k in range(len(source[0]))]
           avg_tmax = (np.mean(df, axis=0)).tolist()
           fig.line(source[0][0]['tgen'], avg_tmax, 
                    legend=var +" = "+ str(source[j][0]['tparam'][indx]),
                    line_color= d3['Category20'][20][j], 
                    muted_alpha=0.2, alpha=1)

    #  Loops over the 10-value parameter sweeps (independent mutation rates)        

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

    for j in range(20):
        df = [swp2x_C_df[j][k]['tmax'] for k in range(len(swp2x_C_df[0][0]))]
        avg_tmax = (np.mean(df, axis=0)).tolist()
        swp2xB_fig.line(swp2x_C_df[0][0]['tgen'], avg_tmax, 
                       legend="Cx. = "+str(swp2x_C_df[j][0]['tparam'][3])
                       +" Mut. = "+str(swp2x_C_df[j][0]['tparam'][2])+
                       " Id. Mut. = "+str(swp2x_C_df[j][0]['tparam'][4]),
                       line_color= d3['Category20'][20][j],
                       muted_alpha=0.2, alpha=1)

    for j in range(12):
        df = [swp2x_D_df[j][k]['tmax'] for k in range(len(swp2x_D_df[0][0]))]
        avg_tmax = (np.mean(df, axis=0)).tolist()
        swp2xC_fig.line(swp2x_D_df[0][0]['tgen'], avg_tmax, 
                       legend="Cx. = "+str(swp2x_D_df[j][0]['tparam'][3])
                       +" Mut. = "+str(swp2x_D_df[j][0]['tparam'][2])+
                       " Id. Mut. = "+str(swp2x_D_df[j][0]['tparam'][4]),
                       line_color= d3['Category20'][20][j],
                       muted_alpha=0.2, alpha=1)

    for fig in figs_full+figs_comm+figs_2x:
        #Set a whole slew of things to make pretty pictures!
        fig.xaxis.axis_label = "Iteration"
        fig.yaxis.axis_label = "Fitness"
        fig.legend.location = "bottom_right"
        fig.legend.click_policy = "hide"
        fig.x_range = Range1d(0, 3000)
        fig.xaxis.bounds = (0, 2500)
        fig.title.align = "center"
        fig.xaxis.axis_label_text_font_size = '16pt'
        fig.yaxis.axis_label_text_font_size = '16pt'
        fig.yaxis.major_label_text_font_size = '16pt'
        fig.xaxis.major_label_text_font_size = '16pt'
        fig.legend.label_text_font_size = '12pt'
        fig.legend.border_line_width = 2
        fig.legend.border_line_alpha = 0.3
        fig.title.text_font_size = '14pt'

    for fig in figs_full:
        fig.y_range = Range1d(0.5, 0.65)
        
    for fig in figs_comm:
        fig.y_range = Range1d(0.65, 0.95)
    
    for fig in figs_2x:
        fig.y_range = Range1d(0.55, 0.8)
        fig.x_range = Range1d(0,4500)

    lay1 = column([cx_comfig, mut_comfig, mutid_comfig])
    lay2 = column([cxfig, mut_fig, mutidfig])
    lay3 = column([swp2x_fig, swp2xB_fig, swp2xC_fig])
    tab1 = Panel(child=lay1, title = "By Communities")
    tab2 = Panel(child=lay2, title = "All Mate-Pairs")
    tab3 = Panel(child=lay3, title = "2 Changing Parameters")
    
    tabs = Tabs(tabs = [tab2, tab1, tab3])

    show(tabs)
    
#    show(gridplot([[cxfig, cx_comfig], [mut_fig, mut_comfig], 
#                   [mutidfig, mutid_comfig],
#                   [swp2x_fig, swp2xB_fig]]) )  
