# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 08:50:14 2017

@author:Karl R. B. Schmitt

This module seperates out the reading of a tally file and initialization of the
graph data-structure and initial scores. 
"""

import networkx as nx

## NetworkX Initializations and data reading.

#Have networkX use an ordereddict to track the nodes, to guarantee a consistent
# mapping of orientations from individuals to contigs.

def load_data(input_filename):
    from collections import OrderedDict
    class OrderedNodeGraph(nx.Graph):
        node_dict_factory=OrderedDict
    
    #Read in tally file, assuming a 6-column format, no duplicated edges
    G=nx.read_edgelist(input_filename,nodetype=int, data=(('w1',int),('w2',int),('w3',int),('w4',int)), create_using=OrderedNodeGraph())
    
    #Set all edges to 'flipable' and give each node a fixed index for referencing.
    nx.set_node_attributes(G, 'flippable', True)
    i=0
    for n in nx.nodes_iter(G):
        G.node[n]['idx']=i
        i=i+1
    
    #Compute basic data about the tally file/input data.
    Initial_matepairs=list([0,0,0]) #(Total, Good, Bad) Total_matepairs=0
    for (u,v,d) in G.edges(data=True):
        Initial_matepairs[0]=Initial_matepairs[0]+d['w1']+d['w2']+d['w3']+d['w4']
        Initial_matepairs[1]=Initial_matepairs[1]+d['w1']
        Initial_matepairs[2]=Initial_matepairs[2]+d['w2']
    
    return G, Initial_matepairs

def fix_orient(input_filename):
    """
    Load in a list of contigs to not allow flipping of.
    
    This will also require a modification or dedicated evaluation function to 
    account for these fixed orientations.
    """


if __name__=='__main__':
    #Test code for loading data.
    import os
    
    basic_content = '0 1 7 0 0 0\n' \
                    '1 2 4 0 0 0\n' \
                    '1 3 0 6 0 0\n' \
                    '1 4 1 8 0 0\n' \
                    '2 3 1 0 0 0\n' \
                    '3 4 5 1 0 0\n' \
                    '4 2 0 1 0 0\n' \
                    '5 6 2 0 0 0\n'
                    
    fh=open("test.edgelist",'w')
    d=fh.write(basic_content)
    fh.close()
    
    G, Init = load_data('test.edgelist')
    
    assert nx.number_of_nodes(G)==7
    assert nx.number_of_edges(G)==8
    assert Init[0]==36
    assert Init[1]==20
    assert Init[2]==16
    
    #Get edge attributes to test.
    good1=nx.get_edge_attributes(G,'w1')
    bad1=nx.get_edge_attributes(G,'w2')
    good2=nx.get_edge_attributes(G,'w3')
    bad2=nx.get_edge_attributes(G,'w4')
    
    assert good1[(1,2)]==4
    assert bad1[(1,3)]==6
    assert good2[(1,2)]==0
    assert bad2[(1,3)]==0
    
    #Get node attributes to test.
    flip=nx.get_node_attributes(G,'flippable')
    indx=nx.get_node_attributes(G,'idx')
    
    assert flip[0]==True
    assert indx[0]==0
    assert indx[5]==5
                       
    #Test on a default data-file. 
    os.chdir("..\\..\\data\\test")
    
    G, Init = load_data('basic.tally')
    
    assert nx.number_of_nodes(G)==5
    assert nx.number_of_edges(G)==6
    assert Init[0]==32
    assert Init[1]==24
    assert Init[2]==8
    
    print("Tested load_data function successfully")
    
    
    
    