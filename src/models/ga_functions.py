# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 08:57:38 2017

@author: Karl R.B. Schmitt

This script contains definitions of all the functions used within the DEAP 
package for GA's. This includes unique implementations of crossovers, evaluate
etc.
"""

#Can try commenting this out and see if function still works!
import networkx as nx  #Defines G.edges and G.node 

#Define an fitness-evaluation function.
def evaluate(individual, G, Init_pairs):
    #Inputs are one individual (with array of flip/no-flip orientations), and 
    # a Graph structure (with nodes/edges from tally file). It returns the 
    # fitness of the individual.       
    Orient_pairs=list(Init_pairs)
    for (u,v,d) in G.edges(data=True):
        if individual[G.node[u]['idx']]!=individual[G.node[v]['idx']]:
            Orient_pairs[1]=Orient_pairs[1]-d['w1']+d['w2']
            Orient_pairs[2]=Orient_pairs[2]-d['w2']+d['w1']
    
    return (Orient_pairs[1]/Orient_pairs[0], )

if __name__=='__main__':
    """
    Test Evaluate
    
    Depends on module: load_data being
    """
    
    import sys
    import os
    
    #Inserts the folder with the load-data module into the sys-path
    #Note this is NOT good practice! I should fix it to be an actual package...
    sys.path.insert(0,'..\\data') 
    
    import load_data as ld
    
    os.chdir("..\\..\\data\\test")
    
    G, Init = ld.load_data('basic.tally')
    
    #Check to make sure the test file loaded looks like we expect.
    assert nx.number_of_nodes(G)==5
    assert nx.number_of_edges(G)==6
    assert Init[0]==32
    assert Init[1]==24
    assert Init[2]==8
    
    #Basic orientations possible in an individual
    no_flip=list([0,0,0,0,0])
    one_flip=list([0,0,0,0,1])
    mny_flip=list([1,1,1,1,0])
    
    a = evaluate(no_flip, G, Init)
    
    assert a[0]==0.75
    #assert a[1]==None
    
    a = evaluate(one_flip, G, Init)
    
    assert a[0]==0.84375
    #assert b==None
    
    a = evaluate(mny_flip, G, Init)
    
    assert a[0]==0.84375
    #assert b==None
    
    #Check for correct clearing of repeated evaluates
    a = evaluate(one_flip, G, Init)
    
    assert a[0]==0.84375
    #assert b==None
    
    a = evaluate(mny_flip, G, Init)
    
    assert a[0]==0.84375
    #assert b==None
    
    #Confirm our initial data hasn't gotten modified
    assert nx.number_of_nodes(G)==5
    assert nx.number_of_edges(G)==6
    assert Init[0]==32
    assert Init[1]==24
    assert Init[2]==8
    
    print("Testing of evaluate in ga_functions has succeeded")