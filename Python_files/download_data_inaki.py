#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 14:53:11 2020

@author: urzaini1
"""

import os
import networkx as nx
import numpy as np
import pandas as pd

directories = os.listdir(os.path.abspath('../Networks/'))
directories = [name for name in directories if name != '.DS_Store']

networks = {}

for mydir in directories:
    sub_directories = os.listdir(os.path.abspath('../Networks/'+mydir))
    sub_directories = [name for name in sub_directories if name != '.DS_Store']
    for subdir in sub_directories:
        print(mydir, subdir)
        
        for filename in sorted(os.listdir(os.path.abspath('../Networks/'+mydir+'/'+subdir))):    
            if filename.endswith(".graphml"):
                print ("filename ............", filename)
                
                G = nx.read_graphml(os.path.abspath('../Networks/'+mydir+'/'+subdir+'/'+ filename))
                G.remove_edges_from(nx.selfloop_edges(G))                
                n_nodes = len(list(G.nodes))
                n_edges = len(list(G.edges))
                
#####
                ## if network does not have weights, add a weight of one to all edges
                if len(nx.get_edge_attributes(G,'weight'))==0:
                    for (n1,n2) in list(G.edges): G[n1][n2]['weight']=1
                ####

                ####################################################
                #if no edges then return NAs
                if n_edges==0:continue
                ########################################################
                ## remove edges with weight zero
                for (n1, n2) in list(G.edges):
                    if G[n1][n2]['weight']==0: 
                        #print ("Removing NULL edge!!!"), (n1, n2), len(G.edges()),
                        G.remove_edge(n1,n2)
                        #print len(G.edges())
        
                #########################################################
                        
                networks[filename] = G

"""
Baboons
"""

baboon_fn = os.path.abspath('../baboons/RFID_data/RFID_data.txt') 

with open(baboon_fn, 'rb'):
    df_bab = pd.read_csv(baboon_fn, delim_whitespace=True)
# aggregated network. No weights  
G_bab_aggr = nx.from_pandas_edgelist(df=df_bab, source='i', target='j')

# Network with weights. Each interaction corresponds to a weight.
G_bab = nx.Graph()
G_bab.add_nodes_from(list(G_bab_aggr.nodes))
# G_bab.add_edges_from(list(G_bab_aggr.edges))
for index, row in df_bab.iterrows():
    if G_bab.has_edge(row.i, row.j):
        G_bab[row.i][row.j]['weight'] += 1.
    else:
        G_bab.add_edge(row.i, row.j, t=row.t, Date=row.Date, Time=row.Time, weight=1.)

# Normalize weights
all_weights = [data['weight'] for i,j,data in G.edges(data=True)]
