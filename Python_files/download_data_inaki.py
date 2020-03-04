#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 14:53:11 2020

@author: urzaini1
"""

import os
import networkx as nx
import numpy as np

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
