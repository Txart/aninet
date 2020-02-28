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
for mydir in directories:
    num_nodes_list = []
    num_edges_list = []
    net_density_list = []
    avg_degree_list = []
    cv_degree_list = []
    modularity_list = []
    qmax_list = []
    qrel_list = []
    cohesion_list = []
    diameter_list = []
    asrt_list = []
    avg_betw_list=[]
    betw_wt_list = []
    clstr_list = []
    clstr_wt_list = []
		
    val_list = [num_nodes_list, num_edges_list, net_density_list, avg_degree_list, cv_degree_list, asrt_list, avg_betw_list, betw_wt_list, clstr_list, clstr_wt_list, modularity_list, qmax_list, qrel_list, cohesion_list, diameter_list]

    sub_directories = os.listdir(os.path.abspath('../Networks/'+mydir))
    sub_directories = [name for name in sub_directories if name != '.DS_Store']
    for subdir in sub_directories:
        print mydir, subdir
        
        for filename in sorted(os.listdir(os.path.abspath('../Networks/'+mydir+'/'+subdir))):    
            if filename.endswith(".graphml"):
                print ("filename ............"), filename
                
                G= nx.read_graphml(os.path.abspath('../Networks/'+mydir+'/'+subdir+'/'+ filename))
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
                num_nodes_list.append(n_nodes)
                num_edges_list.append(n_edges)
                is_connect = nx.is_connected(G) 
                num_comp = nx.number_connected_components(G)
                density = round(nx.density(G),3)
                net_density_list.append(density)
                deg_list = [G.degree(node) for node in list(G.nodes)]
                avg_deg = round(np.mean(deg_list),3)
                avg_degree_list.append(avg_deg)
                std_deg = round(np.std(deg_list),3)

                if avg_deg>0:
                    cv_deg = round(float(std_deg)/avg_deg,3)
                    cv_degree_list.append(cv_deg)
                else: cv_deg = "NA"
        
                try:
                    asrt = round(nx.degree_assortativity_coefficient(G),3)
                    asrt_list.append(asrt)
                except: "asrt not calculated"
                #asrt_wt = round(nx.degree_assortativity_coefficient(G, weight="weight"),3)
                betw_list  = nx.betweenness_centrality(G).values()
                betw_wt_list  = nx.betweenness_centrality(G, weight="weight").values()
                avg_betw =  round(np.mean(betw_list) ,3)
                avg_betw_list.append(avg_betw)
                betw_wt =  round(np.mean(betw_wt_list) ,3)
                betw_wt_list.append(avg_betw)
                clstr = round(nx.average_clustering(G),3)
                clstr_list.append(clstr)
                clstr_wt = round(nx.average_clustering(G, weight="weight"),3)
                clstr_wt_list.append(clstr_wt)
        
                #for the rest of the computations, network is required to be connected
                if not nx.is_connected(G): G = max(nx.connected_component_subgraphs(G), key=len)
                G1 = nx.Graph()
                G1.add_nodes_from(G.nodes)
                G1.add_edges_from(G.edges)
                try:
                    partition = community.best_partition(G1)
                    Q = round(community.modularity(partition, G1),3)
                    modularity_list.append(Q)
                    modules = list(set(partition.values()))
                    mod_nodes= {}
                    for mod in modules: mod_nodes[mod] = [node for node in list(G1.nodes) if partition[node]==mod]
                    Qmax = round(calculate_Qmax(G1, mod_nodes),3)
                    qmax_list.append(Qmax)
                    coh = calculate_avg_wd(G1, partition, n_nodes)/(1.*avg_deg)
                    cohesion_list.append(coh)
                except: "modularity not calculated"
                
                
                diam = nx.diameter(G)
                avg_modsize = float(len(list(G1.nodes)))/len(modules)
                if Qmax>0:Qrel = round(float(Q)/Qmax,3)
                else: Qrel="NA"    
                qrel_list.append(Qrel)