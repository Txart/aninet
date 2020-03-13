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
import matplotlib.pyplot as plt


def read_animals_database():
    directories = os.listdir(os.path.abspath('../Networks/'))
    directories = [name for name in directories if name != '.DS_Store']
    
    df_info = pd.read_csv('Network_summary_master_file.csv', engine='python')
    
    networks = []
    
    for mydir in directories:
        sub_directories = os.listdir(os.path.abspath('../Networks/'+mydir))
        sub_directories = [name for name in sub_directories if name != '.DS_Store']
        for subdir in sub_directories:
            print(mydir, subdir)
            
            for filename in sorted(os.listdir(os.path.abspath('../Networks/'+mydir+'/'+subdir))):
                di = {}
                if filename.endswith(".graphml"):
                    print ("filename ............", filename)
                    
                    G = nx.read_graphml(os.path.abspath('../Networks/'+mydir+'/'+subdir+'/'+ filename))
                    G.remove_edges_from(nx.selfloop_edges(G))                
                    n_nodes = len(list(G.nodes))
                    n_edges = len(list(G.edges))
                    
                    ## if network does not have weights, add a weight of one to all edges
                    if len(nx.get_edge_attributes(G,'weight'))==0:
                        for (n1,n2) in list(G.edges): G[n1][n2]['weight']=1
    
                    #if no edges then return NAs
                    if n_edges==0:continue

                    ## remove edges with weight zero
                    for (n1, n2) in list(G.edges):
                        if G[n1][n2]['weight']==0: 
                            G.remove_edge(n1,n2)
                            
                    di['graph'] = G
                    di['taxa'] = mydir
                    di['species'] = subdir
                    di['filename'] = filename
                    di['ego_net'] = get_ego_matrix(G) # Some are empty, so normalization throws error. Check!!
                    networks.append(di)
     

    """
    Add baboon separate data
    """
    baboon_fn = os.path.abspath('../baboons/RFID_data/RFID_data.txt') 
    
    with open(baboon_fn, 'rb'):
        df_bab = pd.read_csv(baboon_fn, delim_whitespace=True)
    # aggregated network. No weights  
    G_bab_aggr = nx.from_pandas_edgelist(df=df_bab, source='i', target='j')
    
    # Network with weights. Each interaction corresponds to a weight. The weight is saved as an attribute.
    G_bab = nx.Graph()
    G_bab.add_nodes_from(list(G_bab_aggr.nodes))
    for index, row in df_bab.iterrows():
        if G_bab.has_edge(row.i, row.j):
            G_bab[row.i][row.j]['weight'] += 1.
        else:
            G_bab.add_edge(row.i, row.j, t=row.t, Date=row.Date, Time=row.Time, weight=1.)
    
    networks.append({
        'graph': G_bab,
        'taxa': 'Mammalia',
        'species' : 'baboon_separate_dataset',
        'filename': 'baboon_separate_dataset',
        'ego_net': get_ego_matrix(G_bab)})               
    
    df_networks = pd.DataFrame(networks)
                    
    return df_info, df_networks


def _check_symm_matrix(A, atol=1e-5, rtol=1e-8):
    return np.allclose(A, A.transpose(), rtol=rtol, atol=atol)

def get_ego_matrix(graph):
    """
    Parameters
    ----------
    graph : NetworkX graph

    Raises
    ------
    ValueError
        if the adjacency of the matrix is not square.
        
    Returns
    -------
    species_sum_array : np.array
        1-d sorted vector. Each entry corresponds to the fraction of the weight of the ego-network of the species.
        
    Note
    -------
    A_norm_ordered is a 2d numpy array whose rows are the individuals ego networks.

    """
    A = nx.to_numpy_array(graph) # weighted adjacency matrix
    if A.shape[0] != A.shape[1]:
        raise ValueError('Not a square adjacency matrix')
    
    if not _check_symm_matrix(A): # matrix is antisymmetric
        print('Non-symmetric matrix. Forcing symmetry.')
        A = A + A.transpose() # i->j interactions are the same as j->i interactions
    
    norm_factor = np.multiply(np.sum(A,axis=0), np.ones(shape=A.shape)).transpose()
    A_norm = np.divide(A, norm_factor) # vector type normalization to save time
    A_norm_ordered = np.flip(np.sort(A_norm), axis=-1) # order each row from highest to lowest. 
    species_sum_array = np.sum(A_norm_ordered, axis=0) /  A.shape[0] # mean
    
    return species_sum_array

def plot_ego_barplot(ego_array, title):
    plt.figure()
    plt.bar(x=np.arange(1,ego_array.shape[0]+1), height=ego_array)
    plt.title(title)


"""
#####
"""
df_info, df_nets = read_animals_database()

# Ego plots aggregated accross species
species_labels = df_nets.species.unique()
plots = []
for s in species_labels:
    print(s)
    ego_nets = df_nets[df_nets.species == s].ego_net.to_numpy()
    try:
        ego_mean = np.sum(ego_nets, axis=0)/len(ego_nets)
        plots.append(plot_ego_barplot(ego_mean, title=s))
    except:
        print('>>>>>>>> ERROR in {} species. Probably, there is a mismatch in number of individuals accross measurements'.format(s))
        continue
        

"""
TODO:
    - Manage plots better. Maybe save figs.
    - Some plots are identical! E.g., asianelephants and baboons. WHY??
    - Right now, not separation of some data. E.g., ant colony_1 and ant_colony_2 are lumped into the same graph.
""" 



