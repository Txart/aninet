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



def _group_segmentation(fn):
    """

    Parameters
    ----------
    fn : string
        Name of the file (.graphml) containing the data.

    Returns
    -------
    g : string
        Group membership of datafile. Groups are same species but different relevant attribute, e.g., male/female, location, community, etc.
        Several graphfiles go often into the same group
        Different groups often involve different individuals, different number of individuals, etc.
    sg : string
        Subgroup. Further compartimentalization into subgroups.

    """
    g = 'NoGroup' ; sg = 'NoSubgroup'
    
    if 'geese' in fn:
        if 'female' in fn: g = 'female'
        else: g = 'male'
    elif 'arnberg_sparrow' in fn:
        if '2009' in fn: g = '2009'
        else: g = '2010'
    elif 'shizuka' in fn:
        if 'Season2' in fn: g = 'Season2'
        else: g = 'Season3'
    elif 'weaver' in fn: g = fn.split(sep='_')[-1][:2]
    elif 'wildbird' in fn: g = fn.split(sep='_')[-1][0]
    elif 'ant_mersch' in fn: g = fn.split(sep='_')[3]; sg = fn.split(sep='_')[4][:5]
    elif 'ant_quevillon' in fn: g = fn.split(sep='_')[3]; sg = fn.split(sep='_')[4][:4]
    elif 'beetle' in fn: g = fn.split(sep='_')[3]; sg = fn.split(sep='_')[4] + fn.split(sep='_')[5][0]
    elif 'baboon_franz' in fn: g = fn.split(sep='_')[2]; sg = fn.split(sep='_')[4][:2]
    elif 'gazda' in fn: g = fn.split(sep='_')[2]
    elif 'elephantseal_' in fn: g = fn.split(sep='_')[5] + '-' + fn.split(sep='_')[6]
    elif 'hyena' in fn: g = fn.split(sep='_')[2][:8]
    elif 'griffin_primate' in fn: g = fn.split(sep='_')[4][:2]
    elif 'raccoon' in fn: g = fn.split(sep='_')[3][:2]
    elif 'Macaque' in fn: g = fn.split(sep='_')[1]; sg = fn.split(sep='_')[2]
    elif 'voles' in fn: g = fn.split(sep='_')[2] + fn.split(sep='_')[3]; sg = fn.split(sep='_')[4] + fn.split(sep='_')[5] + fn.split(sep='_')[6][:2]
    elif 'tortoise' in fn: g = fn.split(sep='_')[3]; sg = fn.split(sep='_')[6] + fn.split(sep='_')[7] + fn.split(sep='_')[8]
    
    return g, sg


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
                    
                    G_un = G.to_undirected()
                    
                    # Group and subgroup differentiation                    
                    group, subgroup = _group_segmentation(filename)
                    ego_net, ego_net_std = get_ego_matrix(G) # Some are empty, so normalization throws error. Check!!
                    cl = nx.clustering(G, weight='weight') # Weighted clustering coeff. by Onnela
                    cl_un = nx.average_clustering(G_un)
                    
                    di['graph'] = G
                    di['taxa'] = mydir
                    di['species'] = subdir
                    di['group'] = group
                    di['subgroup'] = subgroup # further group segmentation, mentioned above
                    di['filename'] = filename
                    di['undirected_clustering_mean'] = cl_un
                    di['clustering'] = cl
                    di['clustering_mean'] = np.mean(list(cl.values()))
                    di['clustering_std'] = np.std(list(cl.values()))
                    di['ego_net'] = ego_net
                    di['ego_net_std'] = ego_net_std
                    di['exponent'] = None
                    di['intercept'] = None
                    di['r2'] = None
                    if 'tortoiseR_sah'  not in filename: # takes too long (maybe error)
                        di['koadrila_tensor'] = None
                        di['koadrila_indices'] = None
                        di['koadrila_mean'] = None
                    else:
                        k_t, k_i = koadrila(G) # it is not efficient.
                        di['koadrila_tensor'] = k_t
                        di['koadrila_indices'] = k_i
                        di['koadrila_mean'] = np.mean(k_i)
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
    
    ego_net, ego_net_std = get_ego_matrix(G_bab)
    networks.append({
        'graph': G_bab,
        'taxa': 'Mammalia',
        'species' : 'baboon_separate_dataset',
        'filename': 'baboon_separate_dataset',
        'group': 'NoGroup',
        'subgroup': 'NoSubgroup',
        'ego_net': None, # these are set up later
        'ego_net_std':None,
        'exponent':None})               
    
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
    group_sum_array : np.array
        1-d sorted vector. Each entry corresponds to the fraction of the weight of the ego-network of the species.
    group_std_array: np.array
        1-d vector of std for the species_sum_array
        
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
    A_norm = np.nan_to_num(A_norm)
    A_norm_ordered = np.flip(np.sort(A_norm), axis=-1) # order each row from highest to lowest. 
    group_mean_array = np.mean(A_norm_ordered, axis=0)
    group_std_array = np.std(A_norm_ordered, axis=0)
    
    return group_mean_array, group_std_array

def plot_ego_barplot(ego_array, std, title):
    plt.figure()
    plt.bar(x=np.arange(1,ego_array.shape[0]+1), height=ego_array, yerr=std)
    plt.title(title)
    
def plot_ego_log_scale(ego_array, label, taxa):
    if taxa == 'Actinopterygii':
        positioner = (0,0)
    elif taxa == 'Aves':
        positioner = (0,1)
    elif taxa == 'Reptilia':
        positioner = (0,2)
    elif taxa == 'Insecta':
        positioner = (1,0)
    elif taxa == 'Mammalia':
        positioner = (1,1)
        
    taxa_color = {'Actinopterygii': 'blue', 'Aves': 'orange', 'Insecta':'black', 'Mammalia':'red', 'Reptilia':'green' }
    ax[positioner].semilogy(ego_array, 'x', label=taxa, color=taxa_color[taxa], alpha=0.9 if taxa=='Actinopterygii' else 0.3)
    ax[positioner].set_ylim([1e-6, 1])
    ax[positioner].set_xlim([0,110])
    ax[positioner].set_title(taxa)
    
def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))
    
def plot_exponents(df_nets):
    plt.figure()
    plt.title('exponents')
    plt.scatter(x = df_nets[df_nets['exponent'] < 0]['taxa'], y = df_nets[df_nets['exponent'] < 0]['exponent'], alpha = 0.7, marker='x')

    
    
def exponent_linear_regression(ego_values):
    """
    Assumes exponential decrease in ego_values and performs linear regression of the log
    (i.e., computes the exponent)
    
    Parameters
    ----------
    ego_values : np.array
        Values of the ego barplot for each alter

    Returns
    -------
    slope, intercept, r_squared: float
        Usual linear regression coeffs. 'slope' is the best-fitting exponent. 

    """

    from sklearn.linear_model import LinearRegression
    
    x = np.array(range(0, ego_values.size)).reshape((-1,1)) 
    y = np.log(ego_values)
    
    model = LinearRegression().fit(x,y)
    
    slope = model.coef_[0]
    intercept = model.intercept_
    r_squared = model.score(x, y)
    
    return slope, intercept, r_squared


# Made-up koadrila coefficient tensor and index
def koadrila(G):
    
    W = nx.to_numpy_array(G) # Weighted adjacency matrix
    K = np.zeros(shape=(W.shape[0], W.shape[0], W.shape[0])) # 3d tensor
    tol = 1e-3
    for (i,j,k), _ in np.ndenumerate(K):
        if W[i,j] < tol or W[i,k] < tol:
            K[i,j,k] = -1 # Not defined for those 
        if i != j and i != k:
            K[i,j,k] = 2 * W[j,k] / (W[i,j] + W[i,k])
        else:
            continue # defaults to 0
    
    koadrila_index = [0] * W.shape[0]
    for i, Ki in enumerate(K):
        koadrila_index[i] = np.mean(Ki[Ki>0])
    
    return np.array(K), np.array(koadrila_index)

"""
#####
"""
# Read data
df_info, df_nets = read_animals_database()

# Generate ego matrices and store them in dataframe
# for index, row in df_nets.iterrows():
    
#     df_nets['ego_net'].values[index] = ego_net
#     df_nets['ego_net_std'].values[index] = ego_net_std
#     try:
        
#     except:
#         print('problem with row number: ',index)
#     # df_nets['koadrila_tensor'].values[index] = k_t
#     df_nets['koadrila_indices'].values[index] = k_i

# Loop over species. Implemented computations:
##### - Ego network plots (over a hunderd in total)
##### - Linear regression of exponent (assuming ego net plots are exponential)

PLOT_EGO = True

plt.close('all')
species_labels = df_nets.species.unique()
fig,ax = plt.subplots(2,3)

for s in species_labels:
    if 'ants' in s or 'beetle' in s or 'baboon_association' in s or 'voles' in s or 'tortoise' in s: # Skip ants, beetles, baboons, voles and tortoises for now. Problem: measurements of same colony (group) for different days (subgroups) have different amount of individuals
        continue
    print(s)
    df_species = df_nets[df_nets.species == s]
    group_labels = df_species.group.unique()
    for group in group_labels:
        df_group = df_species[df_species.group == group]
        ego_nets = df_group.ego_net.to_numpy()
        ego_net_stds = df_group.ego_net_std.to_numpy() # individual stds
        ego_net_stds_vstack = np.vstack(ego_net_stds) # These 3 lines solve problem with nparray formatting of dataframe
        if len(ego_net_stds_vstack.shape) > 1:
            ego_net_stds_vstack = ego_net_stds_vstack[0]
        try:
            ego_mean = np.mean(ego_nets, axis=0)
            ego_std_group = np.std(ego_nets, axis=0) # Careful! There are two std, 1) at the level of single-measurement, 2) at the level of groups (more than one measurements, more than one graphfiles)!
            
            # Linear regression
            expo, inter, r2 = exponent_linear_regression(np.trim_zeros(ego_mean,trim='b'))
            
            df_nets.loc[df_nets['group'] == group, ['exponent']] = expo
            df_nets.loc[df_nets['group'] == group, ['intercept']] = inter
            df_nets.loc[df_nets['group'] == group, ['r2']] = r2
            
            # Plot
            if PLOT_EGO:
                
                if not np.any(ego_std_group): # if std at the group level is zero, i.e., if there are no groups
                    # plot_ego_barplot(ego_mean, std=ego_net_stds_vstack, title=s + group)
                    plot_ego_log_scale(ego_mean, label=s+group, taxa=df_species.taxa.values[0] )
                else: # if more than one group, plot group-wise std. Here group means  a bunch of measurements of a same species with same relevant attributes.
                    # plot_ego_barplot(ego_mean, std=ego_std_group, title=s + group)
                    plot_ego_log_scale(ego_mean, label=s+group, taxa=df_species.taxa.values[0])
            else:
                continue
        except:
            print('>>>>>>>> ERROR in {} species. Probably, there is a mismatch in number of individuals accross measurements'.format(s))
            continue

# Linear regression graph by graph for excluded species
# And plot of insecta regression
f, ax = plt.subplots()
ax.set_ylim([1e-6, 1])
ax.set_xlim([0,110])
ax.set_title('Reptilia') 
for index, row in df_nets.iterrows():
    if row.taxa == 'Reptilia':
        ax.semilogy(row.ego_net, 'x', label='Reptilia', color='green', alpha=0.3)
    s = row.species
    if 'ants' in s or 'beetle' in s or 'baboon_association' in s or 'voles' in s or 'tortoise' in s:
        print(s)
        ego_net = row.ego_net
        expo, inter, r2 = exponent_linear_regression(np.trim_zeros(ego_net, trim='b'))
        df_nets.loc[index,'exponent'] = expo
        df_nets.loc[index,'intercept'] = inter
        df_nets.loc[index,'r2'] = r2

# legend_without_duplicate_labels(ax)
        
# plot the ones excluded above one by one
# ant_colony = df_nets[df_nets['group'] == 'colony6']
# elephant_seal = df_nets[df_nets['species']=='elephantseal_dominance_weighted']
# beetle_C1 = df_nets[df_nets['group']=='C1']
# a = ant_colony
# for day, i in enumerate(a.ego_net.to_numpy()):
#     plot_ego_barplot(i, std=0, title='Ant_colony6'.format(day))

# 


# plot exponent values
plot_exponents(df_nets)


     


 


"""
TODO:
    - Manage plots better. Maybe save figs.
    - Sometimes the normalization is a division by zero. Check!
    - Problems with ego nets of ants, beetles, baboons, voles and tortoises: measurements of same colony (group) for different days (subgroups) have different amount of individuals, and ego networks cannot be aggregated

""" 


