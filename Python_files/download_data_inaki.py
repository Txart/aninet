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
                    
                    # centrality measures
                    eig_c_dict = nx.eigenvector_centrality_numpy(G, weight='weight')
                    eig_c = np.array([eig_c_dict[node] for node in eig_c_dict])
                    bet_c_dict =  nx.betweenness_centrality(G, weight='weight')
                    bet_c = np.array([bet_c_dict[node] for node in bet_c_dict])
                    
                    # exponents
                    expo, inte, r2 = exponent_linear_regression(np.trim_zeros(ego_net, trim='b'))
                    
                    di['graph'] = G
                    di['taxa'] = mydir
                    di['species'] = subdir
                    # name of the species in the other dataset 
                    di['species_name'] = df_info[df_info.graphfile==filename].Species.values[0] if len(df_info[df_info.graphfile==filename].Species.values) != 0 else None
                    di['Interaction type'] = df_info[df_info.graphfile==filename]['Interaction type'].values[0] if len(df_info[df_info.graphfile==filename]['Interaction type'].values) != 0 else None
                    di['group'] = group
                    di['subgroup'] = subgroup # further group segmentation, mentioned above
                    di['filename'] = filename
                    di['undirected_clustering_mean'] = cl_un
                    di['clustering'] = cl
                    di['clustering_mean'] = np.mean(list(cl.values()))
                    di['clustering_std'] = np.std(list(cl.values()))
                    di['ego_net'] = ego_net
                    di['ego_net_std'] = ego_net_std
                    di['ego_net_gini'] = gini(ego_net)
                    di['exponent'] = expo
                    di['intercept'] = inte
                    di['r2'] = r2
                    # if 'tortoiseR_sah'  not in filename: # takes too long (maybe error)
                    #     di['koadrila_tensor'] = None
                    #     di['koadrila_indices'] = None
                    #     di['koadrila_mean'] = None
                    # else:
                    #     k_t, k_i = koadrila(G) # it is not efficient.
                    #     di['koadrila_tensor'] = k_t
                    #     di['koadrila_indices'] = k_i
                    #     di['koadrila_mean'] = np.mean(k_i)
                    mdc = modified_degree_centrality(G)
                    di['mod_deg_centrality'] = mdc
                    di['mod_deg_mean'] = np.mean(mdc)
                    di['deg_cen_gini'] = gini(mdc)
                    di['eigenvector_centrality'] = eig_c
                    di['betweenness_centrality'] = bet_c
                    di['bet_cen_mean'] = np.mean(bet_c)
                    di['bet_cen_gini'] = gini(bet_c)
                    networks.append(di)
     

    """
    Add baboon and primary school children separate data
    """
    # baboons
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
    
    # children
    df_child = read_primary_children_data()
    G_child = get_primary_children_graph(df_child)
    
    # measures
    bab_ego_net, bab_ego_net_std = get_ego_matrix(G_bab)
    child_ego_net, child_ego_net_std = get_ego_matrix(G_child)
    
    cbab = nx.eigenvector_centrality_numpy(G_bab, weight='weight')
    bab_eig_c = np.array([cbab[node] for node in cbab]) # eigenvector centrality
    cchild = nx.eigenvector_centrality_numpy(G_child, weight='weight')
    child_eig_c = np.array([cchild[node] for node in cchild]) # eigenvector centrality
    
    bab_bet_c_dict =  nx.betweenness_centrality(G_bab, weight='weight')
    bab_bet_c = np.array([bab_bet_c_dict[node] for node in bab_bet_c_dict])
    child_bet_c_dict =  nx.betweenness_centrality(G_child, weight='weight')
    child_bet_c = np.array([child_bet_c_dict[node] for node in child_bet_c_dict])
    
    clbab = nx.clustering(G_bab, weight='weight') # Weighted clustering coeff. by Onnela
    clchild = nx.clustering(G_child, weight='weight') # Weighted clustering coeff. by Onnela
    
    mdc_bab = modified_degree_centrality(G_bab)
    mdc_child = modified_degree_centrality(G_child)
    
    expo_bab, inte_bab, r2_bab = exponent_linear_regression(np.trim_zeros(bab_ego_net, trim='b'))
    expo_child, inte_child, r2_child = exponent_linear_regression(np.trim_zeros(child_ego_net, trim='b'))
    
    # baboons 
    networks.append({
        'graph': G_bab,
        'taxa': 'Mammalia',
        'species' : 'baboon_separate_dataset',
        'filename': 'baboon_separate_dataset',
        'group': 'NoGroup',
        'subgroup': 'NoSubgroup',
        'Interaction type' : 'spatial proximity',
        'ego_net': bab_ego_net, # these are set up later
        'ego_net_std':bab_ego_net_std,
        'exponent':expo_bab,
        'intercept' : inte_bab,
        'r2' : r2_bab,
        'eigenvector_centrality': bab_eig_c,
        'betweenness_centrality': bab_bet_c,
        'clustering' : clbab,
        'clustering_mean' : np.mean(list(clbab.values())),
        'clustering_std' : np.std(list(clbab.values())),
        'ego_net_gini' : gini(bab_ego_net),
        'mod_deg_centrality' : mdc_bab,
        'mod_deg_mean' : np.mean(mdc_bab),
        'deg_cen_gini' : gini(mdc_bab),
        'bet_cen_gini' : gini(bab_bet_c),
        'bet_cen_mean' : np.mean(bab_bet_c)
        })
    
    # children
    networks.append({
        'graph': G_child,
        'taxa': 'Mammalia',
        'species' : 'primary_school_children',
        'filename': 'primary_school_children',
        'group': 'NoGroup',
        'subgroup': 'NoSubgroup',
        'Interaction type' : 'spatial proximity',
        'ego_net': child_ego_net, # these are set up later
        'ego_net_std':child_ego_net_std,
        'exponent':expo_child,
        'intercept' : inte_child,
        'r2' : r2_child,
        'eigenvector_centrality': child_eig_c,
        'betweenness_centrality': child_bet_c,
        'clustering' : clchild,
        'clustering_mean' : np.mean(list(clchild.values())),
        'clustering_std' : np.std(list(clchild.values())),
        'ego_net_gini' : gini(child_ego_net),
        'mod_deg_centrality' : mdc_child,
        'mod_deg_mean' : np.mean(mdc_child),
        'deg_cen_gini' : gini(mdc_child),
        'bet_cen_gini' : gini(child_bet_c),
        'bet_cen_mean' : np.mean(child_bet_c)
        })    
    
    
    
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
    
def single_ego_logplot(ego_array, std, title):
    fig,ax = plt.subplots()
    ax.semilogy(ego_array, 'x')
    ax.set_ylim([1e-6, 1])
    ax.set_xlim([0,ego_array.size])
    ax.set_title(title)
    
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

def modified_degree_centrality(G):
    """
    Usual degree centrality (at least the one implemented in networkx)
    for multigraphs is:
        # edges of node / # total possible edges in unweighted network.
    I think that definition is not particularly useful for some of the
    networks we have because if measured for a long time, then every 
    node will have a large degree centrality. The quantity this function
    computes is:
        # edges of node / # total edges in network.
    """
    A = nx.to_numpy_array(G)
    return 2. * np.sum(A, axis=1)/np.sum(A)

def read_primary_children_data():  
    df_child = pd.read_csv('../primary_school/primaryschool.csv', sep='\t',engine='python')
    return df_child

def get_primary_children_graph(df_child):
     # aggregated network. No weights  
    G_child_aggr = nx.from_pandas_edgelist(df=df_child, source='i', target='j')
    
    # Network with weights. Each interaction corresponds to a weight. The weight is saved as an attribute.
    G_child = nx.Graph()
    G_child.add_nodes_from(list(G_child_aggr.nodes))
    for index, row in df_child.iterrows():
        if G_child.has_edge(row.i, row.j):
            G_child[row.i][row.j]['weight'] += 1.
        else:
            G_child.add_edge(row.i, row.j, t=row.t, weight=1.)
    
    return G_child

def gini(c):
    """
    Gini coef. calculation https://en.wikipedia.org/wiki/Gini_coefficient
    Gini = 0 -> most equal distribution
    Gini = 1 -> most unequal distribution

    Parameters
    ----------
    c : 1d numpy array
        centrality measure for each individual

    Returns
    -------
    gini : float
        Gini coefficient

    """
    c_row = np.tile(c, (len(c), 1)) # an array where each row is the c vector
    c_col = c_row.T
    numerator = np.sum(abs(c_col - c_row))
    denominator = 2 * len(c)**2 * np.mean(c)
    gini = numerator / denominator
    
    return gini

#%%

"""
### READ DATABASE AND ADD COLUMNS FROM df_info
"""

df_info, df_nets = read_animals_database()

# Add measurements from df_info
set_info_names = set(df_info.graphfile.to_list())
graphs_in_nets_but_not_in_info = set([row.filename  for i,row in df_nets.iterrows() if row.filename not in set_info_names])

# column names to copy from df_info
di = {'Transitivity': [],
      'Degree assortativity': [],
      'Newman modularity (Q)': [],
      'Network cohesion': []}

for i, row in df_nets.iterrows():
    
    if row.filename in graphs_in_nets_but_not_in_info:
        di['Transitivity'].append(None)
        di['Degree assortativity'].append(None)
        di['Newman modularity (Q)'].append(None)
        di['Network cohesion'].append(None)
    else:
        di['Transitivity'].append(df_info[df_info.graphfile == row.filename].Transitivity.values[0])
        di['Degree assortativity'].append(df_info[df_info.graphfile == row.filename]['Degree assortativity'].values[0])
        di['Newman modularity (Q)'].append(df_info[df_info.graphfile == row.filename]['Newman modularity (Q)'].values[0])
        di['Network cohesion'].append(df_info[df_info.graphfile == row.filename]['Network cohesion'].values[0])


for key in di:
    df_nets.insert(loc=0, column=key, value=di[key], allow_duplicates=False)
    
#%%
"""
SEPARATE BY SPECIES! USEFUL
"""
EGO_SEP_SPECIES = False
# Loop over species. Implemented computations:
##### - Ego network plots (over a hunderd in total)
##### - Linear regression of exponent (assuming ego net plots are exponential)

if EGO_SEP_SPECIES:
    # Read data
    
    
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

#%%
# Baboon time evolution ego plot
def read_baboon_df():
    baboon_fn = os.path.abspath('../baboons/RFID_data/RFID_data.txt') 
    with open(baboon_fn, 'rb'):
        df_bab = pd.read_csv(baboon_fn, delim_whitespace=True)
    return df_bab

def time_evolution_baboon_ego_plot(n_partitions, df_bab):
    """
    Divides baboon data into n_partitions temporal partitions. Returns ego array of each.

    Parameters
    ----------
    n_partitions : int
        1st partition contains the ego array corresponding to the baboon network
        recorded during the first total_length/n_partitions
        2nd partition contains 1st + second

    Returns
    -------
    list with ego array of all partitions

    """
    # aggregated network. No weights  
    G_bab_aggr = nx.from_pandas_edgelist(df=df_bab, source='i', target='j')
        
    experiment_time = df_bab.t.max() - df_bab.t.min()
    partition_time = int(experiment_time / n_partitions)
    end_times = [df_bab.t.min() + partition_time*i for i in range(1, n_partitions+1)]
    ego_net_partitions = []
    
    for end_time in end_times:
        df_bab_partitioned = df_bab[df_bab.t < end_time]
        G_bab = nx.Graph()
        G_bab.add_nodes_from(list(G_bab_aggr.nodes))
        for index, row in df_bab_partitioned.iterrows():
            if G_bab.has_edge(row.i, row.j):
                G_bab[row.i][row.j]['weight'] += 1.
            else:
                G_bab.add_edge(row.i, row.j, t=row.t, Date=row.Date, Time=row.Time, weight=1.)
    
        ego_net, _ = get_ego_matrix(G_bab)
        ego_net_partitions.append(ego_net)
    
    return ego_net_partitions



df_bab = read_baboon_df()
bab_ego_net_time_evolution = time_evolution_baboon_ego_plot(n_partitions=100, df_bab=df_bab)


#%%
fig, ax = plt.subplots(4,4)
for i, axis in enumerate(ax.flatten()):
    axis.bar(x=np.arange(1,len(bab_ego_net_time_evolution[0])+1), height=bab_ego_net_time_evolution[i])
    axis.set_ylim([0, .4])
    axis.set_title(f'cum_part {i+1} baboons')
plt.tight_layout()
    



#%%
# # Ego network and plots
# ego_child, ego_std_child = get_ego_matrix(G_child)
# plot_ego_barplot(ego_child, ego_std_child, 'primary school children')
# single_ego_logplot(ego_child, ego_std_child, 'primary school children logplot')
# exponent_linear_regression(np.trim_zeros(ego_child, trim='b'))

# #%%
# # Baboon entry on the df_nets, for future use:
# dfbab = df_nets[df_nets.species=='baboon_separate_dataset']
# gbab = dfbab.graph.values[0]
# ego_bab, ego_std_bab = get_ego_matrix(gbab)
# exponent_linear_regression(np.trim_zeros(ego_bab, trim='b'))


#%%
"""
centrality measures
"""

# select only spatial proximity and physical contact networks
df_prox = df_info.loc[(df_info['Interaction type']=='physical contact') |
                       (df_info['Interaction type']=='spatial proximity')]
gfn = set(df_prox.graphfile)

# only proximity nets
df_nets_prox = df_nets[df_nets.filename.isin(gfn)]

axes_dict = {'Actinopterygii':(0,0), 'Aves':(0,1), 'Reptilia':(0,2),
             'Insecta':(1,0), 'Mammalia':(1,1)}
# degree centrality
fig, ax = plt.subplots(2,3)
fig.suptitle('degree centrality')
fig_dnied, ax_dnied = plt.subplots(2,3) # dnied: different number of individuals each day
fig_dnied.suptitle('degree centrality')

#eigenvector centrality
fig_ec, ax_ec = plt.subplots(2,3)
fig_ec.suptitle('eigenvector centrality')
fig_ec_dnied, ax_ec_dnied = plt.subplots(2,3) # dnied: different number of individuals each day
fig_ec_dnied.suptitle('eigenvector centrality')

# Betweenness centrality
fig_bc, ax_bc = plt.subplots(2,3)
fig_bc.suptitle('betweenness centrality')
fig_bc_dnied, ax_bc_dnied = plt.subplots(2,3) # dnied: different number of individuals each day
fig_bc_dnied.suptitle('betweenness centrality')

def avg_std_numpy_arrays(a):
    try:
        vsa = np.vstack(a)
    except:
        raise ValueError('array does not have same number of elements, and'
                          + 'cannot be converted into matrix')
     
    vsa_sort = np.sort(vsa, axis=1)
    mean = vsa_sort.mean(axis=0)[::-1]
    std = vsa_sort.std(axis=0)[::-1]
     
    return mean, std
                    
species_prox = set(df_nets_prox.species_name)
for s in species_prox:
    df = df_nets_prox[df_nets_prox.species_name == s]
    if s=='Camponotus fellah' or s=='Papio cynocephalus' or s=='Tursiops truncatus' or s=='Bolitotherus cornutus' or s=='Procyon lotor': # some are special
        for index, row in df.iterrows():
            mod_deg_sorted = np.sort(row.mod_deg_centrality)[::-1]
            eig_sorted = np.sort(row.eigenvector_centrality)[::-1]
            bet_sorted = np.sort(row.betweenness_centrality)[::-1]
            
            ax_dnied[axes_dict[row.taxa]].plot(mod_deg_sorted)
            ax_dnied[axes_dict[row.taxa]].set_title(row.taxa)
            ax_dnied[axes_dict[row.taxa]].set_ylim([0, 1])
            
            ax_ec_dnied[axes_dict[row.taxa]].plot(mod_deg_sorted)
            ax_ec_dnied[axes_dict[row.taxa]].set_title(row.taxa)
            ax_ec_dnied[axes_dict[row.taxa]].set_ylim([0, 1])
            
            ax_bc_dnied[axes_dict[row.taxa]].plot(mod_deg_sorted)
            ax_bc_dnied[axes_dict[row.taxa]].set_title(row.taxa)
            ax_bc_dnied[axes_dict[row.taxa]].set_ylim([0, 1])
            
    else:
        mean_mod_deg, std_mod_deg = avg_std_numpy_arrays(df.mod_deg_centrality.values)
        mean_eig, std_eig = avg_std_numpy_arrays(df.eigenvector_centrality.values)
        mean_bet, std_bet = avg_std_numpy_arrays(df.betweenness_centrality.values)
        
        ax[axes_dict[df.taxa.values[0]]].errorbar(x=list(range(len(mean_mod_deg))), y=mean_mod_deg, yerr=std_mod_deg, fmt='-x', label=df.species.values[0])
        ax[axes_dict[df.taxa.values[0]]].set_title(df.taxa.values[0])
        ax[axes_dict[df.taxa.values[0]]].set_ylim([0, 1])        
        
        ax_ec[axes_dict[df.taxa.values[0]]].errorbar(x=list(range(len(mean_eig))), y=mean_eig, yerr=std_eig, fmt='-x', label=df.species.values[0])
        ax_ec[axes_dict[df.taxa.values[0]]].set_title(df.taxa.values[0])
        ax_ec[axes_dict[df.taxa.values[0]]].set_ylim([0, 1])
        
        ax_bc[axes_dict[df.taxa.values[0]]].errorbar(x=list(range(len(mean_bet))), y=mean_bet, yerr=std_bet, fmt='-x', label = df.species.values[0])
        ax_bc[axes_dict[df.taxa.values[0]]].set_title(df.taxa.values[0])
        ax_bc[axes_dict[df.taxa.values[0]]].set_ylim([0, 1])
        
fig.legend(); fig_ec.legend(); fig_bc.legend()
plt.show()
#%%
# Persistence of signatures
# Following Jari et al.
def H(p):
    """
    Computes the pseudo shannon entropy of p

    Parameters
    ----------
    p : 1d numpy array
        p(r) is fraction of contacts between the ego and the alter of rank r

    Returns
    -------
    H(P) = - sum(from r=1 to k) of p(r) log p(r) 

    """
    return - np.sum(p * np.log(p))
    
    
def jensen_shannon_divergence(P1, P2):
    """
    Parameters
    ----------
    P1 and P2 : numpy arrays
        Social signatures with P = {p(r)} as defined in H.

    """
    return H(.5*(P1 + P2)) - .5*(H(P1) + H(P2))


#%%
# 2D plots - Choose what to plot from each. IMPROVED  BELOW WITH SEABORN PAIRPLOT
    
# import matplotlib.patches as mpatches
# fig,ax = plt.subplots(1)

# taxa_color = {'Actinopterygii': 'blue', 'Aves': 'orange', 'Insecta':'black', 'Mammalia':'red', 'Reptilia':'green' }

# for i, row in df_nets.iloc[::-1].iterrows(): # df.iloc[::-1] reverses iteration order. I do that to make fishes visible.
#     gini_deg = gini(row.mod_deg_centrality)
#     gini_bet = gini(row.betweenness_centrality)
    
#     # get ego graphs and slopes just by row, without grouping in species
#     G = row.graph
#     ego_mean, ego_std = get_ego_matrix(G)
#     slope, intercept, r2 = exponent_linear_regression(np.trim_zeros(ego_mean, trim='b'))
    
#     if r2<1e-3:
#         continue
    
#     ax.scatter(x=slope, y=r2, c=taxa_color[row.taxa], marker='x')
    
# ax.set_xlabel('slope of log ego barplot', fontsize='xx-large')
# ax.set_ylabel('r squared', fontsize='xx-large')

# # Computations for primary school children graph
# child_slope, child_intercept, child_r2 = exponent_linear_regression(np.trim_zeros(ego_child, trim='b'))
# ec = nx.eigenvector_centrality_numpy(G_child, weight='weight')
# eig_centr_child = np.array([ec[node] for node in ec]) # eigenvector centrality
# clust_child = nx.clustering(G_child, weight='weight')
# clust_child_mean = np.mean(list(clust_child.values()))
# bet_c_dict =  nx.betweenness_centrality(G_child, weight='weight')
# bet_centr_child = np.array([bet_c_dict[node] for node in bet_c_dict])
# mod_deg_centr_child = modified_degree_centrality(G_child)

# x_child = child_slope; y_child = child_r2
# ax.scatter(x=x_child, y=y_child, c=taxa_color['Mammalia'])

# # line
# ax.axhline(y=y_child, xmin=0, xmax=1, alpha=0.3, color='grey')
# ax.axvline(x=x_child, ymin=0, ymax=1, alpha=0.3, color='grey')

# # legend
# mam = mpatches.Patch(color=taxa_color['Mammalia'], label='Mammalia')
# act = mpatches.Patch(color=taxa_color['Actinopterygii'], label='Actinopterygii')
# ave = mpatches.Patch(color=taxa_color['Aves'], label='Aves')
# ins = mpatches.Patch(color=taxa_color['Insecta'], label='Insecta')
# rep = mpatches.Patch(color=taxa_color['Reptilia'], label='Reptilia')
# ax.legend(handles=[mam, act, ave, ins, rep])

#%%

#%%
# SEABORN PAIRPLOT
import seaborn as sns
sns.set(style="ticks", color_codes=True)

# select only spatial proximity and physical contact networks
df_prox = df_nets.loc[(df_nets['Interaction type']=='physical contact') |
                       (df_nets['Interaction type']=='spatial proximity')]

columns_to_plot = ['Transitivity', 'clustering_mean', 'Degree assortativity', 'Newman modularity (Q)', 'Network cohesion', 'taxa']
plotdf = df_prox[columns_to_plot].copy()
plotdf.at[787, 'taxa'] = 'School children' # children are not mammals for this plot!

colors = sns.xkcd_palette(["windows blue", "amber", "light brown", "dusty purple", "faded green",  'red'])

g = sns.pairplot(plotdf, hue='taxa', markers=['o','o','o','o','o','s'], palette=colors, corner=True, dropna=True)

#%%
"""
TODO:
    - Problems with ego nets of ants, beetles, baboons, voles and tortoises: measurements of same colony (group) for different days (subgroups) have different amount of individuals, and ego networks cannot be aggregated
    - Children dataset: identify classes
    - Show fitted lines when doing linear regression
    - Evolution of baboon's <-> Jari's paper
    - Centrality measures difference: Gini coefficient
""" 


