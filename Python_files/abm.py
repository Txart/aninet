# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 14:57:50 2020

@author: 03125327
"""
import numpy as np
import scipy.spatial as spa
import matplotlib.pyplot as plt

"""
PARAMETERS
"""
N_AGENTS = 20
SQUARE_DIM = 10

DETECTION_R = 0.5 # Every contact with agent below this radius is recorded 
MAX_V = 2. # max speed of agents

ITER = 10000

def _boundary_conditions(x, y, d=SQUARE_DIM):
    if x < 0:
        x_new = x + d
    elif x > d:
        x_new = x - d
    else:
        x_new = x
    if y < 0:
        y_new = y + d
    elif y > d:
        y_new = y - d
    else:
        y_new = y
        
    return x_new, y_new

def get_ego(A):
    norm_factor = np.multiply(np.sum(A,axis=0), np.ones(shape=A.shape)).transpose()
    A_norm = np.divide(A, norm_factor) # vector type normalization to save time
    A_norm_ordered = np.flip(np.sort(A_norm), axis=-1) # order each row from highest to lowest. 
    group_mean_array = np.mean(A_norm_ordered, axis=0)
    group_std_array = np.std(A_norm_ordered, axis=0)
    
    return group_mean_array, group_std_array

def plot_ego_barplot(ego_array, std, title):
    plt.figure()
    plt.bar(x=np.arange(1,ego_array.shape[0]+1), height=ego_array, yerr=std)
    plt.title(title)
    
    
class Agent:
    def __init__(self, speed, posx, posy, max_turn_angle):
        self.v = speed
        self.x = posx
        self.y = posy
        self.theta_max = max_turn_angle # In radians
    
    def move_one_step(self):
        theta = self.theta_max * np.random.random()
        self.x, self.y = _boundary_conditions(self.x + self.v * np.cos(theta), self.y + self.v * np.sin(theta))
    
# Create fishes
rand = np.random.random
fishes = [Agent(speed=rand() * MAX_V, posx=rand()*SQUARE_DIM, posy=rand()*SQUARE_DIM, max_turn_angle= rand()*2*np.pi) for i in range(N_AGENTS)]

# Create measurement observation matrices
proximity_adj_matrix = np.zeros(shape=(N_AGENTS, N_AGENTS))

for i in range(ITER):
    # move
    for fish in fishes:
        fish.move_one_step()
    
    positions_array = np.array([np.array([fish.x, fish.y]) for fish in fishes])

    dist_matrix = spa.distance_matrix(positions_array, positions_array)
    
    proximity_adj_matrix += (dist_matrix < DETECTION_R).astype(int)
    

np.fill_diagonal(proximity_adj_matrix, 0) # Remove self-loops

        
# ego data
ego_mean, ego_std = get_ego(proximity_adj_matrix)
plot_ego_barplot(ego_mean, ego_std, title='random agents')

            