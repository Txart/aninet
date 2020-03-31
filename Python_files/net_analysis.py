import networkx as nx
from networkx.algorithms import community
import matplotlib.pyplot as plt
import re 
import sys
import os
import subprocess


def analyse_graph(path):
    G = nx.Graph()   
    with open(path, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines)):
            if ("node" in lines[i]):
                try:
                    node = re.search('\"(.+?)\"', lines[i]).groups()[0]
                except:
                    node = "unknown"
                    pass
                #print(node)
                G.add_node(node)
                pass
            if ("source" in lines[i]):
                result = re.findall('\"(.+?)\"', lines[i])
                source = result[0]
                target = result[1]
                weighted_line = lines[i+1]
                #print(weighted_line)
                try:
                    weight = re.search('>(.+?)<', weighted_line).groups()[0]
                except:
                    weight = 1
                    pass
                #print(source+" "+target)
                G.add_edge(source, target, weight=weight)
                pass
        
        commu = community.label_propagation_communities(G)
        sorted_commu = sorted(map(sorted, commu))
        number_of_nodes = G.number_of_nodes()
        number_of_edges = G.number_of_edges()
        number_of_commu = len(sorted_commu)
        ratio_nc = number_of_nodes/number_of_commu
        ratio_ec = number_of_edges/number_of_commu
        ratio_ne = number_of_nodes/number_of_edges
        print("Number of nodes : "+str(number_of_nodes))
        print("Number of edges : "+str(number_of_edges))
        print("Number of communities : "+str(number_of_commu))
        print("Nodes / Communities ratio = "+str(ratio_nc))
        print("Edges / Communities ratio = "+str(ratio_ec))
        print("Nodes / Edges = "+str(ratio_ne))
        with open("data.log", "a") as log:
            log.write(path.split('/')[-1]+" "+str(number_of_nodes)+" "+str(number_of_edges)+" "+str(number_of_commu)+" "+str(ratio_nc)+" "+str(ratio_ec)+" "+str(ratio_ne)+"\n")

        #top_level_commu = next(commu)
        #next_level_commu = next(commu)
        print("\nCommunities:")
        for com in sorted_commu:
            print(sorted(com))
            pass


        # Get current size
        fig_size = plt.rcParams["figure.figsize"]
        # Set figure width to 12 and height to 9
        fig_size[0] = 12
        fig_size[1] = 9
        plt.rcParams["figure.figsize"] = fig_size
        plt.subplot(111)
        plt.title(path.split('/')[-1])
        color_map = []
        for node in G:
            for i, com in enumerate(sorted_commu):
                if node in com:
                    color_map.append(20*i)
        nx.draw(G, node_color=color_map, with_labels=True, font_weight='bold')
        
        #plt.show()

root = subprocess.run(['pwd'], stdout=subprocess.PIPE).stdout.decode('utf-8')
print(root.split('/')[1:-1])
file_path = ""
for w in root.split('/')[1:-1]:
    file_path = file_path+"/"+w
print(file_path)
try:
    file_path = sys.argv[1]
except:
    pass
files_in_dir = []
for r, d, f in os.walk(file_path):
    #print(d)
    for item in f:
        if 'graphml' in item:
            files_in_dir.append(os.path.join(r, item))

with open("data.log", "w") as log:
    pass
for item in sorted(files_in_dir):
    print("Graph under dir: ", item)
    analyse_graph(item)

