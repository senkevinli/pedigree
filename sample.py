#!/usr/bin/env python3
import csv
import networkx as nx
import matplotlib.pyplot as plt

from constructor.pedigree import Node
from os import path

def main():

    # Getting file paths.
    dirname = path.dirname(__file__)
    bios_csv = path.join(dirname, './samples/sample_bio.csv')
    degrees_csv = path.join(dirname, './samples/sample_degrees.csv')

    node_list = []
    with open(bios_csv, 'r') as bios_file:
        reader = csv.reader(bios_file)

        # Assume we know the ordering for now, so skip the first line.
        # TODO: should change later.
        next(reader)

        for row in reader:
            female = row[1] == 'F'
            node_list.append(Node(
                row[0], female, row[2], row[3] if not female else None, row[4]
            ))
    for node in node_list:
        print(node)

def visualize_graph(nodes):
    
    G = nx.Graph()
    G.add_nodes_from(['A', 'B', 'C'])
    G.add_edge('A', 'B')
    G.add_edge('B', 'C')
    values = [0.25, 0.25, 0.25]
    edge_colors = [0.3, 1]
    pos = {
        'A': [0, 0],
        'B': [0.5, 0],
        'C': [1, 0]
    }
    nx.draw(G, pos, cmap=plt.get_cmap('jet'), edge_color=edge_colors, node_color=values, node_size=900)

    plt.axis('off')
    plt.show()
if __name__ == '__main__':
    # main()
    visualize_graph(None)

