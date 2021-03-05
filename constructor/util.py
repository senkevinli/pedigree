#!/usr/bin/env python3

""" For visualizing the constructed graph through networkx. """

import networkx as nx
import matplotlib.pyplot as plt
import csv
import glob

from networkx.drawing.nx_pydot import graphviz_layout
from typing import List
from .pedigree import Node
from os import path, remove

# Configurations.

OCCUPIED_COLOR = 'cyan'
FREE_COLOR = 'green'
LABEL_SIZE = 5
NODE_SIZE = 100

def parse_data(bios_csv: str, degrees_csv: str):
    """
        Args:
            bios_csv: path to csv file of biographies for each individual.
            degrees_csv: path to csv file for pairwise relations.
        
        Parses through CSV files to construct the given nodes and a mapping
        of pairwise relations.
    """
    node_list = []
    with open(bios_csv, 'r') as bios_file:
        reader = csv.reader(bios_file)

        # Assume we know the ordering for now, so skip the first line.
        # TODO: should change later.
        next(reader)

        for row in reader:
            female = row[1] == 'F'
            node_list.append(Node(
                row[0], female, row[2], row[3] if not female else None, occupied=True, age=row[4]
            ))
    # Dictionary with key as degree and the value as a List of Tuples.
    pairwise_relations = {}
    with open(degrees_csv, 'r') as degrees_file:
        reader = csv.reader(degrees_file)

        #TODO: same issue as above.
        next(reader)
        for row in reader:
            cur = pairwise_relations.get(int(row[2]), [])
            cur.append((row[0], row[1]))
            pairwise_relations.update({int(row[2]) : cur})
    
    print(pairwise_relations)
    return node_list, pairwise_relations

def _format_label(node):
    return f'ID: {node.id}\n' \
           f'MtDna: {node.mt_dna}\n' \
           f'YChrom: {node.y_chrom}'

def visualize_graph(nodes: List[Node], filename):
    """
        Constructs a networkx graph from the given nodes
        and then displays it through matplotlib.
    """
    plt.figure(3)
    G = nx.DiGraph()
    
    ids = [node.id for node in nodes]
    females = [node for node in nodes if node.female]
    males = [node for node in nodes if not node.female]

    colormap_m = [OCCUPIED_COLOR if node.occupied else FREE_COLOR for node in males]
    colormap_f = [OCCUPIED_COLOR if node.occupied else FREE_COLOR for node in females]

    info = {node.id: _format_label(node) for node in nodes}

    G.add_nodes_from(ids)

    for node in nodes:
        # Only draw connections for children.
        for child in node.children:
            G.add_edge(node.id, child.id)

    pos = graphviz_layout(G, prog='dot')
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=list(map(lambda node: node.id, males)), node_shape='s',
        node_color=colormap_m,
        node_size=NODE_SIZE
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=list(map(lambda node: node.id, females)), node_shape='o',
        node_color=colormap_f,
        node_size=NODE_SIZE
    )

    nx.draw_networkx_edges(G, pos)
    nx.draw_networkx_labels(G, pos, labels=info, font_size=LABEL_SIZE)

    #nx.draw(G, pos, node_color=color_map, with_labels=True)
    plt.axis('off')
    
    cur = path.dirname(__file__)
    plt.savefig(path.join(cur, f'../output/{filename}'), dpi=300)
    plt.clf()
