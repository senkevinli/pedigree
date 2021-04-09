#!/usr/bin/env python3

""" For visualizing the constructed graph through networkx. """

import networkx as nx
import matplotlib.pyplot as plt
import csv
import glob
import re

from networkx.drawing.nx_pydot import graphviz_layout
from typing import List
from .pedigree import Node
from os import path, remove

# Configurations.

OCCUPIED_COLOR = 'cyan'
FREE_COLOR = 'green'
KNOWN_COLOR = 'yellow'

LABEL_SIZE = 5
NODE_SIZE = 100

def post_process(graph: List[Node]) -> List[Node]:
    """
        Post process graph to remove extraneous nodes.
    """
    ret = []
    for node in graph:
        #  if re.search('[0-9]+$', node.id):
        #      rels = node.get_first_degree_rel()
        #      if len(rels) == 1:
        #          rels[0].parents = ()
        #          continue
         ret.append(node)
    return ret
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
                row[0], female, row[2], row[3] if not female else None, occupied=True, original=True, age=row[4]
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
    
    colormap_f = []
    for node in females:
        if re.search('^[^0-9]+$', node.id):
            colormap_f.append(KNOWN_COLOR)
        elif node.occupied:
            colormap_f.append(OCCUPIED_COLOR)
        else:
            colormap_f.append(FREE_COLOR)

    colormap_m = []
    for node in males:
        if re.search('^[^0-9]+$', node.id):
            colormap_m.append(KNOWN_COLOR)
        elif node.occupied:
            colormap_m.append(OCCUPIED_COLOR)
        else:
            colormap_m.append(FREE_COLOR)

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
