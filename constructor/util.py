#!/usr/bin/env python3

""" For visualizing the constructed graph through networkx. """

import networkx as nx
import matplotlib.pyplot as plt
import csv

from networkx.drawing.nx_pydot import graphviz_layout
from typing import List
from .pedigree import Node
from os import path

OCCUPIED_COLOR = 'cyan'
FREE_COLOR = 'green'

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
            cur = pairwise_relations.get(row[2], [])
            cur.append((row[0], row[1]))
            pairwise_relations.update({row[2] : cur})
    
    print(pairwise_relations)
    return node_list, pairwise_relations


def visualize_graph(nodes: List[Node]):
    """
        Constructs a networkx graph from the given nodes
        and then displays it through matplotlib.
    """
    G = nx.DiGraph()
    
    ids = [node.id for node in nodes]
    color_map = [OCCUPIED_COLOR if node.occupied else FREE_COLOR for node in nodes]
    G.add_nodes_from(ids, node_shape='o')

    for node in nodes:
        # Only draw connections for children.
        for child in node.children:
            G.add_edge(node.id, child.id)

    pos = graphviz_layout(G, prog='dot')
    nx.draw(G, pos, node_color=color_map, with_labels=True)
    plt.axis('off')
    
    cur = path.dirname(__file__)
    plt.savefig(path.join(cur, '../output/ex.png'))

def visit_nodes(node_list: List[Node]) -> List[Node]:
    """
        Returns a complete list of nodes. Found through
        a BFS search of children and parents.
    """
    visited = set()

    def visit_edges(relations: List[Node]):
        for relative in relations:
            if relative not in visited:
                visited.add(relative)
                node_list.append(relative)

    # BFS search to get all the nodes in the visited set.
    while len(node_list) > 0:
        node = node_list.pop()
        visited.add(node)

        # Sufficient to visit only parents and children.
        visit_edges(node.parents)
        visit_edges(node.children)
    return list(visited)

