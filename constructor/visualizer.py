#!/usr/bin/env python3

""" For visualizing the constructed graph through networkx. """

import networkx as nx
import matplotlib.pyplot as plt
import pydot

from networkx.drawing.nx_pydot import graphviz_layout
from typing import List
from .pedigree import Node

OCCUPIED_COLOR = 'cyan'
FREE_COLOR = 'green'

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
    plt.show()