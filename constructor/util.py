#!/usr/bin/env python3

""" For visualizing the constructed graph through graphviz. """

import networkx as nx
import matplotlib.pyplot as plt
import csv
import glob
import re
import math

from networkx.drawing.nx_pydot import graphviz_layout
from distinctipy import distinctipy
from graphviz import Graph, Digraph, Source
from typing import List
from .pedigree import Node
from os import path, remove
from copy import deepcopy
from colormap import rgb2hex
from .graph import Graph

# Configurations.

LABEL_SIZE = 5
NODE_SIZE = 100


def _gender_top_sort(graph: Graph) -> List[Node]:
    """
        Takes in a graph as a list of nodes. Performs
        gender topological sort on the family tree.
    """

    # Preprocess
    graph = graph.node_list
    mapping = {}
    for node in graph:
        mapping.update({node.id : node})

    ids = [node.id for node in graph if node.children is None or len(node.children) == 0]
    ids.sort()

    num_mapping = {}
    for i, id in enumerate(ids):
        num_mapping.update({id : i})


    # Parallel array
    visited = {}
    for node in graph:
        visited.update({node.id : False})

    result = []

    def visit(node_id):
        node = mapping.get(node_id)
        if visited[node_id]:
            return
        if node.parents is None or len(node.parents) == 0:
            result.append(node_id)
            visited[node_id] = True
        else:
            # Visit mother first
            visit(node.parents[0].id)
            visit(node.parents[1].id)
            result.append(node_id)
            visited[node_id] = True

    for i, ident in enumerate(ids):
        if visited[ident]:
            continue
        visit(ident)
    
    real = []
    for res in result:
        real.append(mapping.get(res).female)
    return real


def compare_isomorph(graphs: List[List[Node]]) -> List[List[Node]]:
    """
        Prunes graphs that are similar/isomorphic.
    """
    ret = []

    tops = [_gender_top_sort(graph) for graph in graphs]
    ok = [True for top in tops]
    

    for i in range(len(tops)):
        if ok[i]:
            ret.append(graphs[i])
        else:
            continue
        for j in range(i + 1, len(tops)):
            if tops[i] == tops[j]:
                ok[j] = False
                continue
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
    ret = Graph(node_list)
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

    return ret, pairwise_relations

def parse_prob(prob_csv: str, graph: Graph):
    if prob_csv is None:
        return None
    with open(prob_csv, 'r') as prob_file:
        reader = csv.reader(prob_file)

        ret = {}
        next(reader)
        for row in reader:
            key = (row[0], row[1])
            val = [float(prob_val) for prob_val in row[2:]]
            ret.update({key: val})
    return ret


def visualize_graph_graphviz(nodes: Graph, name) -> None:
    """
        Visualizes graph using Graphvis instead of NetworkX. Different
        Nodes and their colors are represented in a key.
    """

    # Colors for mitochondrial
    nodes = nodes.node_list
    colors = distinctipy.get_colors(14, pastel_factor=1)
    mt_map = {}
    i = 0
    for node in nodes:
        if node.mt_dna not in mt_map:
            mt_map.update({node.mt_dna : colors[i]})
            i += 1

    # Colors for y chromosome
    colors = distinctipy.get_colors(14, pastel_factor=0.1)
    y_map = {}
    i = 0
    for node in nodes:
        if node.y_chrom not in y_map:
            y_map.update({node.y_chrom: colors[i]})
            i += 1

    # Constructing graph
    u = Digraph('Nodes', format='png')
    u.attr('node', shape='circle')
    for i, node in enumerate(nodes):
        color1 = mt_map.get(node.mt_dna)
        hex_code = rgb2hex(math.floor(color1[0]*255), math.floor(color1[1]*255), math.floor(color1[2]*255))

        u.attr('node', color=f'{hex_code}')
        if not node.female:
            u.attr('node', shape='square', penwidth='1')
            color2 = y_map.get(node.y_chrom)
            if not node.is_given():
                u.attr('node', style='dashed')
            else:
                u.attr('node', style='filled')
        else:
            u.attr('node', shape='circle', penwidth='1')
            if not node.is_given():
                u.attr('node', style='dashed')
            else:
                u.attr('node', style='filled')

        color2 = y_map.get(node.y_chrom)
        hex_code = rgb2hex(math.floor(color2[0]*255), math.floor(color2[1]*255), math.floor(color2[2]*255)) \
                   if not node.female else "#000000"
        u.node(node.id, label=f'<<font color="{hex_code}"> {node.id} </font>>')

    known = [node]
    for node in nodes:
        # Only draw connections for children.
        for child in node.children:
            u.edge(node.id, child.id)

    with u.subgraph(name='cluster_0') as c:
        c.attr(label=f'MtDNA Key')
        c.attr('node', shape='plaintext', color='white')
        mt_dna_table = '<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">'
        for mt, color in mt_map.items():
            hex_code = rgb2hex(math.floor(color[0]*255), math.floor(color[1]*255), math.floor(color[2]*255))
            mt_dna_table += f'<TR><TD> {mt if mt is not None else "Unknown"} </TD><TD BGCOLOR="{hex_code}"></TD></TR>'

        mt_dna_table += '</TABLE>>'
        c.node('table', label=mt_dna_table)

    with u.subgraph(name='cluster_1') as c:
        c.attr(label=f'yChrom Key')
        c.attr('node', shape='plaintext', color='white')
        y_table = '<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">'
        for y, color in y_map.items():
            hex_code = rgb2hex(math.floor(color[0]*255), math.floor(color[1]*255), math.floor(color[2]*255))
            y_table += f'<TR><TD> {y if y is not None else "Unknown"} </TD><TD BGCOLOR="{hex_code}"></TD></TR>'

        y_table += '</TABLE>>'
        c.node('table2', label=y_table)
    
    u.render(filename=name)
