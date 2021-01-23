#!/usr/bin/env python3

""" Runs pedigree construction on sample data. """
import csv

from constructor.pedigree import Node
from constructor.visualizer import visualize_graph
from typing import List
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
                row[0], female, row[2], row[3] if not female else None, occupied=True, age=row[4]
            ))
    # a = node_list[0]
    # b = node_list[1]
    # c = node_list[2]
    # d = node_list[3]
    # e = node_list[4]
    # f = node_list[5]


    # a.children.append(b)
    # c.children.append(b)
    # d.children.append(a)
    # e.children.append(a)
    # f.children.append(b)
    for node in node_list:
        node.extrapolate_node()
    node_list = visit_nodes(node_list)
    print(len(node_list))
    visualize_graph(node_list)

def visit_nodes(node_list: List[Node]) -> List[Node]:
    """
        Returns a complete list of nodes.
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

if __name__ == '__main__':
    main()

