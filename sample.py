#!/usr/bin/env python3

""" Runs pedigree construction on sample data. """
import glob

from constructor.graph import construct_all_graphs
# from constructor.pedigree import construct_graph, _visit_nodes
from constructor.util import parse_data, visualize_graph, compare_isomorph, visualize_graph_graphviz
from os import path, remove
from copy import deepcopy

DEGREES = 'first-degrees'
INPUT_BIO = f'{DEGREES}/disjoint_bio.csv'
INPUT_DEGREE = F'{DEGREES}/disjoint_degrees.csv'

def main():

    # Getting file paths.
    dirname = path.dirname(__file__)
    bios_csv = path.join(dirname, f'./samples/{INPUT_BIO}')
    degrees_csv = path.join(dirname, f'./samples/{INPUT_DEGREE}')
    
    node_list, mappings = parse_data(bios_csv, degrees_csv)
    
    results = []
    construct_all_graphs(node_list, mappings, results, deepcopy(mappings), 1, 2)

    # Clean up leftover .png files.
    files = glob.glob(path.join(dirname, f'./output2/*'))
    for f in files:
        remove(f)
    
    results = compare_isomorph(results)
    print(f'Generating: {len(results)} graphs')
    for i, node_list in enumerate(results):
        # visualize_graph(node_list, path.join(dirname, f'./output2/graph{i}'))
        visualize_graph_graphviz(node_list, path.join(dirname, f'./output2/graph{i}'))

if __name__ == '__main__':
    main()

