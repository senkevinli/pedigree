#!/usr/bin/env python3

""" Runs pedigree construction on sample data. """
import glob

from constructor.pedigree import construct_graph, _visit_nodes
from constructor.util import parse_data, visualize_graph, post_process
from os import path, remove
from copy import deepcopy

DEGREES = 'third-degrees'
INPUT_BIO = f'{DEGREES}/simple_great_grand_bio.csv'
INPUT_DEGREE = F'{DEGREES}/simple_great_grand_degrees.csv'

def main():

    # Getting file paths.
    dirname = path.dirname(__file__)
    bios_csv = path.join(dirname, f'./samples/{INPUT_BIO}')
    degrees_csv = path.join(dirname, f'./samples/{INPUT_DEGREE}')
    
    node_list, mappings = parse_data(bios_csv, degrees_csv)
    
    results = []
    construct_graph(node_list, mappings, results, deepcopy(mappings), 1)
    real_results = []
    for result in results:
        real_results.append(post_process(_visit_nodes(result)))

    # Clean up leftover .png files.
    files = glob.glob(path.join(dirname, f'./output/*.png'))
    for f in files:
        remove(f)

    print(f'Generating: {len(real_results)} graphs')
    for i, node_list in enumerate(real_results):
        visualize_graph(node_list, f'graph{i}')
    

if __name__ == '__main__':
    main()

