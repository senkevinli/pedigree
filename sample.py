#!/usr/bin/env python3

""" Runs pedigree construction on sample data. """
import glob

from constructor.graph import construct_all_graphs, construct_all_known
from constructor.util import parse_data, compare_isomorph, visualize_graph_graphviz, parse_prob
from os import path, remove
from copy import deepcopy

DEGREES = 'real'
FILE = 'hazelton'
MAX = 3
PROB = True

INPUT_BIO = f'{DEGREES}/{FILE}_bio.csv'
INPUT_DEGREE = f'{DEGREES}/{FILE}_degrees.csv'

INPUT_PROB = None
if PROB:
    INPUT_PROB = f'{DEGREES}/{FILE}_prob.csv'

OUTPUT_DIR = './output3' 

def main():

    # Getting file paths.
    dirname = path.dirname(__file__)
    bios_csv = path.join(dirname, f'./samples/{INPUT_BIO}')
    degrees_csv = path.join(dirname, f'./samples/{INPUT_DEGREE}')
    prob_csv = None

    if PROB:
        prob_csv = path.join(dirname, f'./samples/{INPUT_PROB}')
    
    node_list, mappings = parse_data(bios_csv, degrees_csv)
    probabilities = parse_prob(prob_csv, node_list)
    
    results = []

    construct_all_graphs(node_list, mappings, results, deepcopy(mappings), 1, MAX, probabilities)

    # Clean up leftover .png files.
    files = glob.glob(path.join(dirname, f'{OUTPUT_DIR}/*'))
    for f in files:
        remove(f)
    
    results = compare_isomorph(results)
    print(f'Generating: {len(results)} graphs')
    for i, node_list in enumerate(results):
        visualize_graph_graphviz(node_list, path.join(dirname, f'{OUTPUT_DIR}/graph{i}'))

if __name__ == '__main__':
    main()

