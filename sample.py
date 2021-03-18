#!/usr/bin/env python3

""" Runs pedigree construction on sample data. """
import glob

from constructor.pedigree import construct_graph
from constructor.util import parse_data, visualize_graph
from os import path, remove

DEGREES = 'second-degrees'
INPUT_BIO = f'{DEGREES}/grandparents_bio.csv'
INPUT_DEGREE = F'{DEGREES}/grandparents_degrees.csv'

def main():

    # Getting file paths.
    dirname = path.dirname(__file__)
    bios_csv = path.join(dirname, f'./samples/{INPUT_BIO}')
    degrees_csv = path.join(dirname, f'./samples/{INPUT_DEGREE}')
    
    node_list, mappings = parse_data(bios_csv, degrees_csv)
    
    results = []
    construct_graph(node_list, mappings, results)

    # Clean up leftover .png files.
    files = glob.glob(path.join(dirname, f'./output/*.png'))
    for f in files:
        remove(f)

    for i, node_list in enumerate(results):
        visualize_graph(node_list, f'graph{i}')
    

if __name__ == '__main__':
    main()

