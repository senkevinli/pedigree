#!/usr/bin/env python3

""" Runs pedigree construction on sample data. """
import csv

from constructor.pedigree import construct_graph
from constructor.util import parse_data, visualize_graph
from typing import List
from os import path

INPUT_BIO = 'simple_bio.csv'
INPUT_DEGREE = 'simple_degrees.csv'

def main():

    # Getting file paths.
    dirname = path.dirname(__file__)
    bios_csv = path.join(dirname, f'./samples/{INPUT_BIO}')
    degrees_csv = path.join(dirname, f'./samples/{INPUT_DEGREE}')
    
    node_list, mappings = parse_data(bios_csv, degrees_csv)
    
    complete_nodes = construct_graph(node_list, mappings)
    for i, node_list in enumerate(complete_nodes):
        visualize_graph(node_list, f'graph{i}')

if __name__ == '__main__':
    main()

