#!/usr/bin/env python3

""" Runs pedigree construction on sample data. """
import csv

from constructor.pedigree import construct_graph
from constructor.util import parse_data, visualize_graph
from typing import List
from os import path

OUTPUT_NAME = 'ex.png'

def main():

    # Getting file paths.
    dirname = path.dirname(__file__)
    bios_csv = path.join(dirname, './samples/sample_bio.csv')
    degrees_csv = path.join(dirname, './samples/sample_degrees.csv')
    
    node_list, mappings = parse_data(bios_csv, degrees_csv)
    
    complete_nodes = construct_graph(node_list, mappings)
    visualize_graph(complete_nodes)

if __name__ == '__main__':
    main()

