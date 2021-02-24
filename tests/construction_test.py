#!/usr/bin/env python3

""" Tests for construction of graph only. """

from os import path
from constructor.pedigree import construct_graph
from constructor.util import parse_data

def test_simple():
    dirname = path.dirname(__file__)
    simple_bio = path.join(dirname, './fixtures/simple_bio.csv')
    simple_degrees = path.join(dirname, './fixtures/simple_degrees.csv')

    node_list, mappings = parse_data(simple_bio, simple_degrees)
    complete_nodes = construct_graph(node_list, mappings)