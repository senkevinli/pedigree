#!/usr/bin/env python3

from .pedigree import Node
from typing import List
from copy import deepcopy

class Graph:
    def __init__(self, given_nodes: List[Node]) -> None:
        """
            Constructs a graph that is a deep copy of the given node list.
            3 data structures => node list, node dictionary, node set.
        """

        self.node_list = deepcopy(given_nodes)

        # Construct mapping.
        self.node_mapping = {}
        for node in self.node_list:
            if node.id not in self.node_mapping.keys():
                self.node_mapping.update({node.id : node})
        
        self.node_set = set(self.node_mapping.keys())
    
    def __str__(self) -> str:
        """
            Returns string representation (all the node ids).
        """
        ret = ''
        for node in self.node_list:
            ret += node.id
        return ret

    def update_nodes(self) -> None:
        """
            Updates all the nodes via breadth first search. Updates
            all auxiliary data structures as well.
        """
        visited = set()
        copy_list = [node for node in self.node_list if not node.occupied]
        self.node_list = []
        self.node_mapping = {}

        def visit_edges(relations: List[Node], copy_list: List[Node]):
            for relative in relations:
                if relative not in visited:
                    visited.add(relative)
                    self.node_list.add(relative)
                    self.node_mapping.update({node.id : node})
                    copy_list.append(relative)

        # BFS search to get all the nodes in the visited set.
        while len(copy_list) > 0:
            node = copy_list.pop()
            if node in visited:
                continue
            visited.add(node)
            self.node_list.add(node)
            self.node_mapping.update({node.id : node})
            # Sufficient to visit only parents and children.
            visit_edges(node.parents)
            visit_edges(node.children)
        
        self.node_set = self.node_mapping.keys()

    def __deepcopy__(self, memo):
        
        # Update before deepcopying.
        self.update_nodes()

        copy = type(self)([])
        memo[id(self)] = copy
        
        node_list = []
        node_mapping = {}

        # First create mapping with correctly copied nodes.
        for node in self.node_list:
            copied = deepcopy(node)
            assert(node.id not in node_mapping.keys())
            node_mapping.update({node.id : copied})
            node_list.append(copied)
        
        for node in node_mapping.values():
            node.children = [node_mapping.get(rel.id) for rel in node.children]
            if len(node.parents) > 0:
                node.parents = (
                    node_mapping.get(node.parents[0].id),
                    node_mapping.get(node.parents[1].id)
                )
    
        copy.node_list = node_list
        copy.node_mapping = node_mapping
        copy.node_set = set(copy.node_mapping.keys())

        return copy
