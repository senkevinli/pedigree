#!/usr/bin/env python3

""" For traversing the pedigree and other related operations. """

from __future__ import annotations
from copy import deepcopy
from typing import Tuple, List, Dict, Optional

AGE_DEFAULT = 100
class Node:
    filler_id = 0
    def __init__(
        self,
        id: str,
        female: bool,
        mt_dna: Optional[str] = None,
        y_chrom: Optional[str] = None,
        occupied: Optional[bool] = False,
        age: Optional[int] = AGE_DEFAULT,
        parents: Optional[Tuple[Node, Node]] = None,
        children: Optional[List[Node]] = None,
        partners: Optional[List[Node]] = None,
        siblings: Optional[List[Node]] = None
    ) -> None:
        """
            Constructs a node representing an individual in the pedigree.
            Args:
                `female`: represents gender of the individual, true is female, false is male.
                `mt_dna`: represents the type of mitochondrial dna.
                `y_chrom`: represents the type of y chromosome -- only if male.
                `age`: age of the individual, if not given then assumes individual is an adult.
                `parents`: biological parent nodes, both must be of the opposite sex.
                `children`: list of children that this node fathers/mothers.
                `partners`: node qualifies as a partner if there exists a child where this node and
                            the node in question are parents of the child.
                `siblings`: full siblings of the current node.
        """
        self.female = female
        self.mt_dna = mt_dna
        self.age = age
        self.occupied = occupied

        # Tuple, first element is mother, second is father.
        if parents is not None:
            assert parents[0].female
            assert not parents[1].female

        self.parents = parents if parents is not None else ()
        self.children = children if children is not None else []

        if female:
            assert y_chrom is None

        self.y_chrom = y_chrom
        self.siblings = siblings if siblings is not None else []
        self.partners = partners if partners is not None else []
        self.id = id
    
    def __str__(self):
        """
            String representation for debugging.
        """
        return f'id: {self.id}\n' + \
               f'gender: {"F" if self.female else "M"}\n' + \
               f'mtDna: {self.mt_dna}\n' + \
               f'yChrom: {self.y_chrom}\n' + \
               f'age: {self.age}\n'

    def extrapolate(self):
        """
            Constructs surrounding nodes of a father, mother, two siblings of each gender,
            a partner, and two childrens of each gender. Extrapolate for all generations
            until the last one.
        """
        # TODO: Make this clean.
        # Add parents, don't add if already present.
        if len(self.parents) == 0:
            Node.filler_id += 1
            father = Node(
                str(Node.filler_id),
                False,
                y_chrom=self.y_chrom if not self.female else None,
                occupied=False,
                children=[self]
            )
            Node.filler_id += 1
            mother = Node(
                str(Node.filler_id),
                True,
                mt_dna=self.mt_dna,
                occupied=False,
                children=[self]
            )
            self.parents = (mother, father)

        assert len(self.parents) == 2

        # Give two full siblings of each gender.
        Node.filler_id += 1
        fem_sibling = Node(
            str(Node.filler_id),
            True,
            mt_dna=self.mt_dna,
            parents=self.parents,   
        )
        Node.filler_id += 1
        male_sibling = Node(
            str(Node.filler_id),
            False,
            mt_dna=self.mt_dna,
            y_chrom=self.parents[1].y_chrom,
            parents=self.parents
        )

        fem_sibling.siblings += [male_sibling, self]
        male_sibling.siblings += [male_sibling, self]

        for parent in self.parents:
            parent.children += [male_sibling, fem_sibling]

        # Make filler partner.
        Node.filler_id += 1
        partner = Node(
            str(Node.filler_id),
            not self.female,
            partners=[self]
        )
        self.partners.append(partner)

        # Give two children of each gender.
        Node.filler_id += 1
        fem_child = Node(
            str(Node.filler_id),
            True,
            mt_dna=self.mt_dna if self.female else None,
            parents=(self, partner) if self.female else (partner, self) 
        )
        Node.filler_id += 1
        male_child = Node(
            str(Node.filler_id),
            False,
            mt_dna=self.mt_dna if self.female else None,
            y_chrom=self.y_chrom if not self.female else None,
            parents=(self, partner) if self.female else (partner, self)
        )

        partner.children += [fem_child, male_child]
        self.children += [fem_child, male_child]

def _visit_nodes(node_list: List[Node]) -> List[Node]:
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
        visit_edges(node.partners)
    return list(visited)

def construct_graph(
        node_list: List[Node],
        pairwise_relations: Dict[int, List[Tuple[str, str]]]
    ) -> List[Node]:
    """
        Constructs a graph from the given information of known
        nodes and their pairwise relationships with each other.
        If no pairwise relation exists between two nodes, we assume
        that the two nodes are not related.
    """
    for node in node_list:
        node.extrapolate()
    ret = _visit_nodes(node_list)
    ret = deepcopy_graph(ret)
    return ret

def deepcopy_graph(node_list: List[Node]) -> List[Node]:
    """
        Constructs a deepcopy of the graph and returns it.
        Does not modify the graph that is passed in.
    """
    all_nodes = _visit_nodes(node_list)
    node_mapping = {}

    for node in all_nodes:
        copied = deepcopy(node)

        # All IDs in the list should be unique.
        assert(node.id not in node_mapping.keys())
        node_mapping.update({node.id: copied})
    
    for node in node_mapping.values():
        # Reset all connections to our copied nodes based on ID.
        node.children = [node_mapping.get(rel.id) for rel in node.children]
        node.partners = [node_mapping.get(rel.id) for rel in node.partners]
        node.siblings = [node_mapping.get(rel.id) for rel in node.siblings]
        if len(node.parents) > 0:
            node.parents = (
                node_mapping.get(node.parents[0].id),
                node_mapping.get(node.parents[1].id)
            )
    ret = [node for node in node_mapping.values()]
    return ret