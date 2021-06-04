#!/usr/bin/env python3

from __future__ import annotations
from .pedigree import Node
from typing import List
from copy import deepcopy
from typing import Tuple, List, Dict, Set, Optional
from contextlib import contextmanager

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

    def get_node (self, id : int) -> Node:
        """
            Returns the specified node by id, if not present, then 
            returns None.
        """
        return self.node_mapping.get(id)
    
    def size(self)->int:
        """
            Returns the length of nodes.
        """
        return len(self.node_list)
    
    def extrapolate_all(self) -> None:
        """
            Extrapolates all nodes.
        """
        for node in self.node_list:
            node.extrapolate()

    def update_nodes(self) -> None:
        """
            Updates all the nodes via breadth first search. Updates
            all auxiliary data structures as well.
        """
        visited = set()
        copy_list = [node for node in self.node_list if node.occupied]
        self.node_list = []
        self.node_mapping = {}

        for node in copy_list:
            self.node_mapping.update({node.id : node})

        def visit_edges(relations: List[Node], copy_list: List[Node]):
            for relative in relations:
                if relative not in visited:
                    copy_list.append(relative)

        # BFS search to get all the nodes in the visited set.
        while len(copy_list) > 0:
            node = copy_list.pop()
            if node in visited:
                continue
            visited.add(node)
            self.node_list.append(node)
            self.node_mapping.update({node.id : node})

            # Sufficient to visit only parents and children.
            visit_edges(node.parents, copy_list)
            visit_edges(node.children, copy_list)

        self.node_set = self.node_mapping.keys()
    
    def validate_nodes(self):
        """
            Validates nodes, only validates starting from one generation.
        """
        for node in self.node_list:
            if not node.female and node.y_chrom is None:
                for child in node.children:
                    if not child.female:
                        node.y_chrom = child.y_chrom
            if node.mt_dna is None and len(node.parents) != 0:
                node.mt_dna = node.parents[0].mt_dna
            

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

        for node in node_list:
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

@contextmanager
def _assign_parental (child: Node, parent: Node) -> None:
    """
        Helper function for assigning parental relationships between
        `child` and `parent`. Returns the Node tuple of the original
        parents of `child`. Returned value is used to backtrack.
    """
    assert(child is not None)
    assert(parent is not None)
    assert(child.parents is not None)
    assert(parent.children is not None)

    orig_mother = child.parents[0]
    orig_father = child.parents[1]

    if child.parents[0] == parent.parents[0] and \
       child.parents[1] == parent.parents[1]:
        yield False
        return
    if parent is orig_mother or parent is orig_father:
        yield True
        return

    if child is parent.parents[0] or child is parent.parents[1]:
        yield False
        return

    if (orig_mother.occupied and parent.female) or \
       (orig_father.occupied and not parent.female):
       yield False
       return

    # Begin assignment.
    to_replace = orig_father if not parent.female else orig_mother
    orig_parent_children = [node for node in parent.children]

    # Cycle detection first.
    for child in to_replace.children:
        if child.search_descendants([parent]):
            yield False
            return

    for child in to_replace.children:
        child.parents = (parent, child.parents[1]) if parent.female else (child.parents[0], parent)
        parent.children.append(child)
    
    orig_mt = child.mt_dna
    orig_ychrom = child.y_chrom

    if child.mt_dna is None and parent.female:
        child.mt_dna = parent.mt_dna
    if not child.female and child.y_chrom is None and not parent.female:
        child.y_chrom = parent.y_chrom

    yield True

    child.mt_dna = orig_mt
    child.y_chrom = orig_ychrom

    parent.children = orig_parent_children
    for child in to_replace.children:
        child.parents = (to_replace, child.parents[1]) if to_replace.female else (child.parents[0], to_replace)

@contextmanager
def _assign_sibling (sib1: Node, sib2: Node) -> None:
    """
        Assigns `sib1` and `sib2` as siblings.
    """
    sib1_parents = sib1.parents
    sib2_parents = sib2.parents

    mother = None
    father = None

    assert(sib1_parents is not None)
    assert(sib2_parents is not None)

    if sib1.search_descendants([sib2]) or sib2.search_descendants([sib1]):
        yield False
        return
    # Confirming existing relationship
    if sib1_parents[0] == sib2_parents[0] and \
    sib1_parents[1] == sib2_parents[1]:
        yield True
        return

    all_parents = sib1_parents + sib2_parents
    all_parents = list({node for node in all_parents if node.occupied})

    if len(all_parents) > 2:
        # More than two unique occupied nodes, impossible to merge.
        yield False
        return

    if len(all_parents) == 2:
        # Exactly two occupied parents.
        if all_parents[0].female == all_parents[1].female:
            # Same gender, wrong configuration.
            yield False
            return

        father = all_parents[0] if all_parents[1].female else all_parents[1]
        mother = all_parents[0] if all_parents[0].female else all_parents[1]

    elif len(all_parents) == 1:
        if all_parents[0].female:
            mother = all_parents[0]
            # Should use the father of the male child if there is one.
            father = sib1_parents[1] if not sib1.female else sib2_parents[1]
        else:
            father = all_parents[0]
            # Using any mother should be ok.
            mother = sib1_parents[0]
    else:
        # No occupied parents, anything goes as long as we reserve one father
        # from the male sibling.
        father = sib1_parents[1] if not sib1.female else sib2_parents[1]
        mother = sib1_parents[0]
    
    orig_mother_children = [child for child in mother.children]
    orig_father_children = [child for child in father.children]

    father_to_delete = sib1_parents[1] if sib1_parents[1] is not father else sib2_parents[1]
    mother_to_delete = sib1_parents[0] if sib1_parents[0] is not mother else sib2_parents[0]

    to_d_father_children = [child for child in father_to_delete.children]
    to_d_mother_children = [child for child in mother_to_delete.children]

    # Check for cycles first.
    if father_to_delete is not father:
        if father_to_delete.search_descendants([father, mother]):
            yield False
            return
        for child in father_to_delete.children:
            if child.search_descendants([father]):
                yield False
                return

    if mother_to_delete is not mother:
        if mother_to_delete.search_descendants([father, mother]):
            yield False
            return
        for child in mother_to_delete.children:
            if child.search_descendants([mother]):
                yield False
                return

    sibs = [sib1, sib2]
    for sib in sibs:
        if not sib.female and sib.y_chrom is not None:
            if sib.y_chrom != father.y_chrom:
                yield False
                return

    if father_to_delete is not father:
        for child in father_to_delete.children:
            father.children.append(child)
            child.parents = (child.parents[0], father)

    if mother_to_delete is not mother:
        for child in mother_to_delete.children:
            mother.children.append(child)
            child.parents = (mother, child.parents[1])
 
    sib1_orig_mt = sib1.mt_dna
    sib2_orig_mt = sib2.mt_dna

    sib1_orig_ychrom = sib1.y_chrom
    sib2_orig_ychrom = sib2.y_chrom

    to_assign_mt = sib1_orig_mt if sib1_orig_mt is not None else sib2_orig_mt
    sib1.mt_dna = to_assign_mt
    sib2.mt_dna = to_assign_mt

    if not sib1.female and not sib2.female:
        to_assign_ychrom = sib1_orig_ychrom if sib1_orig_ychrom is not None else sib2_orig_ychrom
        sib1.y_chrom = to_assign_ychrom
        sib2.y_chrom = to_assign_ychrom
    
    father_orig_ychrom = father.y_chrom
    if father.y_chrom is None:
        to_assign = sib1.y_chrom if not sib1.female else sib2.y_chrom
        father.y_chrom = to_assign
    yield True

    sib1.mt_dna = sib1_orig_mt
    sib2.mt_dna = sib2_orig_mt

    sib1.y_chrom = sib1_orig_ychrom
    sib2.y_chrom = sib2_orig_ychrom

    father.y_chrom = father_orig_ychrom
    father.children = orig_father_children
    mother.children = orig_mother_children

    father_to_delete.children = to_d_father_children
    mother_to_delete.children = to_d_mother_children

    for child in to_d_father_children:
        child.parents = (child.parents[0], father_to_delete)
    for child in to_d_mother_children:
        child.parents = (mother_to_delete, child.parents[1])


def _assign_helper(
        relation: List[Tuple[str, str]],
        graph: Graph,
        all_possible: List[Graph],
        idx: int,
        degree: int
    ) -> None:
    """
        Assigns according to all possible relationships specified by the
        `relation` list. All possible graphs are put into the `all_possible` list.
    """
    if relation is None or idx == len(relation):
        all_possible.append(deepcopy(graph))
        return

    src, dest = relation[idx]

    src = graph.get_node(src)
    dest = graph.get_node(dest)

    share_mt_dna = src.mt_dna == dest.mt_dna
    one_is_none = src.mt_dna is None or dest.mt_dna is None

    if degree == 2:
        second_rel = src.get_first_degree_rel()
        if dest in second_rel:
             _assign_helper(relation, graph, all_possible, idx + 1, degree)
             return

    # ------ ASSIGNMENTS ------

    # Case that source and dest are both male.
    if not src.female and not dest.female:

        share_y = src.y_chrom == dest.y_chrom
        one_chrom_none = src.y_chrom is None or dest.y_chrom is None

        if (share_y or one_chrom_none) and (share_mt_dna or one_is_none):
            # Must be siblings.
            with _assign_sibling(src, dest) as ok:
                if ok:
                    _assign_helper(relation, graph, all_possible, idx + 1, degree)
        if (share_y or one_chrom_none) and (one_is_none or not share_mt_dna):
            # Either father/son or son/father.
            # Father/son first.
            with _assign_parental(dest, src) as ok:
                if ok:
                    _assign_helper(relation, graph, all_possible, idx + 1, degree)
            # Son/father.
            with _assign_parental(src, dest) as ok:
                if ok:
                    _assign_helper(relation, graph, all_possible, idx + 1, degree)
        else:
            # No configuration works here.
            return

    # Case that one is female and the other is male.
    elif (not src.female and dest.female) or \
         (src.female and not dest.female):

            male_node = src if dest.female else dest
            female_node = src if src.female else dest
            if share_mt_dna or one_is_none:

                # Either siblings or son/mother.

                # Case 1 siblings.
                with _assign_sibling(male_node, female_node) as ok:
                    if ok:
                        _assign_helper(relation, graph, all_possible, idx + 1, degree)

                # Case 2 parental.
                with _assign_parental(male_node, female_node) as ok:
                    if ok:
                        _assign_helper(relation, graph, all_possible, idx + 1, degree)
                
            if not share_mt_dna:
                # Don't share mtDNA. Must be father/daughter.
                with _assign_parental(female_node, male_node) as ok:
                    if ok:
                        _assign_helper(relation, graph, all_possible, idx + 1, degree)

    # Case that source and dest are both females.
    else:
        if share_mt_dna or one_is_none:
            # May be siblings or daughter/mother or mother/daughter.
            
            # Case 1 siblings.
            with _assign_sibling(src, dest) as ok:
                if ok:
                    _assign_helper(relation, graph, all_possible, idx + 1, degree)
            
            # Case 2 daughter/mother.
            with _assign_parental(src, dest) as ok:
                if ok:
                    _assign_helper(relation, graph, all_possible, idx + 1, degree)

            # Case 3 mother/daugther.
            with _assign_parental(dest, src) as ok:
                if ok:
                    _assign_helper(relation, graph, all_possible, idx + 1, degree)
        else:
            # No configuration works here.
            return

# -------------
# GRAPH PRUNERS
#--------------


# def _prune_graphs3(
#     third_degrees: List[Tuple[str, str]],
#     node_map: Dict[str, Node],
#     occupied_nodes: List[Node],
#     all_possible: List[List[Node]]
# ) -> List[List[Node]]:
#     """
#         Prunes graphs that assign third degree relationships to nodes
#         that were not described in original pairwise relationships.
#     """

#     if third_degrees is None:
#         return all_possible

#     mapping = {}
#     for node in occupied_nodes:
#         lst = []
#         mapping.update({node.id : lst})
    
#     for rel in third_degrees:
#         lst = mapping.get(rel[0])
#         lst.append(rel[1])
#         mapping.update({rel[0] : lst})

#         lst1 = mapping.get(rel[1])
#         lst1.append(rel[0])
#         mapping.update({rel[1] : lst1})
#     ret = []

#     def _check_graph(graph: List[Node]) -> bool:
#         for node in graph:
#             if node.is_given():
#                 first_relatives = set(node.get_first_degree_rel())
#                 second_relatives = set(node.get_second_degree_rel())

#                 for rel in second_relatives:
#                     layer_second_relatives = set(rel.get_first_degree_rel())
#                     for third_rel in layer_second_relatives:
#                         if third_rel.is_given() and third_rel != node and third_rel not in first_relatives \
#                            and third_rel not in second_relatives and third_rel not in rel.parents:
#                            if third_rel.id not in mapping.get(node.id):
#                                return False
#         return True

#     # Begin pruning graphs.
#     for graph in all_possible:
#         if _check_graph(graph):
#             ret.append(graph)

#     return ret


def _prune_graphs2(
    second_degrees: List[Tuple[str, str]],
    graph: Graph,
    all_possible: List[Graph]
) -> List[List[Node]]:
    """
        Prunes graphs that assign third degree relationships to nodes
        that were not described in original pairwise relationships.
    """
    if second_degrees is None:
        return all_possible

    mapping = {}
    for node in graph.node_list:
        lst = []
        mapping.update({node.id : lst})
    
    for rel in second_degrees:
        lst = mapping.get(rel[0])
        lst.append(rel[1])
        mapping.update({rel[0] : lst})

        lst1 = mapping.get(rel[1])
        lst1.append(rel[0])
        mapping.update({rel[1] : lst1})
    ret = []

    def _check_graph(graph: Graph) -> bool:
        for node in graph.node_list:
            if node.is_given():
                first_relatives = set(node.get_first_degree_rel())
                for rel in first_relatives:
                    layer_first_relatives = set(rel.get_first_degree_rel())
                    for second_rel in layer_first_relatives:
                        if second_rel.is_given() and second_rel != node and second_rel not in first_relatives \
                           and second_rel not in rel.parents:
                            if second_rel.id not in mapping.get(node.id):
                                return False
        return True

    # Begin pruning graphs.
    for graph in all_possible:
        if _check_graph(graph):
            ret.append(graph)
    for graph in ret:
        for node in graph.node_list:
            if node.y_chrom is None:
                for child in node.children:
                    if not child.female:
                        node.y_chrom = child.y_chrom

    return ret

def _prune_graphs(
    first_degrees: List[Tuple[str, str]],
    graph: Graph,
    all_possible: List[Graph]
) -> List[Graph]:
    """
        Prunes graphs that assign first degree relationships to nodes
        that were not described in original pairwise relationships.
    """
    if first_degrees is None:
        return all_possible

    # Sort the first degrees.
    mapping = {}
    for node in graph.node_list:
        lst = []
        mapping.update({node.id : lst})
    
    for rel in first_degrees:
        lst = mapping.get(rel[0])
        lst.append(rel[1])
        mapping.update({rel[0] : lst})

        lst1 = mapping.get(rel[1])
        lst1.append(rel[0])
        mapping.update({rel[1] : lst1})
    ret = []

    def _check_graph(graph: Graph) -> bool:
        for node in graph.node_list:
            if node.occupied:
                first_relatives = node.get_first_degree_rel()
                for rel in first_relatives:
                    first, second = (rel.id, node.id) if rel.id < node.id else (node.id, rel.id)
                    if rel.occupied and second not in mapping.get(first):
                        return False
        return True

    # Begin pruning graphs.
    for graph in all_possible:
        if _check_graph(graph):
            ret.append(graph)

    for graph in ret:
        graph.validate_nodes()
    return ret

# ----------------
# MARK/EXTRAPOLATE
#-----------------
def _mark_and_extrapolate(graphs: List[Graph], extrap: bool) -> List[Node]:
    """ Function for marking the unoccupied nodes and then extrapolating them. """
    ret = []
    for graph in graphs:
        for node in graph.node_list:
            if not node.occupied:
                node.occupied = True
                if extrap:
                    node.extrapolate()
        graph.update_nodes()
        ret.append(graph)
    return ret

# -----------------
# DEGREE RELAXATION
#------------------
def _relax_helper(
    buffer: List[List[Tuple[str, str]]], 
    idx: int,
    temp: List[Tuple[str, str]],
    results: List[List[Tuple[str, str]]]
) -> None:

    """
        Recursive helper for generating all possible pairwise combinations.
    """

    if idx == len(buffer):
        results.append(deepcopy(temp))
        return
    
    current = buffer[idx]
    for rel in current:
        temp.append(rel)
        _relax_helper(buffer, idx + 1, temp, results)
        temp.pop()

def _relax_helper2 (
    buffer,
    idx: int,
    temp,
    results
) -> None:
    """
        Recursive helper for generating all pairwise relation
        dictionaries. Results stored with dictionaries
    """
    if idx == len(buffer):
        results.append(deepcopy(temp))
        return
    current = buffer[idx]
    for possibility in current:
        temp.update({idx + 1 : possibility})
        _relax_helper2(buffer, idx + 1, temp, results)
        temp.pop(idx + 1, None)

def _reduce_relation (first: Node, second: Node) -> List[Tuple[str, str]]:
    """
        Reduces relationship of first and second node by one degree. Returns
        all possible pairwise arrangements.
    """
    ret = []
    first_rel = first.get_first_degree_rel()
    for node in first_rel:
        if node.id is second.id:
            continue
        assert(node.id is not second.id)
        ret.append((second.id, node.id))
    second_rel = second.get_first_degree_rel()
    for node in second_rel:
        if node.id is first.id:
            continue
        assert(node.id is not first.id)
        ret.append((first.id, node.id))
    return ret

def _relax_degree(
        graph: List[Node],
        pairwise_relations: Dict[int, List[Tuple[str, str]]]
    ) -> List[List[Tuple[str, str]]]:
    """
        Decrements degrees by one, assigns a new relationship as well
        based on possible configurations.
    """
    known = {}
    for node in graph.node_list:
        assert(node.id not in known.keys())
        known.update({node.id: node})

    pairwise_relations.pop(1, None)

    possibilities = []
    for degree in pairwise_relations.keys():
        degree_possibilities = []
        buffer = []
        for rel in pairwise_relations.get(degree):
            first, second = known.get(rel[0]), known.get(rel[1])
            relaxed = _reduce_relation(first, second)
            buffer.append(relaxed)
        _relax_helper(buffer, 0, [], degree_possibilities)
        possibilities.append(degree_possibilities)

    ret = []
    _relax_helper2(possibilities, 0, {}, ret)
    return ret

# ------------------
# GRAPH CONSTRUCTION
#-------------------
def construct_all_graphs(
        current: Graph,
        pairwise_relations: Dict[int, List[Tuple[str, str]]],
        results: List[Graph],
        original_pairwise,
        degree: int,
        max: int
    ) -> None:
    """
        Constructs all graph from the given information of known
        nodes and their pairwise relationships with each other.
        If no pairwise relation exists between two nodes, we assume
        that the two nodes are not related.
    """
    assert (degree >= 1 and max <= 4)

    if degree == max:
        if current is not None and current.size() != 0:
            results.append(current)
        return

    # Edge case for extrapolation (first round)
    if degree == 1:
        current.extrapolate_all()
    

    # Pipeline: assign => prune => mark => relax.
    valid = []
    _assign_helper(pairwise_relations.get(1), current, valid, 0, degree)

    if degree == 1:
        valid = _prune_graphs(original_pairwise.get(1), current, valid)
    elif degree == 2:
        valid = _prune_graphs2(original_pairwise.get(2), current, valid)
    elif degree == 3:
        # valid = _prune_graphs3(original_pairwise.get(3), current, valid)
        pass

    # Don't extrapolate if we've hit the end.
    valid = _mark_and_extrapolate(valid, degree + 1 != max)
    i = 0
    for graph in valid:
        print(f'graph: {i}', degree)
        i += 1
        pairwise_copy = deepcopy(pairwise_relations)
        dicts = _relax_degree(graph, pairwise_copy)
        if degree == max - 1 or dicts is None or len(dicts) == 0:
            pairwise_map = deepcopy(pairwise_relations)
            pairwise_map.pop(1, None)
            construct_all_graphs(graph, pairwise_map, results, original_pairwise, degree + 1, max)
            continue
        for dict_pairs in dicts:
            pairwise_map = deepcopy(dict_pairs)
            construct_all_graphs(graph, pairwise_map, results, original_pairwise, degree + 1, max)

    return