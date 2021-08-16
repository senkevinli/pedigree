#!/usr/bin/env python3

from __future__ import annotations
from .pedigree import Node, Graph
from .util import visualize_graph_graphviz
from os import sys
from os import path
from typing import List
from copy import deepcopy
from typing import Tuple, List, Dict, Set, Optional
from contextlib import contextmanager

global_graph = 0
dirname = path.dirname(__file__)

# 1. prioritize rows where both samples are from the same cluster
# 2. prioritize for which only one sample is in the largest cluster (sort on largest cluster length in pair)
# 3. sort on # of times particular sample appears in 2nd or 3rd degree.
@contextmanager
def _assign_parental (child: Node, parent: Node, half_ok: Optional[bool] = False) -> None:
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

    # if (orig_mother.original and parent.female) or \
    #    (orig_father.original and not parent.female):
    #    yield False
    #    return

    if (orig_mother.occupied and parent.female) or \
       (orig_father.occupied and not parent.female):
       yield False
       return
    # Begin assignment.
    to_replace = orig_father if not parent.female else orig_mother
    orig_parent_children = [node for node in parent.children]

    if parent in child.get_nodes_in_cluster():
        yield False
        return

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
    
    if father_to_delete is not father and father.modified_search_ancestors(father_to_delete.children):
        yield False
        return

    if mother_to_delete is not mother and mother.modified_search_ancestors(mother_to_delete.children):
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


def reverse_ordering(rel_assignment):
    """
        Reverses ordering for assignments.
    """
    return lambda *args : rel_assignment(*args[::-1])

# def _assign_avuncular(src: Node, avunc: Node) -> None:
#     if src is None or len(src.parents) == 0:
#         return False
    
#     mom = src.parents[0]
#     dad = src.parents[1]

#     one_is_none = mom.mtDna is None or avunc.mtDna is None
#     if one_is_none or mom.mtDna == avunc.mtDna:
#         with _assign_sibling(mom, avunc) as result:
#             yield result

#     one_is_none = dad.mtDna is None or avunc.mtDna is None
#     if avunc.female:
#         if one_is_none or dad.mtDna == avunc.mtDna:
#             with _assign_sibling(dad, avunc) as result:
#                 yield result
#     else:
#         mt_ok = one_is_none or dad.mtDna == avunc.mt_dna
#         y_ok = (dad.y_chrom is None or avunc.y_chrom is None) or dad.y_chrom == avunc.y_chrom

#         if mt_ok and y_ok:
#             with _assign_sibling(dad, avunc) as result:
#                 yield result
#     return False

# def _assign_grandparents(src: Node, grandparent: Node) -> None:
#     if src is None or len(src.parents) == 0:
#         return False

#     mom = src.parents[0]
#     dad = src.parents[1]

#     if grandparent.female:
#         one_is_none = mom.mtDna is None or grandparent.mtDna is None
#         if one_is_none or mom.mtDna == grandparent.mtDna:
#             with _assign_parental(mom, grandparent, True) as result:
#                 yield result
#         if one_is_none or dad.mtDna == grandparent.mtDna:
#             with _assign_parental(dad, grandparent, True) as result:
#                 yield result
#     else:
#         with _assign_parental(mom, grandparent, True) as result:
#             yield result
        
#         one_is_none = dad.y_chrom is None or grandparent.y_chrom is None
#         if one_is_none or dad.y_chrom == grandparent.y_chrom:
#             with _assign_parental(dad, grandparent, True):
#                 yield result
#     return False

# def _assign_grandchildren(src: Node, grandchild: Node) -> None:
#     if src is None or len(src.parents) == 0:
#         return False
    
#     for child in src.children:
        
#     return False
            

def _assign_verification(
        graph: Graph,
        original: Dict[int, List[Tuple[str, str]]]
    ) -> bool:
    
    for node in graph.node_list:
        if node.is_given():
            first_rels = node.get_first_degree_rel()
            second_rels = node.get_second_degree_rel()
            first_rels = [rel for rel in first_rels if rel.is_given()]
            second_rels = [rel for rel in second_rels if rel.is_given()]

            for relative in first_rels:
                if (relative.id, node.id) not in original[1] and (node.id, relative.id) not in original[1]:
                    return False
            for relative in second_rels:
                if (relative.id, node.id) not in original[2] and (node.id, relative.id) not in original[2]:
                    return False
    return True

def _assign_helper(
        relation: List[Tuple[str, str]],
        graph: Graph,
        all_possible: List[Graph],
        idx: int,
        degree: int,
        original
    ) -> None:
    """
        Assigns according to all possible relationships specified by the
        `relation` list. All possible graphs are put into the `all_possible` list.
    """
    if relation is None or idx == len(relation):
        if not _assign_verification(graph, original):
            return
        all_possible.append(deepcopy(graph))
        return

    src, dest = relation[idx]

    src = graph.get_node(src)
    dest = graph.get_node(dest)

    share_mt_dna = src.mt_dna == dest.mt_dna
    one_is_none = src.mt_dna is None or dest.mt_dna is None

    if degree >= 2:
        second_rel = src.get_first_degree_rel()
        if dest in second_rel:
             _assign_helper(relation, graph, all_possible, idx + 1, degree, original)
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
                    _assign_helper(relation, graph, all_possible, idx + 1, degree, original)
        if (share_y or one_chrom_none) and (one_is_none or not share_mt_dna):
            # Either father/son or son/father.
            # Father/son first.
            with _assign_parental(dest, src) as ok:
                if ok:
                    _assign_helper(relation, graph, all_possible, idx + 1, degree, original)
            # Son/father.
            with _assign_parental(src, dest) as ok:
                if ok:
                    _assign_helper(relation, graph, all_possible, idx + 1, degree, original)
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
                        _assign_helper(relation, graph, all_possible, idx + 1, degree, original)

                # Case 2 parental.
                with _assign_parental(male_node, female_node) as ok:
                    if ok:
                        _assign_helper(relation, graph, all_possible, idx + 1, degree, original)
                
            if not share_mt_dna:
                # Don't share mtDNA. Must be father/daughter.
                with _assign_parental(female_node, male_node) as ok:
                    if ok:
                        _assign_helper(relation, graph, all_possible, idx + 1, degree, original)

    # Case that source and dest are both females.
    else:
        if share_mt_dna or one_is_none:
            # May be siblings or daughter/mother or mother/daughter.
            
            # Case 1 siblings.
            with _assign_sibling(src, dest) as ok:
                if ok:
                    _assign_helper(relation, graph, all_possible, idx + 1, degree, original)
            
            # Case 2 daughter/mother.
            with _assign_parental(src, dest) as ok:
                if ok:
                    _assign_helper(relation, graph, all_possible, idx + 1, degree, original)

            # Case 3 mother/daugther.
            with _assign_parental(dest, src) as ok:
                if ok:
                    _assign_helper(relation, graph, all_possible, idx + 1, degree, original)
        else:
            # No configuration works here.
            return

global_graph = 0

def _assign_helper_evolved(
        relation: List[List[Tuple[str, str]]],
        graph: Graph,
        all_possible: List[Graph],
        idx: int,
        degree: int,
        original: Dict[int, List[Tuple[str, str]]]
    ) -> None:
    """
        Assigns according to all possible relationships specified by the
        `relation` list. All possible graphs are put into the `all_possible` list.
    """
    global global_graph
    if relation is None or idx == len(relation):
        if not _assign_verification(graph, original):
            return
        # print('got one!')
        # visualize_graph_graphviz(graph, filename)
        # sys.exit(0)
        print(len(all_possible))
        all_possible.append(deepcopy(graph))
        return
    print("ASSIGNING THIS INDEX", idx, original[2][idx])
    if (idx == 1):
        global_graph = global_graph + 1
        filename = path.join(dirname, f'emergency/graph{global_graph}')
        visualize_graph_graphviz(graph, filename)
    choices = relation[idx]

    for current in choices:
        src, dest = current
        src = graph.get_node(src)
        dest = graph.get_node(dest)

        share_mt_dna = src.mt_dna == dest.mt_dna
        one_is_none = src.mt_dna is None or dest.mt_dna is None

        if degree >= 2:
            second_rel = src.get_first_degree_rel()
            if dest in second_rel:
                _assign_helper_evolved(relation, graph, all_possible, idx + 1, degree, original)
                continue

        # ------ ASSIGNMENTS ------

        # Case that source and dest are both male.
        if not src.female and not dest.female:

            share_y = src.y_chrom == dest.y_chrom
            one_chrom_none = src.y_chrom is None or dest.y_chrom is None

            if (share_y or one_chrom_none) and (share_mt_dna or one_is_none):
                # Must be siblings.
                with _assign_sibling(src, dest) as ok:
                    original2 = deepcopy(original)

                    if ok and _remove_thirds(graph, original2):
                        _assign_helper_evolved(relation, graph, all_possible, idx + 1, degree, original2)
            if (share_y or one_chrom_none) and (one_is_none or not share_mt_dna):
                # Either father/son or son/father.
                # Father/son first.
                with _assign_parental(dest, src) as ok:
                    original2 = deepcopy(original)

                    if ok and _remove_thirds(graph, original2):
                        _assign_helper_evolved(relation, graph, all_possible, idx + 1, degree, original2)
                # Son/father.
                with _assign_parental(src, dest) as ok:
                    original2 = deepcopy(original)

                    if ok and _remove_thirds(graph, original2):
                        _assign_helper_evolved(relation, graph, all_possible, idx + 1, degree, original2)
            else:
                # No configuration works here.
                continue

        # Case that one is female and the other is male.
        elif (not src.female and dest.female) or \
            (src.female and not dest.female):

                male_node = src if dest.female else dest
                female_node = src if src.female else dest
                if share_mt_dna or one_is_none:

                    # Either siblings or son/mother.

                    # Case 1 siblings.
                    with _assign_sibling(male_node, female_node) as ok:
                        original2 = deepcopy(original)

                        if ok and _remove_thirds(graph, original2):
                            _assign_helper_evolved(relation, graph, all_possible, idx + 1, degree, original2)

                    # Case 2 parental.
                    with _assign_parental(male_node, female_node) as ok:
                        original2 = deepcopy(original)

                        if ok and _remove_thirds(graph, original2):
                            _assign_helper_evolved(relation, graph, all_possible, idx + 1, degree, original2)
                    
                if not share_mt_dna:
                    # Don't share mtDNA. Must be father/daughter.
                    with _assign_parental(female_node, male_node) as ok:
                        original2 = deepcopy(original)

                        if ok and _remove_thirds(graph, original2):
                            _assign_helper_evolved(relation, graph, all_possible, idx + 1, degree, original2)

        # Case that source and dest are both females.
        else:
            if share_mt_dna or one_is_none:
                # May be siblings or daughter/mother or mother/daughter.
                
                # Case 1 siblings.
                with _assign_sibling(src, dest) as ok:
                    original2 = deepcopy(original)

                    if ok and _remove_thirds(graph, original2):
                        _assign_helper_evolved(relation, graph, all_possible, idx + 1, degree, original2)
                
                # Case 2 daughter/mother.
                with _assign_parental(src, dest) as ok:
                    original2 = deepcopy(original)

                    if ok and _remove_thirds(graph, original2):
                        _assign_helper_evolved(relation, graph, all_possible, idx + 1, degree, original2)

                # Case 3 mother/daugther.
                with _assign_parental(dest, src) as ok:
                    original2 = deepcopy(original)

                    if ok and _remove_thirds(graph, original2):
                        _assign_helper_evolved(relation, graph, all_possible, idx + 1, degree, original2)
            else:
                # No configuration works here.
                continue


def _check_ok(
        relation: Tuple[str, str],
        graph: Graph
    ) -> bool:
    """
        Assigns according to all possible relationships specified by the
        `relation` list. All possible graphs are put into the `all_possible` list.
    """

    src, dest = relation

    src = graph.get_node(src)
    dest = graph.get_node(dest)

    share_mt_dna = src.mt_dna == dest.mt_dna
    one_is_none = src.mt_dna is None or dest.mt_dna is None

    # ------ CHECKS ------

    # Case that source and dest are both male.
    if not src.female and not dest.female:

        share_y = src.y_chrom == dest.y_chrom
        one_chrom_none = src.y_chrom is None or dest.y_chrom is None

        if (share_y or one_chrom_none) and (share_mt_dna or one_is_none):
            # Must be siblings.
            with _assign_sibling(src, dest) as ok:
                if ok:
                   return True
        if (share_y or one_chrom_none) and (one_is_none or not share_mt_dna):
            # Either father/son or son/father.
            # Father/son first.
            with _assign_parental(dest, src) as ok:
                if ok:
                    return True
            # Son/father.
            with _assign_parental(src, dest) as ok:
                if ok:
                    return True
        else:
            # No configuration works here.
            return False

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
                        return True

                # Case 2 parental.
                with _assign_parental(male_node, female_node) as ok:
                    if ok:
                       return True
                
            if not share_mt_dna:
                # Don't share mtDNA. Must be father/daughter.
                with _assign_parental(female_node, male_node) as ok:
                    if ok:
                        return True
            return False

    # Case that source and dest are both females.
    else:
        if share_mt_dna or one_is_none:
            # May be siblings or daughter/mother or mother/daughter.
            
            # Case 1 siblings.
            with _assign_sibling(src, dest) as ok:
                if ok:
                    return True
            
            # Case 2 daughter/mother.
            with _assign_parental(src, dest) as ok:
                if ok:
                   return True

            # Case 3 mother/daugther.
            with _assign_parental(dest, src) as ok:
                if ok:
                    return True
        else:
            # No configuration works here.
            return False

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
        print(len(results))
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

def _reduce_relation (first: Node, second: Node, graph: Graph) -> List[Tuple[str, str]]:
    """
        Reduces relationship of first and second node by one degree. Returns
        all possible pairwise arrangements.
    """

    ret = []
    first_rel = set(first.get_first_degree_rel())
    second_rel = set(second.get_first_degree_rel())
    

    for node in first_rel.difference(second_rel):
        if node.id is second.id or (node.original and second.original):
            continue
        if not _check_ok((node.id, second.id), graph):
            continue
        assert(node.id is not second.id)
        ret.append((second.id, node.id))
    
    for node in second_rel.difference(first_rel):
        if node.id is first.id or ((node.original and first.original)):
            continue
        if not _check_ok((node.id, second.id), graph):
            continue
        assert(node.id is not first.id)
        ret.append((first.id, node.id))
    
    for node in first_rel.intersection(second_rel):
        if node.original and first.original:
            continue
        if not _check_ok((node.id, second.id), graph):
            continue
        ret.append((first.id, node.id))
    return ret

def _relax_degree(
        graph: Graph,
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
        print(degree)
        degree_possibilities = []
        buffer = []
        for rel in pairwise_relations.get(degree):
            first, second = known.get(rel[0]), known.get(rel[1])
            relaxed = _reduce_relation(first, second, graph)
            buffer.append(relaxed)
        # from pprint import pprint
        # iterations = 1
        # for elem in buffer:
        #     iterations *= len(elem)
        #     pprint(elem)
        #     print(len(elem))
        # print(iterations)
        possibilities = buffer
        break
        # _relax_helper(buffer, 0, [], degree_possibilities)
        # possibilities.append(degree_possibilities)

    # ret = []
    # _relax_helper2(possibilities, 0, {}, ret)
    from pprint import pprint
    pprint(possibilities)
    ret = possibilities
    return ret

# ------------------
# GRAPH CONSTRUCTION
#-------------------

def prob_construct_helper(
        current: Graph,
        id_pairs: List[Tuple[str, str]],
        probs: List[List[int]],
        relationship_arr: List[function],
        current_prob: float,
        prob_results: List[float],
        graph_results: List[Graph],
        idx: int,
        original
    ):

    if idx == len(id_pairs):
        if not _assign_verification(current, original):
            return
        prob_results.append(current_prob)
        graph_results.append(deepcopy(current))
        return
    
    id1, id2 = id_pairs[idx]
    relation_probs = probs[idx]

    for assigner, prob in zip(relationship_arr, relation_probs):
        if prob == 0:
            continue
        first_node = current.get_node(id1)
        second_node = current.get_node(id2)

        with assigner(first_node, second_node) as ok:
            if ok:
                prob_construct_helper(current, id_pairs, probs, relationship_arr,
                                      current_prob * prob, prob_results, graph_results, idx + 1, original)



def construct_all_known(
        current: Graph,
        first_probs: Dict[Tuple[str, str], List[int]],
        results: List[Graph],
        graph_probabilities: List[float],
        original
    ) -> None:
    """
        Constructs all graphs from the given first degree relationships.
    """
    reverse_parental = reverse_ordering(_assign_parental)
    relationship_arr = [
        reverse_parental, # father-daughter
        _assign_parental, # daughter-father
        reverse_parental, # father-son
        _assign_parental, # son-father
        _assign_parental, # son-mother
        reverse_parental, # mother-son
        _assign_parental, # daughter-mother
        reverse_parental, # mother-daughter
        _assign_sibling   # siblings
    ]

    result_graphs = results
    result_probs = graph_probabilities

    id_pairs = []
    probs = []

    for key, val in first_probs.items():
        id_pairs.append(key)
        probs.append(val)
    
    current.extrapolate_all()
    prob_construct_helper(current, id_pairs, probs, relationship_arr, 1, result_probs,
                          result_graphs, 0, original)
    
    return result_graphs, result_probs

def construct_all_graphs(
        current: Graph,
        pairwise_relations: Dict[int, List[Tuple[str, str]]],
        results: List[Graph],
        original_pairwise,
        degree: int,
        max: int,
        first_probs = None
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
    if not first_probs:
        if degree == 1:
            _assign_helper(pairwise_relations.get(1), current, valid, 0, degree, original_pairwise)
        else:
            _assign_helper_evolved(pairwise_relations, current, valid, 0, degree, original_pairwise)
    else:
        probs = []
        construct_all_known(current, first_probs, valid, probs, original_pairwise)
        print(f'Probabilities for degree: {degree} are: {probs}')

    if degree == 1:
        pass
        valid = _prune_graphs(original_pairwise.get(1), current, valid)
    elif degree == 2:
        pass
        valid = _prune_graphs2(original_pairwise.get(2), current, valid)
    elif degree == 3:
        # valid = _prune_graphs3(original_pairwise.get(3), current, valid)
        pass

    # Don't extrapolate if we've hit the end.
    valid = _mark_and_extrapolate(valid, degree + 1 != max)

    i = 0
    for graph in valid:
        i += 1
        pairwise_copy = deepcopy(pairwise_relations)
        if degree == 1:
            val = _remove_seconds(graph, pairwise_copy)
            # print(pairwise_copy)
            if not val:
                continue
        
        # val = _remove_thirds(graph, pairwise_copy)
        # if not val:
        #     continue

        if degree != max - 1:
            dicts = _relax_degree(graph, pairwise_copy)
        if degree == max - 1 or dicts is None or len(dicts) == 0:
            pairwise_map = deepcopy(pairwise_relations)
            # pairwise_map.pop(1, None)
            pairwise_map = {}
            construct_all_graphs(graph, pairwise_map, results, original_pairwise, degree + 1, max)
            continue
        # for dict_pairs in dicts:
        #     pairwise_map = deepcopy(dict_pairs)
        #     construct_all_graphs(graph, pairwise_map, results, original_pairwise, degree + 1, max)
        construct_all_graphs(graph, dicts, results, original_pairwise, degree + 1, max)


    return

def _remove_thirds(graph: Graph, pairs: dict):
    if not pairs or 3 not in pairs.keys():
        return True
    # print("paris is ", pairs)
    def _find_cluster(src: str, clusters):
        for cluster in clusters:
            if src in cluster:
                return cluster
        return None
    clusters = []

    for node in graph.node_list:
        second_rel = node.get_second_degree_rel()
        first_rel = node.get_first_degree_rel()
        less_two_rel = second_rel + first_rel
        flag = False
        for cluster in clusters:
            if set(set([node.id for node in less_two_rel]).intersection(set(cluster))):
                cluster.append(node.id)
                flag = True
                break
        if not flag:
            clusters.append([node.id])
    
    third_rels = pairs.get(3)
    new_rels = []
    for rel in third_rels:
        src, dest = rel
        src, dest = graph.get_node(src), graph.get_node(dest)

        src_cluster = _find_cluster(src.id, clusters)
        if src_cluster is None:
            return False
        
        # print(src, [node.id for node in src.get_third_degree_rel()])
        if dest.id in src_cluster and dest not in src.get_third_degree_rel():
            # visualize_graph_graphviz(graph, filename)
            print(dest.id, src.id)
            return False
        if dest.id not in src_cluster:
            new_rels.append((src.id, dest.id))
    
    pairs[3] = new_rels
    return True

def _remove_seconds(graph: Graph, pairs: dict):
    if 2 not in pairs.keys():
        return True
    def _find_cluster(src: str, clusters):
        for cluster in clusters:
            if src in cluster:
                return cluster
        return None

    clusters = []
    for node in graph.node_list:
        first_rel = node.get_first_degree_rel()
        flag = False
        for cluster in clusters:
            if set([node.id for node in first_rel]).intersection(set(cluster)):
                cluster.append(node.id)
                flag = True
                break
        if not flag:
            clusters.append([node.id])

    second_rels = pairs.get(2)
    new_rels = []
    for rel in second_rels:
        src, dest = rel
        src, dest = graph.get_node(src), graph.get_node(dest)

        src_cluster = _find_cluster(src.id, clusters)
        if src_cluster is None:
            return False
        
        if dest.id in src_cluster and dest not in src.get_second_degree_rel():
            return False
        if dest.id not in src_cluster:
            new_rels.append((src.id, dest.id))
    
    # print("new rels: ", new_rels)
    pairs[2] = new_rels
    return True
                


