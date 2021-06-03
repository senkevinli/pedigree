#!/usr/bin/env python3

""" For traversing the pedigree and other related operations. """

from __future__ import annotations
import re
from contextlib import contextmanager
#`from constructor.util import visualize_graph
from copy import deepcopy
from typing import Tuple, List, Dict, Set, Optional

AGE_DEFAULT = 100
DEGREE_CAP = 4
MAX = 3
class Node:
    filler_id = 0
    def __init__(
        self,
        id: str,
        female: bool,
        mt_dna: Optional[str] = None,
        y_chrom: Optional[str] = None,
        occupied: Optional[bool] = False,
        original: Optional[bool] = False,
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
        self.original = original

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

    def is_given(self):
        return self.original

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

    def search_descendants(self, nodes: List[Node]) -> bool:
        """
            Searches children for the specified `node`. Returns
            true if found, false otherwise. Requires that there
            must not be a cycle in the graph.
        """
        if self.children is None or len(self.children) == 0:
            return False
        for child in self.children:
            if child in nodes:
                return True
            if child.search_descendants(nodes):
                return True
        return False
    
    def get_first_degree_rel(self) -> List[Node]:
        """
            Returns all first degree relatives of the current node.
        """
        ret = []

        # Parents.
        if (self.parents != None and len(self.parents) == 2):
            ret.append(self.parents[0])
            ret.append(self.parents[1])

        # Children.
        ret += self.children

        # Siblings.
        if (self.parents != None and len(self.parents) == 2):
            for child in self.parents[0].children:
                if child.parents[1] == self.parents[1] and child != self:
                    ret.append(child)

        return ret
    
    def get_second_degree_rel(self) -> List[Node]:
        """
            Returns all second degree relatives of the current node.
        """

        ret = []

        first = set(self.get_first_degree_rel())
        for node in first:
            layer = set(node.get_first_degree_rel())
            for second in layer:
                if second != self and second not in first \
                   and second not in node.parents:
                   ret.append(second)
        return ret
            
    
    def search_entire_tree(self, other: Node, visited):
        if self is other:
            return True
        for child in self.children:
            if child in visited:
                continue
            visited.add(child)
            val = child.search_entire_tree(other, visited)
            if val:
                return True

        if self.parents != None and len(self.parents) == 2:
            if self.parents[0] not in visited:
                visited.add(self.parents[0])
                val = self.parents[0].search_entire_tree(other, visited)
                if val:
                    return True
            if self.parents[1] not in visited:
                visited.add(self.parents[1])
                val2 = self.parents[1].search_entire_tree(other, visited)
                if val2:
                    return True
        return False

def _visit_nodes(node_list: List[Node]) -> List[Node]:
    """
        Returns a complete list of nodes.
    """
    temp = []
    for node in node_list:
        if node.occupied:
            temp.append(node)
    if (len(temp) != len(node_list)):
        return _visit_nodes(temp)
    visited = set()
    copy_list = [node for node in node_list]

    def visit_edges(relations: List[Node]):
        for relative in relations:
            if relative not in visited:
                visited.add(relative)
                copy_list.append(relative)

    # BFS search to get all the nodes in the visited set.
    while len(copy_list) > 0:
        node = copy_list.pop()
        visited.add(node)

        # Sufficient to visit only parents and children.
        visit_edges(node.parents)
        visit_edges(node.children)
    return list(visited)


def _assign_helper(
        relation: List[Tuple[str, str]],
        node_map: Dict[str, Node],
        node_list: List[Node],
        all_possible: List[List[Node]],
        idx: int
    ) -> None:
    if relation is None or idx == len(relation):
        all_possible.append(deepcopy_graph(node_list))
        return

    src, dest = relation[idx]
    # We want to be able to get back to our original state,
    # cache the original relationships for both destination
    # and source.
    src_node = deepcopy(node_map.get(src))
    dest_node = deepcopy(node_map.get(dest))

    src = node_map.get(src)
    dest = node_map.get(dest)

    share_mt_dna = src.mt_dna == dest.mt_dna
    one_is_none = src.mt_dna is None or dest.mt_dna is None

    # ------ ASSIGNMENTS ------

    # Case that source and dest are both male.
    if not src.female and not dest.female:

        share_y = src.y_chrom == dest.y_chrom
        one_chrom_none = src.y_chrom is None or dest.y_chrom is None

        if (share_y or one_chrom_none) and (share_mt_dna or one_is_none):
            # Must be siblings.
            with _assign_sibling(src, dest) as ok:
                if ok:
                    _assign_helper(relation, node_map, node_list, all_possible, idx + 1)
        if (share_y or one_chrom_none) and (one_is_none or not share_mt_dna):
            # Either father/son or son/father.
            # Father/son first.
            with _assign_parental(dest, src) as ok:
                if ok:
                    _assign_helper(relation, node_map, node_list, all_possible, idx + 1)
            # Son/father.
            with _assign_parental(src, dest) as ok:
                if ok:
                    _assign_helper(relation, node_map, node_list, all_possible, idx + 1)
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
                        _assign_helper(relation, node_map, node_list, all_possible, idx + 1)

                # Case 2 parental.
                with _assign_parental(male_node, female_node) as ok:
                    if ok:
                        _assign_helper(relation, node_map, node_list, all_possible, idx + 1)
                
            if not share_mt_dna:
                # Don't share mtDNA. Must be father/daughter.
                with _assign_parental(female_node, male_node) as ok:
                    if ok:
                        _assign_helper(relation, node_map, node_list, all_possible, idx + 1)

    # Case that source and dest are both females.
    else:
        if share_mt_dna or one_is_none:
            # May be siblings or daughter/mother or mother/daughter.
            
            # Case 1 siblings.
            with _assign_sibling(src, dest) as ok:
                if ok:
                    _assign_helper(relation, node_map, node_list, all_possible, idx + 1)
            
            # Case 2 daughter/mother.
            with _assign_parental(src, dest) as ok:
                if ok:
                    _assign_helper(relation, node_map, node_list, all_possible, idx + 1)

            # Case 3 mother/daugther.
            with _assign_parental(dest, src) as ok:
                if ok:
                    _assign_helper(relation, node_map, node_list, all_possible, idx + 1)
        else:
            # No configuration works here.
            return

# ------ ASSIGNMENT SUBHELPER METHODS ------

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
    
    # if sib1.search_entire_tree(sib2, set()):
    #     yield False
    #     return
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

def _validate_graph3(
    third_degrees: List[Tuple[str, str]],
    node_map: Dict[str, Node],
    occupied_nodes: List[Node],
    all_possible: List[List[Node]]
) -> List[List[Node]]:
    """
        Prunes graphs that assign third degree relationships to nodes
        that were not described in original pairwise relationships.
    """

    if third_degrees is None:
        return all_possible

    mapping = {}
    for node in occupied_nodes:
        lst = []
        mapping.update({node.id : lst})
    
    for rel in third_degrees:
        lst = mapping.get(rel[0])
        lst.append(rel[1])
        mapping.update({rel[0] : lst})

        lst1 = mapping.get(rel[1])
        lst1.append(rel[0])
        mapping.update({rel[1] : lst1})
    ret = []

    def _check_graph(graph: List[Node]) -> bool:
        for node in graph:
            if node.is_given():
                first_relatives = set(node.get_first_degree_rel())
                second_relatives = set(node.get_second_degree_rel())

                for rel in second_relatives:
                    layer_second_relatives = set(rel.get_first_degree_rel())
                    for third_rel in layer_second_relatives:
                        if third_rel.is_given() and third_rel != node and third_rel not in first_relatives \
                           and third_rel not in second_relatives and third_rel not in rel.parents:
                           if third_rel.id not in mapping.get(node.id):
                               return False
        return True

    # Begin pruning graphs.
    for graph in all_possible:
        if _check_graph(graph):
            ret.append(graph)

    return ret


def _prune_graphs3(
    third_degrees: List[Tuple[str, str]],
    node_map: Dict[str, Node],
    occupied_nodes: List[Node],
    all_possible: List[List[Node]]
) -> List[List[Node]]:
    """
        Prunes graphs that assign third degree relationships to nodes
        that were not described in original pairwise relationships.
    """

    if third_degrees is None:
        return all_possible

    mapping = {}
    for node in occupied_nodes:
        lst = []
        mapping.update({node.id : lst})
    
    for rel in third_degrees:
        lst = mapping.get(rel[0])
        lst.append(rel[1])
        mapping.update({rel[0] : lst})

        lst1 = mapping.get(rel[1])
        lst1.append(rel[0])
        mapping.update({rel[1] : lst1})
    ret = []

    def _check_graph(graph: List[Node]) -> bool:
        for node in graph:
            if node.is_given():
                first_relatives = set(node.get_first_degree_rel())
                second_relatives = set(node.get_second_degree_rel())

                for rel in second_relatives:
                    layer_second_relatives = set(rel.get_first_degree_rel())
                    for third_rel in layer_second_relatives:
                        if third_rel.is_given() and third_rel != node and third_rel not in first_relatives \
                           and third_rel not in second_relatives and third_rel not in rel.parents:
                           if third_rel.id not in mapping.get(node.id):
                               return False
        return True

    # Begin pruning graphs.
    for graph in all_possible:
        if _check_graph(graph):
            ret.append(graph)

    return ret


def _prune_graphs2(
    second_degrees: List[Tuple[str, str]],
    node_map: Dict[str, Node],
    occupied_nodes: List[Node],
    all_possible: List[List[Node]]
) -> List[List[Node]]:
    """
        Prunes graphs that assign third degree relationships to nodes
        that were not described in original pairwise relationships.
    """
    if second_degrees is None:
        return all_possible

    mapping = {}
    for node in occupied_nodes:
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

    def _check_graph(graph: List[Node]) -> bool:
        for node in graph:
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
        for node in graph:
            if node.y_chrom is None:
                for child in node.children:
                    if not child.female:
                        node.y_chrom = child.y_chrom

    return ret

def _prune_graphs(
    first_degrees: List[Tuple[str, str]],
    node_map: Dict[str, Node],
    occupied_nodes: List[Node],
    all_possible: List[List[Node]]
) -> List[List[Node]]:
    """
        Prunes graphs that assign first degree relationships to nodes
        that were not described in original pairwise relationships.
    """
    if first_degrees is None:
        return all_possible
    # Sort the first degrees.
    mapping = {}
    for node in occupied_nodes:
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
    # mapped = map(lambda x: sorted(x), first_degrees)

    # first_degree_map = {}
    # for rel in mapped:
    #     lst = first_degree_map.get(rel[0], [])
    #     lst.append(rel[1])
    #     first_degree_map.update({rel[0] : lst})

    def _check_graph(graph: List[Node]) -> bool:
        for node in graph:
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
        for node in graph:
            if node.y_chrom is None:
                for child in node.children:
                    if not child.female:
                        node.y_chrom = child.y_chrom
    return ret

def _mark_and_extrapolate(graphs: List[List[Node]], extrap: bool) -> List[Node]:
    """ Function for marking the unoccupied nodes and then extrapolating them. """
    ret = []
    for graph in graphs:
        for node in graph:
            if not node.occupied:
                node.occupied = True
                if extrap:
                    node.extrapolate()
        ret.append(_visit_nodes(graph))
    return ret

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

def _relax_degree(
        graph: List[Node],
        pairwise_relations: Dict[int, List[Tuple[str, str]]]
    ) -> List[List[Tuple[str, str]]]:
    """
        Decrements degrees by one, assigns a new relationship as well
        based on possible configurations.
    """
    known = {}
    for node in graph:
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


def construct_graph(
        node_list: List[Node],
        pairwise_relations: Dict[int, List[Tuple[str, str]]],
        results: List[List[Node]],
        original_pairwise,
        degree: int
    ) -> None:
    """
        Constructs a graph from the given information of known
        nodes and their pairwise relationships with each other.
        If no pairwise relation exists between two nodes, we assume
        that the two nodes are not related.
    """
    if degree == MAX:
        if node_list is not None and len(node_list) != 0:
            results.append(node_list)
        return

    # Construct mapping for known nodes.
    known = {}
    for node in node_list:
        assert(node.id not in known.keys())
        known.update({node.id: node})

    # Edge case for extrapolation (first round)
    if degree == 1:
        for node in node_list:
            node.extrapolate()

    # Pipeline: assign => prune => mark => relax.
    valid = []
    _assign_helper(pairwise_relations.get(1), known, node_list, valid, 0)

    if degree == 1:
        valid = _prune_graphs(original_pairwise.get(1), known, node_list, valid)
    elif degree == 2:
        valid = _prune_graphs2(original_pairwise.get(2), known, node_list, valid)
    elif degree == 3:
        valid = _prune_graphs3(original_pairwise.get(3), known, node_list, valid)

    # Don't extrapolate if we've hit the end.
    valid = _mark_and_extrapolate(valid, degree + 1 != MAX)
    i = 0
    for graph in valid:
        i += 1
        pairwise_copy = deepcopy(pairwise_relations)
        dicts = _relax_degree(graph, pairwise_copy)
        if degree == MAX - 1 or dicts is None or len(dicts) == 0:
            pairwise_map = deepcopy(pairwise_relations)
            pairwise_map.pop(1, None)
            construct_graph(graph, pairwise_map, results, original_pairwise, degree + 1)
            continue
        for dict_pairs in dicts:
            pairwise_map = deepcopy(dict_pairs)
            construct_graph(graph, pairwise_map, results, original_pairwise, degree + 1)

    return

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
