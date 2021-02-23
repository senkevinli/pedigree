#!/usr/bin/env python3

""" For traversing the pedigree and other related operations. """

from __future__ import annotations
from contextlib import contextmanager
#`from constructor.util import visualize_graph
from copy import deepcopy
from typing import Tuple, List, Dict, Set, Optional

AGE_DEFAULT = 100
DEGREE_CAP = 4
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

def _visit_nodes(node_list: List[Node]) -> List[Node]:
    """
        Returns a complete list of nodes.
    """
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
    if idx == len(relation):
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

    # ------ ASSIGNMENTS ------

    # Case that source and dest are both male.
    if not src.female and not dest.female:
        share_y = src.y_chrom == dest.y_chrom
        if share_y and share_mt_dna:
            # Must be siblings.

            # First, add parental linkages. Change references 
            # only for those that are occupied.
            src_parents = src.parents
            dest_parents = dest.parents

            assert(src_parents is not None)
            assert(dest_parents is not None)

            # Confirming existing relationship
            if src_parents[0] == dest_parents[0] and \
               src_parents[1] == dest_parents[1]:
               _assign_helper(relation, node_map, node_list, all_possible, idx + 1)
               return
            else:
                total = src_parents + dest_parents
                occupied = [node for node in total if node.occupied]

                # Two or more occupied parents from different nodes, wrong configuration.
                if len(occupied) > 2:
                    return
                if len(occupied) == 2:
                    # Same gender, wrong configuration.
                    if occupied[0].female == occupied[1].female:
                        return
                    # Different genders, merge.
                    father = occupied[0] if not occupied[0].female else occupied[1]
                    mother = occupied[1] if occupied[0].female else occupied[1]

                    # Save state to revert.
                    orig_src_parents = src.parents
                    orig_dest_parents = dest.parents
                    orig_father_children = [child for child in father.children]
                    orig_mother_children = [child for child in mother.children]

                    src.parents = (mother, father)
                    dest.parents = (mother, father)

                    for child in mother.children:
                        if child not in father.children:
                            father.children.append(child)
                            child.parents = (mother, father)
                    for child in father.children:
                        if child not in mother.children:
                            mother.children.append(child)
                            child.parents = (mother, father)

                    _assign_helper(relation, node_map, node_list, all_possible, idx + 1)

                    src.parents = orig_src_parents
                    dest.parents = orig_dest_parents
                    for child in mother.children:
                        if child not in orig_mother_children:
                            child


        elif share_y:
            # Either father/son or son/father.

            # Father/son first.
            if dest.parents[1].occupied and \
               dest.parents[1] != src:
               return
            elif dest.parents[1] == src:
                _assign_helper(relation, node_map, node_list, all_possible, idx + 1)
            else:
                orig_dest_parents = dest.parents
                mother = dest.parents[0]
                orig_mother_children = [child for child in mother.children]
                orig_src_children = [child for child in src.children]
                dest.parents = (mother, src)
                src.children.append(dest)

                # for child in mother.children:
                #     if child not in src.children:
                #         child.parents = (mother, src)
                #         src.children.append(child)
                # for child in src.children:
                #     if child not in mother.children:
                #         child.parents = (mother, src)
                #         mother.children.append(child)
                _assign_helper(relation, node_map, node_list, all_possible, idx + 1)

                dest.parents = orig_dest_parents
                src.children.remove(dest)
                # for child in src.children:
                #     if child in orig_mother_children:
                #         child.parents = orig_dest_parents
                #         src.children.remove(child)

                # for child in mother.children:
                #     if child in orig_src_children:
                #         child.parents = (src.children[0].parents[0], src)
                #         mother.children.remove(child)

            pass
        else:
            # No configuration works here.
            return

    # Case that one is female and the other is male.
    elif (not src.female and dest.female) or \
         (src.female and not dest.female):

            male_node = src if dest.female else dest
            female_node = src if src.female else dest
            if share_mt_dna:

                # Either siblings or son/mother.
                # Case 1 siblings.

                male_node_parents = male_node.parents
                female_node_parents = female_node.parents
                mother = None
                father = None

                assert(male_node_parents is not None)
                assert(female_node_parents is not None)

                # Confirming existing relationship
                if male_node_parents[0] == female_node_parents[0] and \
                   male_node_parents[1] == female_node_parents[1]:
                   _assign_helper(relation, node_map, node_list, all_possible, idx + 1)
                   return

                total = male_node_parents + female_node_parents
                occupied = [node for node in total if node.occupied]
                occupied = set(occupied)
                occupied = list(occupied)
                if len(occupied) > 2:
                    return
                if len(occupied) == 2:
                    # Same gender, wrong configuration.
                    if occupied[0].female == occupied[1].female:
                        return
                    father = occupied[0] if occupied[1].female else occupied[1]
                    mother = occupied[0] if occupied[0].female else occupied[1]
                elif len(occupied) == 1:
                    if occupied[0].female:
                        mother = occupied[0]
                        father = male_node_parents[1]
                    else:
                        mother = male_node_parents[0]
                        father = occupied[0]
                else:
                    # Must use father, since this tells us the most information.
                    father = male_node_parents[1]
                    mother = female_node_parents[0]
                    print(mother)
                    print(father)
                orig_mother_children = [child for child in mother.children]
                orig_father_children = [child for child in father.children]
                orig_father_partners = [child.parents[0] for child in father.children]
                orig_mother_partners = [child.parents[1] for child in mother.children]

                for child in mother.children:
                    if child not in father.children:
                        child.parents = (mother, father)
                        father.children.append(child)
                for child in father.children:
                    if child not in mother.children:
                        child.parents = (mother, father)
                        mother.children.append(child)

                _assign_helper(relation, node_map, node_list, all_possible, idx + 1)

                father.children = orig_father_children
                mother.children = orig_mother_children

                for i, child in enumerate(mother.children):
                    child.parents = (mother, orig_mother_partners[i])
                for i, child in enumerate(father.children):
                    child.parents = (father, orig_father_partners[i])
                return

                        
            else:
                # Don't share mtDNA. Must be father/daughter.

                orig_female_node_parents = female_node.parents

                if female_node.parents[1].occupied and \
                   female_node.parents[1] != male_node:

                    # Configuration is impossible.
                    return 
                if female_node.parents[1] == male_node:

                    # Confirming existing relationship, continue.
                    _assign_helper(relation, node_map, node_list, all_possible, idx + 1)
                else:
                    
                    with _assign_parental(female_node, male_node):
                        _assign_helper(relation, node_map, node_list, all_possible, idx + 1)

    # Case that source and dest are both females.
    else:
        if share_mt_dna:
            # May be siblings or daughter/mother or mother/daughter.
            pass
        else:
            # No configuration works here.
            return

# ------ ASSIGNMENT SUBHELPER METHODS ------

@contextmanager
def _assign_parental (child: Node, parent: Node) -> Tuple[bool, Node]:
    """
        Helper function for assigning parental relationships between
        `child` and `parent`. Returns the Node tuple of the original
        parents of `child`. Returned value is used to backtrack.
    """

    try:
        assert(child is not None)
        assert(parent is not None)
        assert(child.parents is not None)
        assert(parent.children is not None)

        orig_mother = child.parents[0]
        orig_father = child.parents[1]

        orig_mother_children = orig_mother.children
        orig_father_children = orig_father.children

        # Begin assignment.
        parent.children.append(child)

        child.parents = (parent, orig_father) if parent.female else (orig_mother, parent)
        yield

    finally:
        child.parents = (orig_mother, orig_father)
        parent.children.remove(child)


    
 
    return (orig_mother, orig_father)

# def _assign_sibling (sib1: Node, sib2: Node) -> None:
#     """
#         Assigns `sib1` and `sib2` as siblings.
#     """

#     sib1_parents = sib1.parents
#     sib2_parents = sib2.parents

#     mother = None
#     father = None

#     assert(sib1_parents is not None)
#     assert(sib2_parents is not None)

#     # Confirming existing relationship
#     if sib1_parents[0] == sib2_parents[0] and \
#        sib1_parents[1] == sib2_parents[1]:
#         return (True, None)

#     all_parents = sib1_parents + sib2_parents
#     all_parents = list({node for node in all_parents if node.occupied})

#     if len(all_parents) > 2:
#         # More than two unique occupied nodes, impossible to merge.
#         return (False, None)

#     if len(all_parents) == 2:
#         # Exactly two occupied parents.
#         if all_parents[0].female == all_parents[1].female:
#             # Same gender, wrong configuration.
#             return (False, None)
#         father = all_parents[0] if all_parents[1].female else all_parents[1]
#         mother = all_parents[0] if all_parents[0].female else all_parents[1]
    
#     elif len(all_parents) == 1:
#         if all_parents[0].female:
#             mother = all_parents[0]
#             # Should use the father of the male child if there is one.
#             father = sib1_parents[1] if not sib1.female else sib2_parents[1]
#         else:
#             father = all_parents[0]
#             # Using any mother should be ok.
#             mother = sib1_parents[0]
#     else:
#         # No occupied parents, anything goes as long as we reserve one father
#         # from the male sibling.
#         father = sib1_parents[1] if not sib1.female else sib2_parents[1]
#         mother = sib1_parents[0]
    
#     orig_mother_children = [child for child in mother.children]
#     orig_father_children = [child for child in father.children]

#     orig_father_partners = [child.parents[0] for child in father.children]
#     orig_mother_partners = [child.parents[1] for child in mother.children]

#     f_children = set(father.children)
#     m_children = set(mother.children)

#     father_to_delete = sib1_parents[1] if sib1_parents[1] is not father else sib2_parents[1]
#     mother_to_delete = sib1_parents[0] if sib1_parents[0] is not mother else sib2_parents[0]

#     for child in 

    # occupied = set(occupied)
    # occupied = list(occupied)

    # if len(occupied) > 2:
    #     # More than two unique occupied nodes, impossible to merge.
    #     return (False, None)

    # if len(occupied) == 2:
        
    #     if occupied[0].female == occupied[1].female:
    #         # Same gender, wrong configuration.
    #         return (False, None)

    #     father = occupied[0] if occupied[1].female else occupied[1]
    #     mother = occupied[0] if occupied[0].female else occupied[1]

    # elif len(occupied) == 1:
    #     if occupied[0].female:
    #         mother = occupied[0]
    #         father = male_node_parents[1]
    #     else:
    #         mother = male_node_parents[0]
    #         father = occupied[0]
    # else:
    #     # Must use father, since this tells us the most information.
    #     father = male_node_parents[1]
    #     mother = female_node_parents[0]
    #     print(mother)
    #     print(father)

    # orig_mother_children = [child for child in mother.children]
    # orig_father_children = [child for child in father.children]

    # orig_father_partners = [child.parents[0] for child in father.children]
    # orig_mother_partners = [child.parents[1] for child in mother.children]

    # f_children = set(father.children)
    # m_children = set(mother.children)

    # for child in mother.children:
    #     if child not in f_children:
    #         child.parents = (mother, father)
    #         father.children.append(child)

    # for child in father.children:
    #     if child not in mother.children:
    #         child.parents = (mother, father)
    #         mother.children.append(child)

    # _assign_helper(relation, node_map, node_list, all_possible, idx + 1)

    # father.children = orig_father_children
    # mother.children = orig_mother_children

    # for i, child in enumerate(mother.children):
    #     child.parents = (mother, orig_mother_partners[i])
    # for i, child in enumerate(father.children):
    #     child.parents = (father, orig_father_partners[i])
    # return


def _construct_helper(
        relations: Dict[int, List[Tuple[str, str]]],
        node_map: Dict[str, Node],
        node_list: List[Node],
        all_possible: List[List[Node]]
    ) -> bool:
    """
        Recursive helper for assigning nodes.
    """
    if len(relations.keys() == 0):
        all_possible.append(deepcopy_graph(node_list))
        # Stop recursing.
        return
    
    # First, assign all the ones that are degree 1.
    to_assign = relations.get(1)

    for i, relation in enumerate(to_assign):
        pass


def construct_graph(
        node_list: List[Node],
        pairwise_relations: Dict[int, List[Tuple[str, str]]]
    ) -> List[List[Node]]:
    """
        Constructs a graph from the given information of known
        nodes and their pairwise relationships with each other.
        If no pairwise relation exists between two nodes, we assume
        that the two nodes are not related.
    """
    # Construct mapping for known nodes.
    known = {}
    for node in node_list:
        assert(node.id not in known.keys())
        known.update({node.id: node})

    # Extrapolate first
    for node in node_list:
        node.extrapolate()
        
        # Begin assigning.
    results = []
    #print(node_list)
    _assign_helper(pairwise_relations.get('1'), known, node_list, results, 0)
    #print(node_list)
    results.append(_visit_nodes(node_list))
    return results

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
