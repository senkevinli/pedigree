#!/usr/bin/env python3

""" For traversing the pedigree and other related operations. """

from __future__ import annotations
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

        # # Give two full siblings of each gender.
        # Node.filler_id += 1
        # fem_sibling = Node(
        #     str(Node.filler_id),
        #     True,
        #     mt_dna=self.mt_dna,
        #     parents=self.parents,   
        # )
        # Node.filler_id += 1
        # male_sibling = Node(
        #     str(Node.filler_id),
        #     False,
        #     mt_dna=self.mt_dna,
        #     y_chrom=self.parents[1].y_chrom,
        #     parents=self.parents
        # )

        # fem_sibling.siblings += [male_sibling, self]
        # male_sibling.siblings += [male_sibling, self]

        # for parent in self.parents:
        #     parent.children += [male_sibling, fem_sibling]

        # # Make filler partner.
        # Node.filler_id += 1
        # partner = Node(
        #     str(Node.filler_id),
        #     not self.female,
        #     partners=[self]
        # )
        # self.partners.append(partner)

        # # Give two children of each gender.
        # Node.filler_id += 1
        # fem_child = Node(
        #     str(Node.filler_id),
        #     True,
        #     mt_dna=self.mt_dna if self.female else None,
        #     parents=(self, partner) if self.female else (partner, self) 
        # )
        # Node.filler_id += 1
        # male_child = Node(
        #     str(Node.filler_id),
        #     False,
        #     mt_dna=self.mt_dna if self.female else None,
        #     y_chrom=self.y_chrom if not self.female else None,
        #     parents=(self, partner) if self.female else (partner, self)
        # )

        # partner.children += [fem_child, male_child]
        # self.children += [fem_child, male_child]

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
            pass
        else:
            # No configuration works here.
            return False

    # Case that one is female and the other is male.
    elif (not src.female and dest.female) or \
         (src.female and not dest.female):
            male_node = src if dest.female else dest
            female_node = src if src.female else dest
            print(male_node)
            print(female_node)

            if share_mt_dna:
                # Either siblings or son/mother.
                
                # Case 1 siblings.
                pass
            else:
                orig_female_node_parents = female_node.parents
                # Must be father/daughter.
                if female_node.parents[1].occupied and \
                   female_node.parents[1] != male_node:
                    return # Configuration is impossible.
                elif female_node.parents[1] == male_node:
                    _assign_helper(relation, node_map, node_list, all_possible, idx + 1)
                else:
                    orig_female_node_parents = female_node.parents
                    mother = female_node.parents[0]
                    orig_mother_children = [child for child in mother.children]
                    orig_father_children = [child for child in male_node.children]
                    female_node.parents = (mother, male_node)

                    # Add children of daughter's siblings.
                    for child in mother.children:
                        if child not in male_node.children:
                            child.parents = (mother, male_node)
                            male_node.children.append(child)
                    for child in male_node.children:
                        if child not in mother.children:
                            child.parents = (mother, male_node)
                            mother.children.append(child)
                    # Recurse.
                    _assign_helper(relation, node_map, node_list, all_possible, idx + 1)

                    female_node.parents = orig_female_node_parents
                    for child in male_node.children:
                        if child in orig_mother_children:
                            child.parents = orig_female_node_parents 
                            male_node.children.remove(child)

                    for child in mother.children:
                        if child in orig_father_children:
                            child.parents = (male_node.children[0].parents[0], male_node)
                            mother.children.remove(child)

    # Case that source and dest are both females.
    else:
        if share_mt_dna:
            # May be siblings or daughter/mother or mother/daughter.
            pass
        else:
            # No configuration works here.
            return
    
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

    # for degree in node_map.keys():
    #     degree_nodes = node_map.get(degree)
    #     for src, dest in degree_nodes:
    #         src = node_map.get(src)
    #         dest = node_map.get(dest)

    #         both_male = not src.female and not dest.female
    #         one_female = (src.female and not dest.female) or \
    #                      (dest.female and not src.female)

    #         if both_male:
    #             pass
    #         elif one_female:
    #             pass
    #         else:
    #             pass


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
    _assign_helper(pairwise_relations.get('1'), known, node_list, results, 0)
    print(results)
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
