#!/usr/bin/env python3

""" For traversing the pedigree and other related operations. """

from __future__ import annotations
from typing import Tuple, List, Optional

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

    def extrapolate_node(self):
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