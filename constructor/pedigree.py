#!/usr/bin/env python3

""" For traversing the pedigree and other related operations. """

from __future__ import annotations
from typing import Tuple, List, Optional

AGE_DEFAULT = 100

class Node:
    def __init__(
        self,
        id: str,
        female: bool,
        mt_dna: str,
        y_chrom: Optional[str] = None,
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

        # Tuple, first element is mother, second is father.
        if parents is not None:
            assert parents[0].female
            assert not parents[1].female

        self.parents = parents
        self.children = children

        if female:
            assert y_chrom is None

        self.y_chrom = y_chrom
        self.siblings = siblings
        self.partners = partners
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
