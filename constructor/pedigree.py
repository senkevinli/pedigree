#!/usr/bin/env python3

from __future__ import annotations
from _typeshed import SupportsNoArgReadline
from typing import Tuple, List, Optional
""" For traversing the pedigree and other related operations. """

AGE_DEFAULT = 100

class Node:
    def __init__(
        self,
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
    