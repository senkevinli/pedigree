#!/usr/bin/env python3

from .pedigree import Node
from typing import List

class Graph:
    def __init__(node_list: List[Node]) -> None:
        """
            Constructs a graph that is a deep copy of the given node list.
        """