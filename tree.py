class PartitionIntervalNode:
    """
    Node in a partition interval search tree

    Partition interval search trees require that each interval except the root
    be a strict subinterval of its parent (the "nestedness" property).
    If this assumption is not satisfied, we will raise an exception.
    """
    def __init__(self, start, end):
        """
        Node constructor

        @param start start point of interval
        @param end end point of interval
        """
        self.start = start
        self.end = end
        self.children = []

    def __repr__(self):
        return "Node ({0}, {1})".format(self.start, self.end)

    def print_tree(self, level=0):
        """
        Print a schematic representation of the tree descending from `self`
        """
        # This is super-ugly right now, but it gets the job done for testing
        print " "*level + u"\u2520\u2500" + str(self)
        if len(self.children) != 0:
            print " "*level + u"\u2517\u2513"
            for child in self.children:
                child.print_tree(level+1)

    def insert(self, start, end):
        """
        Insert new note as child of this node

        @param start start point of interval
        @param end end point of interval
        """
        # If nestedness fails, raise an exception
        if start <= self.start or end >= self.end or start >= end:
            raise ValueError("Nesting failure! Parent ({0}, {1}), child ({2}, {3})".format(self.start, self.end, start, end))

        assigned = False
        for child in self.children:
            # This is where we assume that nestedness holds
            if child.start < start < child.end:
                child.insert(start, end)
                assigned = True
                break

        if not assigned:
            self.children.append( PartitionIntervalNode(start, end) )

    def valency(self):
        return len(self.children)

    def sort(self):
        self.children.sort()
        for child in self.children:
            child.sort()
