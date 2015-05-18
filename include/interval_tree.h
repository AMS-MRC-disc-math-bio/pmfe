#ifndef _INTERVAL_TREE_H
#define _INTERVAL_TREE_H

#include <deque>
#include <stdexcept>

class IntervalTreeNode {
    /*
      Node in an interval search tree

      Each node is an interval in some range, with the requirement
      that each node's interval is a strict subinterval of its parent's
    */
public:
    int start;
    int end;
    std::deque<IntervalTreeNode> children;

IntervalTreeNode(int start, int end): start(start), end(end) {};
IntervalTreeNode(): start(0), end(0) {}; // Default constructor to appease the compiler: will behave badly if used without reinitializing!

    void insert(const int start, const int end, const bool sort = false) {
        /*
          Insert a new node, either as a child here or by handing off to a child
        */
        if (start <= this->start or end >= this->end) {
            throw std::invalid_argument("Could not insert node");
        }

        if (start >= end) {
            throw std::invalid_argument("Invalid interval");
        }

        // Check whether [start, end] is a subinterval of any child and, if so, insert it there
        for (auto& child: children) {
            if ((child.start < start) and (end < child.end)) {
                child.insert(start, end);
                return;
            }
        }

        // If the interval was not contained in a child, we insert it here
        IntervalTreeNode child (start, end);
        children.push_back(child);

        if (sort) {
                this->sort();
            }

        return;
    }

    bool operator<(const IntervalTreeNode& b) const {
        return ((this->start < b.start) or (this->start == b.start and this->end < b.end));
    }

    void sort() {
        /*
          Sort the children in increasing order of start base position
        */
        std::sort(children.begin(), children.end());

        // Recurse! Recurse!
        for (auto& child: children) {
            child.sort();
        }
    }

    const int valency() const {
        /*
          Return the number of children of this node
        */
        return children.size();
    }
};

#endif
