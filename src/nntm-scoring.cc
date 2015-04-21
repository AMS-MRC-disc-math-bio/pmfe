#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdexcept>

#include "nntm.h"
#include "nndb_constants.h"
#include "interval_tree.h"

#include <gmpxx.h>
#include "boost/multi_array.hpp"

#include <vector>

namespace pmfe {
    const mpq_class NNTM::energy(const RNAStructure& structure) const {
        RNAStructureTree tree (structure);
        return scoreTree(tree);
    }

    const mpq_class NNTM::scoreTree(const RNAStructureTree& tree) const {
        // Score the external node
        mpq_class energy = eE(tree);

        // Recurse through the tree
        for (std::deque<IntervalTreeNode>::const_iterator child = tree.root.children.begin(); child != tree.root.children.end(); ++child) {
            energy += scoreInternalNodeRecursively(tree, *child);
        }

        return energy;
    }

    const mpq_class NNTM::scoreInternalNodeRecursively(const RNAStructureTree& tree, const IntervalTreeNode& node) const {
        // TODO: Confirm that node is in tree
        mpq_class energy = 0;

        int i = node.start;
        int j = node.end;

        // Compute contribution from the current loop
        // This depends on the type of loop, which we determine using the number of children
        switch (node.valency()) {
        case 0:
        {
            // Hairpin loop
            energy += eH(i, j, tree.seq);
            break;
        }

        case 1:
        {
            // This could be a stack, a bulge, or an internal loop
            const IntervalTreeNode& child = node.children.front();

            if (child.start == i + 1 and child.end == j - 1) {
                // Stack
                energy += eS(i, j, tree.seq);
            } else {
                // Internal or bulge
                energy += eL(i, j, child.start, child.end, tree.seq);
            }

            break;
        }

        default:
        {
            // At least two children means this is a multiloop
            mpq_class bonus = eM(tree, node);// + auPenalty(i, j, tree.seq);
            energy += bonus;
        }
        }

        // Recurse! Recurse!
        for (std::deque<IntervalTreeNode>::const_iterator child = node.children.begin(); child != node.children.end(); ++child) {
            energy += scoreInternalNodeRecursively(tree, *child);
        }

        return energy;
    }

    const mpq_class NNTM::eMUnpairedRegion(const RNAStructure& structure, int i1, int j1, int i2, int j2, bool is_external) const {
        /*
          Helper method to compute the energy of an unpaired region in a multiloop which lies between pairs (i1, j1) and (i2, j2)
          Note that this requires i1 < i2 < j2 < j1 if (i1, j1) initiates the loop,
          i2 < i1 < j1 < j2 if (i2, j2) initiates the loop,
          and  i1 < j1 < i2 < j2 if neither pair is initial.
        */

        int start, end;
        mpq_class d5, d3;

        if (i1 < i2 and i2 < j2 and j2 < j1) {
            // (i1, j1) is the initiating pair of the loop, so
            // region is between i1 and i2
            start = i1;
            end = i2;
            d5 = Ed5(i1, j1, i1+1, structure.seq); // Dangle at the 5' end of the region
            d3 = Ed3(j2, i2, i2-1, structure.seq); // Dangle at the 3' end of the region
        } else if (i2 < i1 and i1 < j1 and j1 < j2) {
            // (i2, j2) is the initiating pair of the loop, so
            // region is between j1 and j2
            start = j1;
            end = j2;
            d5 = Ed5(j1, i1, j1+1, structure.seq);
            d3 = Ed3(i2, j2, j2-1, structure.seq);
        } else if (i1 < j1 and j1 < i2 and i2 < j2) {
            // neither (i1, j1) and (i2, j2) is the initiating pair, so
            // region is between j1 and i2
            start = j1;
            end = i2;
            d5 = Ed5(j1, i1, j1+1, structure.seq); // Dangle at the 5' end of the region
            d3 = Ed3(j2, i2, i2-1, structure.seq); // Dangle at the 3' end of the region
        } else {
            throw std::invalid_argument("Invalid multiloop nesting");
        }

        mpq_class energy = 0;

        if (not is_external) {
            (end - start - 1) * constants.multConst[2]; // Unpaired base penalty
        }

        // Determine dangling base contribution
        switch (dangles) {
        case BOTH_DANGLE:
        {
            // In BOTH_DANGLE mode, we just compute the dangle energies and add them
            // without any logic to determine whether this is reasonable
            energy += d5 + d3;
            break;
        }

        case CHOOSE_DANGLE:
        {
            // Apply dangling contributions as recorded in the structure
            if (structure.does_d5(end - 1)) {
                energy += d3;
            }

            if (structure.does_d3(start + 1)) {
                energy += d5;
            }

            break;
        }

        case NO_DANGLE:
        default:
            // Otherwise, no dangling contribution
            break;
        }

        return energy;
    }

    const mpq_class NNTM::eM(const RNAStructureTree& tree, const IntervalTreeNode& node) const {
        /*
           Score a multiloop
        */

        mpq_class energy = constants.multConst[0]; // Initiation penalty
        energy += (node.valency() + 1) * constants.multConst[2]; // Branch penalty
        // Unpaired bases are accounted in the eMUnpairedRegion helper method

        // Add the energy contribution for each unpaired region
        // First, consider the regions beginning or ending at the initiating pair
        // TODO: Verify that there are children!
        energy += eMUnpairedRegion(tree, node.start, node.end, node.children.front().start, node.children.front().end);
        energy += eMUnpairedRegion(tree, node.children.back().start, node.children.back().end, node.start, node.end);

        // Then consider the regions between the children
        for (std::deque<IntervalTreeNode>::const_iterator child = node.children.begin(); child != node.children.end() - 1; ++child) {
            std::deque<IntervalTreeNode>::const_iterator neighbor = child +1;
            energy += eMUnpairedRegion(tree, child->start, child->end, neighbor->start, neighbor->end);
        }

        // We also need to account for any non-GC pairs at the branches
        energy += auPenalty(node.start, node.end, tree.seq);

        for (std::deque<IntervalTreeNode>::const_iterator child = node.children.begin(); child != node.children.end(); ++child) {
            energy += auPenalty(child->start, child->end, tree.seq);
        }

        return energy;
    }


    const mpq_class NNTM::eE(const RNAStructureTree& tree) const {
        /*
          Score the external loop
        */
        mpq_class energy = 0;

        // Add the energy contribution for each unpaired region
        // First, consider the dangling regions at the ends of the sequence
        // TODO: Verify that there are children!
        bool has_5d, has_3d;
        mpq_class d3, d5;

        const IntervalTreeNode& firstbranch = tree.root.children.front();
        const IntervalTreeNode& lastbranch = tree.root.children.back();

        if (firstbranch.start > 0) { // If there is a dangling 5' end
            has_5d = tree.does_d5(firstbranch.start - 1);
            d5 = Ed3(firstbranch.end, firstbranch.start, firstbranch.start - 1, tree.seq);
        } else {
            has_5d = false;
            d5 = 0;
        }

        if (lastbranch.end < tree.len() - 1) { // If there is a dangling 3' end
            has_3d = tree.does_d3(lastbranch.end + 1);
            d3 = Ed5(lastbranch.end, lastbranch.start, lastbranch.end + 1, tree.seq);
        } else {
            has_3d = false;
            d3 = 0;
        }

        switch (dangles) {
        case BOTH_DANGLE:
        {
            energy += d3 + d5;
            break;
        }

        case CHOOSE_DANGLE:
        {
            if (has_5d) {
                energy += d5;
            }

            if (has_3d) {
                energy += d3;
            }

            break;
        }

        case NO_DANGLE:
        default:
        {
            break;
        }
        }

        // Then consider the regions between the children
        for (std::deque<IntervalTreeNode>::const_iterator child = tree.root.children.begin(); child != tree.root.children.end() - 1; ++child) {
            std::deque<IntervalTreeNode>::const_iterator neighbor = child + 1;
            energy += eMUnpairedRegion(tree, child->start, child->end, neighbor->start, neighbor->end, true);
        }

        // We also need to account for any non-GC pairs at the branches
        for (std::deque<IntervalTreeNode>::const_iterator child = tree.root.children.begin(); child != tree.root.children.end(); ++child) {
            energy += auPenalty(child->start, child->end, tree.seq);
        }

        return energy;
    }
}
