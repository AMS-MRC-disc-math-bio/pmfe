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
    ScoreVector NNTM::score(const RNAStructure& structure, bool compute_w) const {
        RNAStructureTree tree (structure);

        // Score the structure using these parameters
        ScoreVector result = scoreTree(tree);

        // If requested, compute the w value by re-scoring the structure with the classical parameters
        if (compute_w) {
            Turner99 classical_constants;
            NNTM classical_model(classical_constants, dangles, thread_pool);
            ScoreVector classical_score = classical_model.score(structure, false);
            mpq_class classical_energy = classical_score.energy;
            result.w = classical_energy - (result.multiloops * classical_constants.multConst[0] + result.unpaired * classical_constants.multConst[1] + result.branches * classical_constants.multConst[2]);

            // Check that the computed w is consistent
            mpq_class formula_energy = result.multiloops * constants.params.multiloop_penalty + result.unpaired * constants.params.unpaired_penalty + result.branches * constants.params.branch_penalty + result.w * constants.params.dummy_scaling;
            formula_energy.canonicalize();

            if (result.energy != formula_energy) {
                throw std::logic_error("w calculation was inconsistent!");
            }
        }

        result.canonicalize();
        return result;

    }

    ScoreVector NNTM::scoreTree(const RNAStructureTree& tree) const {
        // Score the external node
        ScoreVector score = scoreE(tree);

        // Recurse through the tree
        for (std::deque<IntervalTreeNode>::const_iterator child = tree.root.children.begin(); child != tree.root.children.end(); ++child) {
            score += scoreInternalNodeRecursively(tree, *child);
        }

        return score;
    }

    ScoreVector NNTM::scoreInternalNodeRecursively(const RNAStructureTree& tree, const IntervalTreeNode& node) const {
        // TODO: Confirm that node is in tree
        ScoreVector score;

        int i = node.start;
        int j = node.end;

        // Compute contribution from the current loop
        // This depends on the type of loop, which we determine using the number of children
        switch (node.valency()) {
        case 0:
        {
            // Hairpin loop
            score.energy += eH(i, j, tree.seq);
            break;
        }

        case 1:
        {
            // This could be a stack, a bulge, or an internal loop
            const IntervalTreeNode& child = node.children.front();

            if (child.start == i + 1 and child.end == j - 1) {
                // Stack
                score.energy += eS(i, j, tree.seq);
            } else {
                // Internal or bulge
                score.energy += eL(i, j, child.start, child.end, tree.seq);
            }

            break;
        }

        default:
        {
            // At least two children means this is a multiloop
            score += scoreM(tree, node);// + auPenalty(i, j, tree.seq);
            break;
        }
        }

        // Recurse! Recurse!
        for (std::deque<IntervalTreeNode>::const_iterator child = node.children.begin(); child != node.children.end(); ++child) {
            score += scoreInternalNodeRecursively(tree, *child);
        }

        return score;
    }

    ScoreVector NNTM::scoreMUnpairedRegion(const RNAStructure& structure, int i1, int j1, int i2, int j2, bool is_external) const {
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
            d5 = Ed5(i1, j1, structure.seq, true); // Dangle at the 5' end of the region
            d3 = Ed5(i2, j2, structure.seq); // Dangle at the 3' end of the region
        } else if (i2 < i1 and i1 < j1 and j1 < j2) {
            // (i2, j2) is the initiating pair of the loop, so
            // region is between j1 and j2
            start = j1;
            end = j2;
            d5 = Ed3(i1, j1, structure.seq);
            d3 = Ed3(i2, j2, structure.seq, true);
        } else if (i1 < j1 and j1 < i2 and i2 < j2) {
            // neither (i1, j1) and (i2, j2) is the initiating pair, so
            // region is between j1 and i2
            start = j1;
            end = i2;
            d5 = Ed3(i1, j1, structure.seq); // Dangle at the 5' end of the region
            d3 = Ed5(i2, j2, structure.seq); // Dangle at the 3' end of the region
        } else {
            throw std::invalid_argument("Invalid multiloop nesting");
        }

        ScoreVector score;

        if (not is_external) {
            score.energy += (end - start - 1) * constants.multConst[1]; // Unpaired base penalty
            score.unpaired += end - start - 1;
        }

        // Determine dangling base contribution
        switch (dangles) {
        case BOTH_DANGLE:
        {
            // In BOTH_DANGLE mode, we just compute the dangle energies and add them
            // without any logic to determine whether this is reasonable
            score.energy += d5 + d3;
            break;
        }

        case CHOOSE_DANGLE:
        {
            // Apply dangling contributions as recorded in the structure
            if (structure.does_d5(end - 1)) {
                score.energy += d3;
            }

            if (structure.does_d3(start + 1)) {
                score.energy += d5;
            }

            break;
        }

        case NO_DANGLE:
        default:
            // Otherwise, no dangling contribution
            break;
        }

        return score;
    }

    ScoreVector NNTM::scoreM(const RNAStructureTree& tree, const IntervalTreeNode& node) const {
        /*
           Score a multiloop
        */

        ScoreVector score;
        // Contribution from initializing the multiloop
        score.energy += constants.multConst[0];
        score.multiloops += 1;

        // Contribution from branches
        score.energy += (node.valency() + 1) * constants.multConst[2]; // Branch penalty
        score.branches += (node.valency() + 1);

        // Unpaired bases are accounted in the eMUnpairedRegion helper method

        // Add the contribution for each unpaired region
        // First, consider the regions beginning or ending at the initiating pair
        if (node.valency() > 0) {
            score += scoreMUnpairedRegion(tree, node.start, node.end, node.children.front().start, node.children.front().end);
            score += scoreMUnpairedRegion(tree, node.children.back().start, node.children.back().end, node.start, node.end);

            // Then consider the regions between the children
            for (std::deque<IntervalTreeNode>::const_iterator child = node.children.begin(); child != node.children.end() - 1; ++child) {
                std::deque<IntervalTreeNode>::const_iterator neighbor = child +1;
                score += scoreMUnpairedRegion(tree, child->start, child->end, neighbor->start, neighbor->end);
            }

            // We also need to account for any non-GC pairs at the branches
            score.energy += auPenalty(node.start, node.end, tree.seq);

            for (std::deque<IntervalTreeNode>::const_iterator child = node.children.begin(); child != node.children.end(); ++child) {
                score.energy += auPenalty(child->start, child->end, tree.seq);
            }
        }

        return score;
    }


    ScoreVector NNTM::scoreE(const RNAStructureTree& tree) const {
        /*
          Score the external loop
        */
        ScoreVector score;

        // Add the contribution for each unpaired region
        // First, consider the dangling regions at the ends of the sequence
        bool has_5d, has_3d;
        mpq_class d3, d5;

        // If there are no children, stop immediately
        if (tree.root.valency() == 0) {
            return score;
        }

        // Otherwise, we first process any dangling ends
        const IntervalTreeNode& firstbranch = tree.root.children.front();
        const IntervalTreeNode& lastbranch = tree.root.children.back();

        if (firstbranch.start > 0) { // If there is a dangling 5' end
            has_5d = tree.does_d5(firstbranch.start - 1);
            d5 = Ed5(firstbranch.start, firstbranch.end, tree.seq);
        } else {
            has_5d = false;
            d5 = 0;
        }

        if (lastbranch.end < tree.len() - 1) { // If there is a dangling 3' end
            has_3d = tree.does_d3(lastbranch.end + 1);
            d3 = Ed3(lastbranch.start, lastbranch.end, tree.seq);
        } else {
            has_3d = false;
            d3 = 0;
        }

        switch (dangles) {
        case BOTH_DANGLE:
        {
            score.energy += d3 + d5;
            break;
        }

        case CHOOSE_DANGLE:
        {
            if (has_5d) {
                score.energy += d5;
            }

            if (has_3d) {
                score.energy += d3;
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
            score += scoreMUnpairedRegion(tree, child->start, child->end, neighbor->start, neighbor->end, true);
        }

        // We also need to account for any non-GC pairs at the branches
        for (std::deque<IntervalTreeNode>::const_iterator child = tree.root.children.begin(); child != tree.root.children.end(); ++child) {
            score.energy += auPenalty(child->start, child->end, tree.seq);
        }

        return score;
    }
}
