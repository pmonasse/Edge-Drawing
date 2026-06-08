// SPDX-License-Identifier: MPL-2.0
/**
 * @file chain_tree.cpp
 * @brief Binary tree of exploration
 * @author Adle Ben Salem
 *         Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2025-2026
 */

#include "chain_tree.h"

/// Constructor
ChainTree::ChainTree()
: len(0), path(-1), dir(LEFT) {
    child[0] = child[1] = nullptr;
}

/// Constructor with direction and parent
ChainTree::ChainTree(Direction direction, ChainTree* p)
: len(0), path(-1), dir(direction) {
    child[0] = child[1] = nullptr;
    if(p) {
        int i = p->child[0]? 1: 0;
        p->child[i] = this;
    }
}

/// Destructor: root destroys all subtree.
ChainTree::~ChainTree() {
    delete child[0];
    delete child[1];
}

/// Length to most distant leaf.
/// \details Compute and store the path to the most distant leaf. The distance
/// is the sum of points of nodes in the path. The fields \a len and \a path are
/// set for this and all descendents. \a path is 0 or 1 indicating the index of
/// child to the most distant leaf.
int ChainTree::length() {
    int l[2];
    l[0] = child[0] ? child[0]->length() : 0;
    l[1] = child[1] ? child[1]->length() : 0;
    path = (l[0]>=l[1])? 0: 1;
    len = pts.size() + l[path];
    return len;
}

/// Follow path to most distant leaf. The children not followed are roots of
/// sub-trees and are gathered in \a oprhans. The method \a length must have
/// been called, so that fields \a len and \a path are set.
void ChainTree::pruneLongestPath(std::stack<ChainTree*>& orphans) {
    if(! child[path])
        return;
    if(child[1-path])
        orphans.push(child[1-path]);
    child[path]->pruneLongestPath(orphans);
    pts.insert(pts.end(), child[path]->pts.begin(), child[path]->pts.end());
}
