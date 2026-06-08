// SPDX-License-Identifier: MPL-2.0
/**
 * @file chain_tree.h
 * @brief Binary tree of exploration
 * @author Adle Ben Salem
 *         Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2025-2026
 */

#ifndef CHAIN_TREE_H
#define CHAIN_TREE_H

#include "image.h"
#include <vector>
#include <stack>

/// Direction of exploration.
enum Direction { LEFT, RIGHT, UP, DOWN };

/// A node of the binary tree of exploration.
struct ChainTree {
    std::vector<Point> pts; /// Successive points
    ChainTree* child[2]; /// Child nodes
    int len; /// Max length to leaf
    int path; /// Child in the path to most distant leaf
    const Direction dir; /// Direction of exploration

    ChainTree();
    ChainTree(Direction direction, ChainTree* parent);
    ~ChainTree();

    int length();
    void pruneLongestPath(std::stack<ChainTree*>& orphans);
};

/// Anchor for exploration, stored in a stack.
struct StackNode {
    Point pos;
    Direction dir;         // Direction of exploration
    ChainTree* parent;

    StackNode(Point p, Direction d, ChainTree* par)
    : pos(p), dir(d), parent(par) {}
};

#endif
