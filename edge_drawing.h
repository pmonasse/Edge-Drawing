// SPDX-License-Identifier: MPL-2.0
/**
 * @file edge_drawing.h
 * @brief edge drawing
 * @author Adle Ben Salem
 *         Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2025-2026
 */

#ifndef EDGE_DRAWING_H
#define EDGE_DRAWING_H

#include "chain_tree.h"
#include "image.h"
class StackNode;

/// Edge drawing algorithm.
class ED {
public:
    ED(const Image<int>& G, const Image<float>& Theta,
       int gradMin=6, int anchorGap=2, int minPathLen=10);
    void validateNFA(float lEpsNFA=0, int segLevel=0);

    std::vector<std::vector<Point>> edges;

    typedef bool Orientation;
    typedef char State;
    static const Orientation HORIZONTAL=true;
    static const Orientation VERTICAL=false;
    static const State ANCHOR=1;
    static const State EDGE=2;
protected:
    Image<int> G;
    Image<Orientation> O;
    Image<State> S;
    int minGrad;
    int minLen;

private:
    void computeAnchors(int anchorThresh);
    std::vector<int> cumulHistoGradAnchors() const;
    std::vector<Point> sortedAnchors() const;
    void joinAnchors();
    void exploreChain(StackNode&, ChainTree*, std::stack<StackNode>&);
    bool nextPixelChain(StackNode& node);
    void buildChainTree(ChainTree* root, Point p);
    void extractEdgesFromTree(ChainTree* root);
    void buildRootEdge(ChainTree* root);
    void validateEdge(const std::vector<Point>& e,
                      const std::vector<float>& lProba,
                      float lTests, float lEpsNFA, bool bSubLines,
                      std::vector<std::vector<Point>>& valid) const;
};

#endif
