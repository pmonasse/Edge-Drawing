// SPDX-License-Identifier: MPL-2.0
/**
 * @file maxTree.cpp
 * @brief Generic implicit max-tree computation
 * @author Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2026
 */

// Linked to IPOL publication:
// [1] Edge Drawing: A Fast Edge Segment Detector,
// Adle Ben Salem and Pascal Monasse, IPOL, 2026

// Algorithm from
// [2] Effective component tree computation with application to pattern
// recognition in astronomical imaging
// Berger et al., ICIP, 2007

#ifndef MAX_TREE_H
#define MAX_TREE_H

#include <vector>
#include <numeric>
#include <algorithm>
#include <utility>

/// Find_root of Union/Find algorithm.
/// [1] Algorithm 9
inline int uf_root(std::vector<int>& zpar, int i) {
    if(zpar[i]==i)
        return i;
    return (zpar[i] = uf_root(zpar, zpar[i]));
}

/// Generic max-tree computation based on Union-Find (Berger et al.)
/// \param[in] n Number of nodes of the graph
/// \param[in] cmp Functor implementing a strict ordering of nodes
/// \param[in] nbh Functor returning the range of neighbors of a node
/// \param[out] iroot Node number of the tree root (if not null pointer)
/// \return Parent map
/// [1] Algorithm 8.
/// Nodes of the graph are coded by integers from 0 to n-1. \a cmp compares
/// values associated to the nodes. \a nbh provides the edges of the graph.
/// Its operator() should return a std::pair of iterators delimiting the range
/// of neighors of the node number passed as argument. There must be a type
/// Neighborhood::iterator. Typically, \a nbh would have an internal container
/// that is filled when operator() is called.
/// \a iroot returns the node number of the tree root. In the unusual case where
/// the graph is disconnected, it returns the root of one tree of the forest.
/// The parent map is interpreted as follows.
/// - The parent of any node is canonical.
/// - If a node and its parent have the same value, they are private points
///   of the same tree node.
/// - A root (could be several if the graph is diconnected) is its own parent.
template <typename Comparator, typename Neighborhood>
std::vector<int> max_tree(int n, const Comparator& cmp, Neighborhood& nbh,
                          int* iroot=0) {
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), cmp);
    if(iroot)
        *iroot = idx[0];
    std::vector<int> par(n,-1);
    std::vector<int> zpar(n,-1);
    // Build tree
    for(int i=(int)n-1; i>=0; i--) {
        int j=idx[i];
        par[j] = zpar[j] = j;
        using It = typename Neighborhood::iterator;
        for(std::pair<It,It> p=nbh(j); p.first!=p.second; ++p.first)
            if(zpar[*p.first]>=0) {
                int k = uf_root(zpar,*p.first);
                par[k] = zpar[k] = j;
            }
    }
    // Canonize
    for(int i=1; i<n; i++) {
        int j=idx[i], k=par[j];
        if(! cmp(par[k], k)) // Test equality
            par[j] = par[k];
    }
    return par;
}

#endif
