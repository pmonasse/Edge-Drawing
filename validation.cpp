// SPDX-License-Identifier: MPL-2.0
/**
 * @file validation.cpp
 * @brief NFA validation of edge drawing
 * @author Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2026
 */

#include "ED.h"
#include "maxTree.h"
#include <cmath>
#include <cassert>

/// Test if path is closed, that is, first and last points are 8-neighbors.
static bool closed(const std::vector<Point>& e) {
    short dx = (short)e.front().x - (short)e.back().x,
          dy = (short)e.front().y - (short)e.back().y;
    return dx*dy*dx*dy==1;
}

/// Max-tree of edge intervals.
/// Stores bounds \c beg and \c end, a value \a v (initially the min gradient
/// on the interval, later log10 of NFA). If the edge is circular, the interval
/// may loop and then end < beg.
struct Interval {
    Interval* parent;
    std::vector<Interval*> child;
    int beg, end;
    float v; ///< First min gradient, then log10(NFA)

    Interval(int i, float v0): parent(0), beg(i), end(i), v(v0) {}
    ~Interval();
    bool loop() const { return end<beg; }
    void addChild(Interval* c);
    void add(int i, int ext);
    void fillBounds(int ext);
    const Interval* findMinValue() const;
    std::vector<Point> segment(const std::vector<Point>& e) const;
};

/// Destructor. Recursive deallocation of all subtree.
Interval::~Interval() {
    std::vector<Interval*>::iterator it, end=child.end();
    for(it=child.begin(); it!=end; ++it)
        delete *it;
}

/// Add \a c as child.
void Interval::addChild(Interval* c) {
    c->parent = this;
    child.push_back(c);
}

/// Add index \a i, adjusting \c beg and \c end of interval. If the interval is
/// looping, it needs \a ext, a point outside the interval.
void Interval::add(int i, int ext) {
    if(i<beg && (!loop() || ext<i))
        beg = i;
    if(i>end && (!loop() || i<ext))
        end = i;
}

/// Integrate sub-intervals to update fields \c beg and \c end.
void Interval::fillBounds(int ext) {
    std::vector<Interval*>::iterator it, end=child.end();
    for(it=child.begin(); it!=end; ++it) {
        (*it)->fillBounds(ext);
        add((*it)->beg, ext);
        add((*it)->end, ext);
    }
}

/// Find interval with minimal value \c v in all subtree.
const Interval* Interval::findMinValue() const {
    const Interval* min = this;
    std::vector<Interval*>::const_iterator it, end=child.end();
    for(it=child.begin(); it!=end; ++it) {
        const Interval* m = (*it)->findMinValue();
        if(m->v < min->v)
            min = m;
    }
    return min;
}

/// Extract segment from full edge.
std::vector<Point> Interval::segment(const std::vector<Point>& e) const {
    std::vector<Point> v;
    if(loop()) {
        v.insert(v.end(), e.begin()+beg, e.end());
        v.insert(v.end(), e.begin(), e.begin()+end+1);
    } else
        v.insert(v.end(), e.begin()+beg, e.begin()+end+1);
    return v;
}

/// Find index of minimum of \a G among points in \a e.
int min(const Image<int>& G, const std::vector<Point>& e) {
    const int n=(int)e.size();
    int min=G(e[0]);
    for(int i=1; i<n; i++)
        if(min<G(e[i]))
            min=G(e[i]);
    return min;
}

/// A contrario validation. \a lEpsNFA is the log10 of detection threshold.
/// Its normal value is 0, or negative for more requiring detection.
/// @param lEpsNFA log10 of max NFA for validation. 
/// @param segLevel Segmentation level (0,1,other).
void ED::validateNFA(float lEpsNFA, int segLevel) {
    for(Point p={1,1}; p.y+1<S.h; p.y++)
        for(p.x=1; p.x+1<S.w; p.x++)
            S(p) = G(p)<minGrad? 0: ANCHOR;
    std::vector<int> H = cumulHistoGradAnchors();
    if(H.empty()) {
        edges.clear();
        return;
    }

    std::vector<float> lProba(H.size(), 0);
    const int n = H.back();
    const float v = std::log10(n);
    for(size_t i=1; i<H.size(); i++)
        lProba[i] = log10(n-H[i-1])-v;

    int nTests = 0;
    std::vector<std::vector<Point>>::const_iterator it, end=edges.end();
    for(it=edges.begin(); it!=end; ++it) {
        int n=(int)it->size();
        nTests += closed(*it)? n*(n-1)+1: n*(n+1)/2;
    }
    const float lTestsWhole=log10(edges.size()), lTests = log10(nTests);

    std::vector<std::vector<Point>> valid;
    for(it=edges.begin(); it!=end; ++it) {
        if(segLevel==0 || segLevel==1) { // Test whole line valid
            int gmin = min(G, *it);
            if(lTestsWhole+it->size()*0.5f*lProba[gmin] <= lEpsNFA) {
                valid.push_back(*it);
                continue;
            }
        }
        if(segLevel!=0)
            validateEdge(*it, lProba, lTests, lEpsNFA, segLevel!=1, valid);
    }
    std::swap(edges, valid);
}

/// Functor for sorting based on gradient along edge.
struct CompareGradEdge {
    const Image<int>& G;
    const std::vector<Point>& E;
    CompareGradEdge(const Image<int>& g, const std::vector<Point>& e)
    : G(g), E(e) {}
    bool operator()(int i, int j) const { return (G(E[i])<G(E[j])); }
};

/// Functor for neighborhood of point inside an edge.
struct NeighborhoodEdge {
    typedef int* iterator;
    int n; bool circular; int nbh[2];
    NeighborhoodEdge(int n0, bool circ): n(n0), circular(circ) {}
    std::pair<iterator,iterator> operator()(int j) {
        std::pair<iterator,iterator> p(nbh,nbh);
        int l = j-1;
        if(circular && l<0)
            l = (int)n-1;
        if(l>=0)
            *p.second++ = l;
        l = j+1;
        if(circular && l>=(int)n)
            l = 0;
        if(l<(int)n)
            *p.second++ = l;
        return p;
    }
};

/// Create nodes of the max-tree. It creates one node for each canonical element
/// in \a par. Index corresponding to non-canonical elements are null pointers.
std::vector<Interval*> build_tree_nodes(const Image<int>& G,
                                        const std::vector<Point>& e,
                                        const std::vector<int>& par, int root) {
    int n = (int)par.size();
    std::vector<Interval*> tree(n, 0);
    tree[root] = new Interval(root, G(e[root]));
    for(int i=0; i<n; i++) { // Build tree nodes
        int v = G(e[i]);
        if(G(e[par[i]]) != v)
            tree[i] = new Interval(i,v);
    }
    return tree;
}

/// Create edges of the max-tree.
void build_tree_edges(std::vector<Interval*>& tree,
                      const std::vector<int>& par, int root) {
    int n = (int)par.size();
    for(int i=0; i<n; i++)  // Build tree edges
        if(i!=root && tree[i])
            tree[par[i]]->addChild(tree[i]);
}

/// Lowest common ancestor on max-tree.
/// Rely on increasing value while going down-tree.
static Interval* lca(Interval* i1, Interval* i2) {
    while(i1 != i2) {
        assert(i1 && i2);
        if(i2->v <= i1->v)
            i1 = i1->parent;
        else
            i2 = i2->parent;
    }
    return i1;
}

/// Adjust field beg or end for circular intervals.
void tag_circular_intervals(std::vector<Interval*>& tree,
                            const std::vector<int>& par, int ext) {
    Interval* i1 = tree[0]? tree[0]: tree[par[0]];
    Interval* i2 = tree.back()? tree.back(): tree[par[tree.size()-1]];
    for(i1 = lca(i1,i2); i1->parent; i1=i1->parent)
        if(i1->beg < ext)
            i1->beg=(int)tree.size()-1;
        else
            i1->end=0;
}

/// Fill bounds beg and end of every interval.
void fill_bounds(std::vector<Interval*>& tree,
                 const std::vector<int>& par, int root, int ext) {
    int n = (int)par.size();
    for(int i=0; i<n; i++) // Fill bounds (without sub-intervals)
        if(! tree[i])
            tree[par[i]]->add(i, ext);
    tree[root]->fillBounds(ext); // Integrate sub-intervals in computing bounds
    tree[root]->beg=0; tree[root]->end=(int)n-1; // Fix root bounds
}

/// Step 3 of algorithm in ED::validateEdge.
void extract_valid_segments(const std::vector<Point>& e,
                            const Interval* r, float lEpsNFA, bool bSubLines,
                            std::vector<std::vector<Point>>& valid) {
    const Interval* m = r->findMinValue();
    if(m->v > lEpsNFA)
        return;
    if(! bSubLines) {
        valid.push_back(e);
        return;
    }
    valid.push_back(m->segment(e));
    for(; m!=r; m = m->parent) {
        std::vector<Interval*>::const_iterator it, end=m->parent->child.end();
        for(it=m->parent->child.begin(); it!=end; ++it)
            if(*it!=m)
                extract_valid_segments(e, *it, lEpsNFA, bSubLines, valid);
    }
}

/// Append to \a valid the maximally contrasted segments of \a e.
/// \a lProba gives the log10 probability of contrast at least index.
/// \a lTests is log10 of the number of tests and \a lEpsNFA is log10 of the
/// upper bound threshold for meaningfulness.
/// Algo:
/// 1. Compute the max-tree of gradients on \a e (Berger algorithm).
/// 2. Find most meaningful segment of tree.
/// 3. If meaningful, validate and go back to 2 for all disjoint segments.
void ED::validateEdge(const std::vector<Point>& e,
                      const std::vector<float>& lProba,
                      float lTests, float lEpsNFA, bool bSubLines,
                      std::vector<std::vector<Point>>& valid) const {
    const int n=(int)e.size();
    const bool circular = closed(e);
    CompareGradEdge cmp(G,e);
    NeighborhoodEdge nbh(n,circular);
    int root;
    std::vector<int> par = max_tree(n, cmp, nbh, &root);

    std::vector<Interval*> tree = build_tree_nodes(G, e, par, root);
    build_tree_edges(tree, par, root);
    const int ext=tree[root]->beg; // Index of point outside any loop
    if(circular)
        tag_circular_intervals(tree, par, ext);
    fill_bounds(tree, par, root, ext);
    for(int i=0; i<n; i++) // Compute log NFA
        if(tree[i]) {
            int len = (tree[i]->loop()? (int)n-tree[i]->beg+tree[i]->end+1:
                       tree[i]->end-tree[i]->beg+1);
            tree[i]->v = lTests+len*0.5f*lProba[tree[i]->v];
        }
    extract_valid_segments(e, tree[root], lEpsNFA, bSubLines, valid);
    delete tree[root];
}
