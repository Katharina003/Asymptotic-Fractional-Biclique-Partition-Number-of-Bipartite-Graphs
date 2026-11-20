#pragma once

#include "CondensedNode.h"
#include "CondensedFormula.h"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

/*
Data structure to represent bicliques.
Top nodes are what we call left nodes in the thesis; bottom nodes are what we call right nodes in the thesis.
This data structure provides biclique witnesses for the vectors of M in DominoGraphs.
For this reason, the top nodes and bottom nodes always have to be set, as they are needed for a conclusive witness.
*/
class Biclique {
private:

    std::vector<CondensedNode> topNodesCondensed; //left nodes
    std::vector<CondensedNode> bottomNodesCondensed; //right nodes
    CondensedFormula condensedFormula; //condensed formula calculated from left and right nodes
    bool isCondensedFormulaSet=false; //variable that is needed for condensedFormula lookup. If condensed formula is not set, calculate it from top and bottom nodes.

    void topBottomToEdges();
public:
    Biclique(std::vector<CondensedNode> topNodesNew, std::vector<CondensedNode> bottomNodesNew);

    const std::vector<CondensedNode> getTop() const;
    const std::vector<CondensedNode> getBot() const;
    const CondensedFormula getCondensedFormula();
    bool operator==(const Biclique& other) const;
    Biclique kroneckerProduct(Biclique& b) const;

    std::string to_string();
};

namespace std {
template<>
struct hash<Biclique> {
    std::size_t operator()(const Biclique& b) const;
};
}
