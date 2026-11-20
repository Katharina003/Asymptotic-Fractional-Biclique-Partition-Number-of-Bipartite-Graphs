#include "Biclique.h"

#include <algorithm>

using namespace std;


/*
Constructs a condensed formula from left (top) and right (bottom) nodes
*/
void Biclique::topBottomToEdges() {
    for (auto& top : topNodesCondensed) {
        for (auto& bot : bottomNodesCondensed) {
            EdgeOrbit result = EdgeOrbit::topBottomToEdge(top, bot);
            condensedFormula.incrementEdgeCount(result, top.getNumOccurences() * bot.getNumOccurences());
        }
    }
}

/*
Constructs a biclique from left (top) and right (bottom) nodes.
Condenses (and therefore also sorts) nodes for deterministic configuration and comparison of bicliques.
*/
Biclique::Biclique(vector<CondensedNode> topNodesNew, vector<CondensedNode> bottomNodesNew) {
    CondensedNode::condenseNodeList(topNodesNew);
    CondensedNode::condenseNodeList(bottomNodesNew);
    topNodesCondensed = topNodesNew;
    bottomNodesCondensed = bottomNodesNew;
}

/*
returns all left (top) nodes
*/
const vector<CondensedNode> Biclique::getTop() const {
    return topNodesCondensed;
}

/*
returns all right (bottom) nodes
*/
const vector<CondensedNode> Biclique::getBot() const {
    return bottomNodesCondensed;
}

/*
If it is not already constructed, constructs the condensed formula, stores it and returns it.
If already constructed, returns the condensed formula.
*/
const CondensedFormula Biclique::getCondensedFormula() {
    if (!isCondensedFormulaSet) topBottomToEdges();
    isCondensedFormulaSet = true;
    return condensedFormula;
}

/*
As changing top and bottom nodes individually after initialization is not possible
and they are sorted at initialization, this comparison always works.
*/
bool Biclique::operator==(const Biclique& other) const {
    return topNodesCondensed == other.getTop() && bottomNodesCondensed == other.getBot();
}

/*
Returns a new biclique that is the Kronecker product of biclique b and this biclique.
*/
Biclique Biclique::kroneckerProduct(Biclique& b) const {
    vector<CondensedNode> top = b.getTop();
    vector<CondensedNode> bot = b.getBot();
    //the method nodeProductUnordered does not reorder the bits within node names during the process of computing the Kronecker product.
    //this way, the biclique stays a representative.
    vector<CondensedNode> newTop = CondensedNode::nodeProductUnordered(top, topNodesCondensed);
    vector<CondensedNode> newBot = CondensedNode::nodeProductUnordered(bot, bottomNodesCondensed);
    return Biclique(newTop, newBot);
}

/*
Format:
Condensed Formula: <condensed formula>
Top Nodes: <top nodes>
Bot Nodes: <bot nodes>
*/
string Biclique::to_string() {
    string s = "{\n";

    s += "Condensed Formula : ";
    s += getCondensedFormula().to_string()+"\n";
    s += "Top Nodes: ";
    for (CondensedNode n : topNodesCondensed) {
        s += n.to_string();
    }
    s += "\nBot Nodes: ";
    for (CondensedNode n : bottomNodesCondensed) {
        s += n.to_string();
    }
    s += "\n}";

    return s;
}

/*
uses top and bottom nodes for hashes, as the condensed formula for two bicliques with the same nodes is the same.
*/
std::size_t std::hash<Biclique>::operator()(const Biclique& b) const {
    size_t seed = 0;

    auto hashCombine = [](size_t& seed, size_t value) {
        seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        };
    for (const auto& n : b.getTop()) {
        hashCombine(seed, std::hash<CondensedNode>()(n));
    }

    for (const auto& n : b.getBot()) {
        hashCombine(seed, std::hash<CondensedNode>()(n));
    }

    return seed;
}
