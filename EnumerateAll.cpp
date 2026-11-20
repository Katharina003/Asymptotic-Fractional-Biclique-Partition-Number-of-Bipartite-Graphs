#include "EnumerateAll.h"
#include "Biclique.h"
#include "DominoGraph.h"

using namespace std;


void enumerate(
    size_t index, //index selects which orbit is currently being inserted
    const vector<CondensedNode>& top, //of maximal biclique
    const vector<CondensedNode>& bottom, //of maximal biclique
    vector<CondensedNode>& subTop,
    vector<CondensedNode>& subBottom,
    DominoGraph& dominoGraph,
    int& added
) {
    //when all orbits have been considered the biclique has to be added.
    if (index == top.size() + bottom.size()) {
        //only non-empty bicliques allowed
        if (!subTop.empty() && !subBottom.empty()) {
            Biclique subbiclique(subTop, subBottom);
            if (dominoGraph.insert(subbiclique)) ++added;
        }
        return;
    }
    
    //check whether index currently points to a left or right node orbit and get the corresponding orbit of the maxibiclique
    const CondensedNode& node = index < top.size() ? top[index] : bottom[index - top.size()];
    //select the side of the subbiclique in which we therefore need to insert:
    vector<CondensedNode>& side = index < top.size() ? subTop : subBottom;

    //insert 0 to |[e]| many nodes of the node orbit and continue recursively
    //insert no node of the current orbit and move to the next orbit
    enumerate(index + 1, top, bottom, subTop, subBottom, dominoGraph, added);
    for (int count = 1; count <= node.getNumOccurences(); ++count) {
        //insert count nodes of the current orbit and move on
        side.push_back(CondensedNode(node.getName(), count, node.getKroneckerPower()));
        enumerate(index + 1, top, bottom, subTop, subBottom, dominoGraph, added);
        side.pop_back();
    }
}

/*
Add all bicliques to the domino graph. The primal problem will then find the best objective value.
*/
int EnumerateAll::calculateBestSubbiclique(Biclique& maxBiclique,
    vector<double>& ignore1, DominoGraph& dominoGraph, vector<double>& ignore2, int) {
    const auto top = maxBiclique.getTop();
    const auto bottom = maxBiclique.getBot();
    vector<CondensedNode> subTop, subBottom;
    int added = 0;
    //index is binary representation of which nodes get inserted
    enumerate(0, top, bottom, subTop, subBottom, dominoGraph, added);
    return added;
}
