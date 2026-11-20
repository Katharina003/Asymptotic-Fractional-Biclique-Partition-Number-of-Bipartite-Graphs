#include "DominoGraph.h"

#include <iostream>
#include <stdexcept>

using namespace std;

/*
Counts the number of edge orbits in A(D^k) by summing 
(instead of using the formula from the thesis)
*/
int DominoGraph::computeNumEdgeOrbits() {
    int num = 0;
    for (int i = 0; i <= power;++i) {//number of a coordinates
        for (int j = 0; j <= power - i;++j) {//number of b coordinates
            num += (power - i - j) + 1; //number of c coodrdinates
            //the number of d coordinates is fixed afterwards
        }
    }
    return num;
}

/*
Takes two vectors of bicliques and computes for every pair in AxB their Kronecker product.

returns a new vector containing the Kronecker product of all possible pairs in AxB
*/
vector<Biclique>DominoGraph::bicliqueProductTopBottomNodes(const vector<Biclique>& A, const vector<Biclique>& B) {
    vector<Biclique> result;
    for (const Biclique B1: A) {
        for (const Biclique B2 : B) {
            vector<CondensedNode> newTop =
                CondensedNode::nodeProductUnordered(B1.getTop(), B2.getTop());
            vector<CondensedNode> newBottom =
                CondensedNode::nodeProductUnordered(B1.getBot(), B2.getBot());
            result.push_back(Biclique( newTop, newBottom ));
        }
    }
    return result;
}

/*
Function that generates all Maximal Bicliques of the k-th Kronecker Power of the Domino Graph.
The result is NOT stored in the field maxBicliques of the DominoGraph class.
*/
vector<Biclique> DominoGraph::generateMaxBicliques(int k) {
    //The three maximal bicliques of D
    Biclique abcd = Biclique({ CondensedNode(0,1,1), CondensedNode(1,1,1) }, { CondensedNode(0,1,1), CondensedNode(1,1,1) });
    Biclique b2c = Biclique({ CondensedNode(0,2,1),CondensedNode(1,1,1) }, { CondensedNode(1,1,1) });
    Biclique d2c = Biclique({ CondensedNode(1,1,1) }, { CondensedNode(0,2,1),CondensedNode(1,1,1) });
    vector<Biclique> DMax = { abcd,b2c,d2c };
    if (k == 1) return DMax;

    //Computing all maximal bicliques of D^k is the same as computing the maximal bicliques
    //of D^{k-1} and D and taking a Kronecker product of all their combinations.
    //therefore we can use the standard method of computing all maximal bicliques of
    //D^{floor(k/2)} and computing all combinations of bicliques in D^{floor(k/2)} x D^{floor(k/2)}
    //if 2xfloor(k/2)<k we then have to compute the product of the result with D once again.
    vector<Biclique> half = generateMaxBicliques(k / 2);
    vector<Biclique> result = bicliqueProductTopBottomNodes(half, half);
    if (k % 2 != 0) {
        result = bicliqueProductTopBottomNodes(result, DMax);
    }
    return result;
}

/*
for iterating over edge orbits it is helpful to generate an index for every edge orbit
Therefore this method computes the lexicographical order of edges to map them to indices

edge orbits are represented as a string a...ab...bc...cd...d, they are lexicographically ordered in themselves.
We use this property in our code.
For Kronecker power k=2 we compute the order 
aa -> ab -> ac -> ad -> bb -> bc -> bd -> cc -> cd -> dd 
*/
void DominoGraph::generateEdgesToIndices() {
    //start with string a...a
    string s(power, 'a');
    int index = 0;
    while (true) {
        // assign the current edge orbit the current index, increase the index afterwards by one.
        edgesToIndices[EdgeOrbit(s)] = index++; 

        // Find right most position that is not d.
        int i = power - 1;
        while (i >= 0 && s[i] == 'd') i--; 

        // If all positions are d we are done
        if (i < 0) break; 

        // otherwise we count the non-d position up by one, getting the lexicographically next letter.
        s[i]++; 
        // finally we reduce all indicies to the right to contain the same letter
        // to keep the invariant that the edges are lexicographically ordered and to 
        // continue with the next smallest orbit.
        for (int j = i + 1; j < power; ++j) s[j] = s[i];
    }
}

/*
we store for every edge how many edges lie in its orbit.
*/
void DominoGraph::generateEdgesToOrbitSizes() {
    for (auto& [edge,index] : edgesToIndices) {
        indicesToOrbitSizes[index] = EdgeOrbit::edgeOrbitSize(edge);
    }
}

/*
M is the Matrix of the restricted master problem of the aggregated LP of the master's thesis
This method allows us to insert condensed formulas = columns into M as we need for the column generation approach.
It is important that we always insert in the back. We use this invariant for consistency
between the Gurobi problem instance and our data structure.

The method returns the index at which the column was inserted
*/
int DominoGraph::insertColumn(vector<int>& column) {
    for (int i = 0;i < M.size();++i) { // i iterates over the edge orbits
        M[i].push_back(column[i]); // insert the count of the condensed formula for that edge orbit
    }
    return M[0].size() - 1;
}

/*
Constructor of Domino graph D^k. Takes power k and computes all edge orbits and their sizes
Additionally it creates matrix M of the restricted master problem of the aggregated LP
*/
DominoGraph::DominoGraph(int k) {
    power = k;
    numEdgeOrbits = computeNumEdgeOrbits();
    generateEdgesToIndices();
    generateEdgesToOrbitSizes();
    M = vector<vector<int>>(numEdgeOrbits);
}

/*
Returns the set of all representatives of maximal bicliques needed for the 
subproblems of the restricted master problem for the aggregated LP
*/
vector<Biclique> DominoGraph::getMaxBicliques() {
    if (!isMaxBicliquesSet) maxBicliques = generateMaxBicliques(power);
    isMaxBicliquesSet = true;
    return maxBicliques;
}

/*
looks up the orbit size of edgeorbit with edge orbit index
*/
int DominoGraph::getOrbitSize(int edgeIndex) {
    return indicesToOrbitSizes[edgeIndex];
}

/*
returns the Kronecker power
*/
int DominoGraph::getPower() {
    return power;
}

/*
returns the number of edge orbits of D^k
*/
int DominoGraph::getNumEdgeOrbits() {
    return numEdgeOrbits;
}

/*
returns the number of bicliques = columns currently part of the 
restricted master problem of the aggregated problem.
*/
int DominoGraph::getNumBicliques() {
    return foundBicliques.size();
}

/*
returns the index belonging to an edge orbit
*/
int DominoGraph::getEdgeIndex(EdgeOrbit& edge) {
    return edgesToIndices[edge];
}

/*
returns M_{[e],<B>}.
*/
const int DominoGraph::getM(int edgeOrbit, int biclique) {
    if (edgeOrbit >= M.size() || biclique >= M[0].size()) {
        throw  std::invalid_argument("zugriff auf M hat falsche dimensionen");
    }
    return M[edgeOrbit][biclique];
}

/*
returns the biclique which has the condensed formula written in the
indexInM-th column of M
*/
Biclique DominoGraph::getBiclique(int indexInM) {
    if(MtoBiclique.count(indexInM)==0)
        throw  std::invalid_argument("biclique nicht in M enthalten");
    return MtoBiclique.at(indexInM);
}

/*
returns a vector containing all edge orbits.
*/
vector<EdgeOrbit> DominoGraph::getEdges() {
    vector<EdgeOrbit> edges(getNumEdgeOrbits(),EdgeOrbit("null"));
    for (auto it = edgesToIndices.begin();it != edgesToIndices.end();++it) {
        edges[it->second] = it->first;
    }
    return edges;
}

/*
inserts a biclique into M if it is not already in there.

returns true, if the biclique was inserted,
false, if it was already part of the found bicliques.
*/
bool DominoGraph::insert(Biclique& b) {
    if (foundBicliques.count(b)>0) return false;

    //insert in found bicliques
    foundBicliques.insert(b);

    //insert in M
    vector<int> columnM(numEdgeOrbits, 0);
    CondensedFormula cf = b.getCondensedFormula();
    for (auto &[edge,index]:edgesToIndices) {
        columnM[index] += cf.getNumOcc(edge);
    }

    //keep mapping of columns of M to bicliques consistent
    int position = insertColumn(columnM);
    MtoBiclique.insert({position,b});
    return true;
}

/*
prints matrix M. Entries are separated by spaces, rows are separated by \n
*/
string DominoGraph::M_to_string() {
    string s;
    for (int i = 0;i < M.size();++i) {
        for (int j = 0;j < M[0].size();++j) {
            s += std::to_string(M[i][j]) + " ";
        } s += "\n";
    }
    return s;
}

/*
Prints all relevant information about the domino graph:
its power,
the number of edge orbits,
which bicliques are known to the restricted master problem (foundBicliques)
What does the corresponding M matrix look like
*/
string DominoGraph::to_string() {
    string s= "Domino Graph: {\n";
    s += "power: " + std::to_string(power) + "\n";
    s += "numEdgeOrbits: " + std::to_string(numEdgeOrbits) + "\n";
    s += "FoundBicliques: \n";
    for (auto b : foundBicliques) {
        s += b.to_string()+"\n";
    }
    s += M_to_string();
    s += "}";
    return s;
}
