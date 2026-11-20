#pragma once

#include "Biclique.h"
#include "EdgeOrbit.h"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>


/*
    class that stores all information about D^k and its corresponding
    restricted master problem of the aggregated LP
*/
class DominoGraph {
private:
    int power;
    int numEdgeOrbits; // number of different edge orbits in D^k using A(D^k)
    std::unordered_map<EdgeOrbit, int> edgesToIndices; //calculated by a deterministic function in the beginning. Mapping is the same for all runs with the same k.
    std::unordered_map<int, int> indicesToOrbitSizes; //orbit size of each edge orbit (referenced by index from function above)
    std::vector<Biclique> maxBicliques;//set of one representative per condensed formula class of maximal bicliques for subproblem.
    bool isMaxBicliquesSet = false;//variable used for determining whether max bicliques need to be computed or are already stored in maxBicliques
    std::unordered_set<Biclique> foundBicliques;//vector of all bicliques found by initial solution and column generation
    std::vector<std::vector<int>> M;//M of the restricted master problem of the aggregated LP
    std::unordered_map<int, Biclique> MtoBiclique;//mapping of columns of M (their index) to bicliques

    int computeNumEdgeOrbits();

    std::vector<Biclique>bicliqueProductTopBottomNodes(const std::vector<Biclique>& A, const std::vector<Biclique>& B);

    std::vector<Biclique> generateMaxBicliques(int k);

    void generateEdgesToIndices();

    void generateEdgesToOrbitSizes();

    int insertColumn(std::vector<int>& column);
public:
    DominoGraph(int k);

    std::vector<Biclique> getMaxBicliques();

    int getOrbitSize(int edgeIndex);

    int getPower();

    int getNumEdgeOrbits();

    int getNumBicliques();

    int getEdgeIndex(EdgeOrbit& edge);

    const int getM(int edgeOrbit, int biclique);

    Biclique getBiclique(int indexInM);

    std::vector<EdgeOrbit> getEdges();

    bool insert(Biclique& b);

    std::string M_to_string();

    std::string to_string();
};
