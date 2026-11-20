#pragma once

#include "EdgeOrbit.h"

#include <string>
#include <unordered_map>

/*
Class to store condensed formulas
Stores a map that returns for any edge orbit, how often it appears in the condensed formula
*/
class CondensedFormula {
private:
    std::unordered_map<EdgeOrbit,int> edgeOrbitToNumOcc;
public:
    int getNumOcc(const EdgeOrbit& e) const;

    void incrementEdgeCount(EdgeOrbit &e,int amount = 1);

    void setEdgeCount(EdgeOrbit& e, int numOcc);

    std::string to_string() const;
};
