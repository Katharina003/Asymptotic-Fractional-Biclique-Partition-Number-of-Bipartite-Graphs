#pragma once

#include "Biclique.h"
#include "DominoGraph.h"
#include "PricingProblem.h"
#include "gurobi_c++.h"

/*
This class linearizes the standard QUBO problem
It does NOT condense node orbits of any form, just adds aditional variables for linearization
*/
class QuboLinearized : public PricingProblem {
    double lowerBound = -GRB_INFINITY;
    double upperBound = GRB_INFINITY;
    bool isTimeOuted = false;
    bool memLimitReached = false;
    bool useObjThreshold = true;

    void setUpQubo(std::vector<double>& dualSolutions, DominoGraph& dominoGraph,
        GRBModel& qubo, std::vector<GRBVar>& xTop, std::vector<GRBVar>& xBot,
        std::vector<CondensedNode>& topExtended, std::vector<CondensedNode>& botExtended);

public:
    QuboLinearized() = default;
    QuboLinearized(bool useThreshold) : useObjThreshold(useThreshold) {}

    int calculateBestSubbiclique(Biclique& maxBiclique, std::vector<double>& dualSolutions,
        DominoGraph& dominoGraph, std::vector<double>& previousSolution, int timeLimit) override;
    std::pair<double, double> getSolutionInterval() override;
    const bool isTimeOut() const override;
    const bool reachedMemLimit() const override;
};
