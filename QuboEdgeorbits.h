#pragma once

#include "Biclique.h"
#include "DominoGraph.h"
#include "PricingProblem.h"
#include "gurobi_c++.h"

#include <utility>
#include <vector>

/*
class for compressing edge orbits in pricing problem as described in thesis
*/
class QuboEdgeorbits : public PricingProblem {
    double lowerBound=-GRB_INFINITY;
    double upperBound=GRB_INFINITY;
    bool isTimeOuted = false;
    bool memLimitReached = false;
    void setUpQubo(Biclique& maxBiclique, std::vector<double>& dualSolutions, DominoGraph& dominoGraph, GRBModel& qubo, std::vector<GRBVar>& x, int timeLimit);

public:
    int calculateBestSubbiclique(Biclique& maxBiclique, std::vector<double>& dualSolutions, DominoGraph& dominoGraph, std::vector<double>& previousSolution, int timeLimit) override;
    std::pair<double, double> getSolutionInterval() override;
    const bool isTimeOut() const override;
    const bool reachedMemLimit() const override;
};
