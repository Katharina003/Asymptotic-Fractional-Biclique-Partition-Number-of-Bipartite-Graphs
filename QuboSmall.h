#pragma once

#include "Biclique.h"
#include "DominoGraph.h"
#include "PricingProblem.h"
#include "gurobi_c++.h"

#include <utility>
#include <vector>

/*
Solves pricing problem by compressing node orbits as described in thesis
*/
class QuboSmall : public PricingProblem {
    double lowerBound = -GRB_INFINITY;
    double upperBound = GRB_INFINITY;
    bool isTimeOuted = false;
    bool memLimitReached = false;
    bool useObjThreashold = true;

    void setUpQubo(Biclique& maxBiclique, std::vector<double>& dualSolutions, DominoGraph& dominoGraph, GRBModel& qubo, std::vector<GRBVar>& x, int timeLimit);

public:

    QuboSmall();

    QuboSmall(bool useThreashold);

    /*
    calculates the subbicliques with the highest reduced cost and adds them to the dominoGraph.
    Return the number of added bicliques.
    */
    int calculateBestSubbiclique(Biclique &maxBiclique, std::vector<double> &dualSolutions,DominoGraph& dominoGraph, std::vector<double> &previousSolution, int timeLimit) override;

    std::pair<double, double> getSolutionInterval() override;

    const bool isTimeOut() const override;
    
    const bool reachedMemLimit() const override;
};
