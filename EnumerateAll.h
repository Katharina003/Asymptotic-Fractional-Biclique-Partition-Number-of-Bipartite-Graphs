#pragma once

#include "PricingProblem.h"
#include <limits>

/*
class that inserts all subbicliques of a maximal biclique into the restricted master problem
*/
class EnumerateAll : public PricingProblem {
public:
    int calculateBestSubbiclique(Biclique& maxBiclique,
        std::vector<double>& dualSolutions, DominoGraph& dominoGraph,
        std::vector<double>& previousSolution, int timeLimit) override;

    // Enumeration does not evaluate the dual objective.
    std::pair<double, double> getSolutionInterval() override {
        return { 0.0, std::numeric_limits<double>::infinity() };
    }
    const bool isTimeOut() const override { return false; }
    const bool reachedMemLimit() const override { return false; }
};
