#pragma once

#include <utility>
#include <vector>

class Biclique;
class DominoGraph;

/*
Class all Pricing Problems implement for easy swapping of implementations
*/
class PricingProblem {
public:
    virtual ~PricingProblem() = default;

    // Adds subbicliques with the highest reduced cost to dominoGraph.
    // Returns the number of added bicliques.
    virtual int calculateBestSubbiclique(
        Biclique& maxBiclique,
        std::vector<double>& dualSolutions,
        DominoGraph& dominoGraph,
        std::vector<double>& previousSolution,
        int timeLimit) = 0;

    //lower and upper bound on solutions found when calculating best subbiclique
    //only call this when calculate best subbiclique has run for meaningful values
    virtual std::pair<double, double> getSolutionInterval() = 0;
    
    //information on the termination reasons of subproblem
    virtual const bool isTimeOut() const = 0;
    virtual const bool reachedMemLimit() const = 0;
};
