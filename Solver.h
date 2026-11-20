#pragma once

#include "BigInt.h"
#include "DominoGraph.h"
#include "gurobi_c++.h"

#include <vector>

BigInt determinant_bareiss(std::vector<std::vector<BigInt>> A);
BigInt absBigInt(const BigInt& x);
BigInt gcdBigInt(BigInt a, BigInt b);
void reduceFraction(BigInt& numerator, BigInt& denominator);

void generateInitialSolution(DominoGraph& dominoGraph,
                             std::vector<Biclique>& MBaseKMinusOne,
                             std::vector<Biclique>& MBaseOne);
std::vector<Biclique> evaluateSolution(
    GRBModel& model, 
    DominoGraph& dominoGraph, 
    bool calculateFractions, 
    std::vector<std::pair<double,double>>& gaps,
    double solutionValueForGaps);
std::vector<Biclique> solveLPColumnGeneration(DominoGraph& dominoGraph, int pricingType = 1);
