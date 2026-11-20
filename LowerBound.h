#pragma once

#include "BigInt.h"

#include <vector>

class DominoGraph;

struct ExactDualCertificate {
    BigInt determinant;

    // integerDual[i] = dual * l
    std::vector<BigInt> integerDual;
    BigInt denominator;

    // Exact lower bound = objectiveNumerator / objectiveDenominator
    BigInt objectiveNumerator;
    BigInt objectiveDenominator;
};

ExactDualCertificate computeExactDual(
    const std::vector<std::vector<BigInt>>& basis,
    const std::vector<int>& orbitSizes
);

double isTooLarge(const BigInt& x);
bool solveDualIntegerGap(
    DominoGraph& dominoGraph,
    const ExactDualCertificate& dual,
    double timeLimit = 0.0
);
