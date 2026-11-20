#include "LowerBound.h"
#include "Solver.h"

#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;

/*
given the basis of the restriced master problem,
compute the exact dual solution of the restricted master problem.
*/
ExactDualCertificate computeExactDual(
    const vector<vector<BigInt>>& basis,
    const vector<int>& orbitSizes
) {
    const int numEdgeOrbits = static_cast<int>(basis.size());

    BigInt detA = determinant_bareiss(basis);

    //as explained in masters thesis, replace ith row with the all ones vector
    //to apply Cramer's rule to get the exact dual values of the restricted master problem
    vector<BigInt> dualNumerators(numEdgeOrbits);
    vector<BigInt> dualDenominators(numEdgeOrbits, detA);
    for (int i = 0; i < numEdgeOrbits; ++i) {
        vector<vector<BigInt>> Adual = basis;
        for (int j = 0; j < numEdgeOrbits; ++j) {
            Adual[i][j] = 1;
        }
        dualNumerators[i] = determinant_bareiss(Adual);
        reduceFraction(dualNumerators[i], dualDenominators[i]);
    }

    //compute least common multiple l of the denominators as described in master's thesis
    BigInt l = 1;
    for (const BigInt& dualDenominator : dualDenominators) {
        l = (l / gcdBigInt(l, dualDenominator)) * dualDenominator;
    }

    // Multiply each dual value by l to obtain the integer vector l * y.
    vector<BigInt> integerDual(numEdgeOrbits);
    for (int i = 0; i < numEdgeOrbits; ++i) {
        integerDual[i] = (dualNumerators[i] * l) / dualDenominators[i];
    }

    //Compute Exact objective value. As described in masters thesis, the objective
    //of the dual problem using N instead of M is sum_[e] (y_[e] |[e]|)
    // when we extend y by multiplying it by L on the top and bottom we get the following
    BigInt objectiveNumerator = 0;

    for (int i = 0; i < numEdgeOrbits; ++i) {
        objectiveNumerator += BigInt(orbitSizes[i]) * integerDual[i];
    }
    BigInt objectiveDenominator = l;
    reduceFraction(objectiveNumerator, objectiveDenominator);

    //return everything
    ExactDualCertificate result;
    result.determinant = detA;
    result.integerDual = integerDual;
    result.denominator = l;
    result.objectiveNumerator = objectiveNumerator;
    result.objectiveDenominator = objectiveDenominator;

    return result;
}

/*
Gurobi uses doubles, allowing to represent integers in range [-2^53,2^53]
This method of "verifying" duals does not work for LPs that use
larger values. This method checks for a given integer, whether it is to 
large.
It returns true if if the number x is too large for exact double representation
false otherwise
*/
double isTooLarge(const BigInt& x)
{
    static const BigInt MAX_EXACT = (BigInt(1) << 53) - 1;
    return absBigInt(x) > MAX_EXACT;
}


/*
checks for a given dual certificate, if the dual can be verified.
if the dual is to large for double it returns false.
it also returns false, if the max{N^Tly}>l+epsilon
*/
bool solveDualIntegerGap(
    DominoGraph& dominoGraph,
    const ExactDualCertificate& dual,
    double timeLimit)
{
    if(isTooLarge( dual.denominator)) return false;

    const double qDouble =  dual.denominator.convert_to<double>();

    GRBEnv env;
    env.set(GRB_StringParam_LogFile, "mip3.log");
    env.set(GRB_IntParam_LogToConsole, 0);

    bool allFeasible = true;

    for (const Biclique& maxBiclique :
        dominoGraph.getMaxBicliques())
    {
        vector<CondensedNode> topNodes;
        vector<CondensedNode> bottomNodes;

        // Convert condensed nodes to binary nodes.
        for (const CondensedNode& n : maxBiclique.getTop()) {
            for (int j = 0; j < n.getNumOccurences(); ++j) {
                topNodes.push_back(CondensedNode( n.getName(),1,n.getKroneckerPower()));
            }
        }

        for (const CondensedNode& n : maxBiclique.getBot()) {
            for (int j = 0; j < n.getNumOccurences(); ++j) {
                bottomNodes.push_back(CondensedNode(n.getName(),1,n.getKroneckerPower()));
            }
        }

        GRBModel model(env);

        // add binary variables for all nodes
        vector<GRBVar> yTop(topNodes.size());
        vector<GRBVar> yBottom(bottomNodes.size());
        for (int i = 0; i < yTop.size(); ++i) {
            yTop[i] = model.addVar(0.0, 1.0, 0.0,GRB_BINARY,"top_" + to_string(i));
        }
        for (int j = 0; j < yBottom.size(); ++j) {
            yBottom[j] = model.addVar(0.0, 1.0, 0.0,GRB_BINARY,"bottom_" + to_string(j));
        }

        //we want to compute max{N^T*l*y}
        //for this we need to check for all subbicliques, what their solution value is
        //construct matrix that contains the scaled dual values as entries
        GRBQuadExpr objective = 0;
        for (int i = 0; i < yTop.size(); ++i) {
            for (int j = 0; j < yBottom.size(); ++j) {
                //get edge between two nodes and assign dual value
                EdgeOrbit edge = EdgeOrbit::topBottomToEdge(topNodes[i],bottomNodes[j]);
                int edgeIndex = dominoGraph.getEdgeIndex(edge);
                if(isTooLarge(dual.integerDual[edgeIndex])) return false;
                double weight = dual.integerDual[edgeIndex].convert_to<double>();
                objective += weight * yTop[i] * yBottom[j];
            }
        }

        model.setObjective(objective, GRB_MAXIMIZE);

        model.set(GRB_DoubleParam_MIPGap, 0.0);
        model.set(GRB_DoubleParam_MIPGapAbs, 0.0);

        if (timeLimit > 0.0) {
            model.set(
                GRB_DoubleParam_TimeLimit,
                timeLimit
            );
        }

        model.optimize();

        int status =
            model.get(GRB_IntAttr_Status);

        if (status == GRB_OPTIMAL) {
            double value = model.get(GRB_DoubleAttr_ObjVal);
            if (value > qDouble + 1e-7) {
                cout
                << "max weighted subbiclique = "
                << value
                << ", scaling factor = "
                << dual.denominator
                << endl;
            cout<< "dual is valid if the max weighted subbiclique is <= scaling factor +epsilon"<<endl;
                cout
                    << "Dual not feasible!"
                    << endl;
                return false;
            }
        }
        else if (status == GRB_TIME_LIMIT) {

            double bound =
                model.get(GRB_DoubleAttr_ObjBound);

            cout<< "TIME LIMIT, objective bound = "<< bound
            << ", scaling factor = " << dual.denominator << endl;

            if (bound > qDouble + 1e-7) {
                cout<< "could not verify dual in time!"<< endl;
                return false;
            }
        }
        else {
            cout << "Error during compution"<<endl;
            return false;
        }
    }

    return true;
}
