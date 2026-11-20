#include "Solver.h"
#include "LowerBound.h"
#include "Printing.h"
#include "QuboEdgeorbits.h"
#include "QuboSmall.h"
#include "QuboLinearized.h"
#include "EnumerateAll.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <ctime>

using namespace std;

/*
 * Bareiss-Algorithm: fraction-free Gaussian elimination
 * Computes determinant exactly using integer arithmetic.
 */
BigInt determinant_bareiss(vector<vector<BigInt>> A) {
    int n = A.size();
    long long sign = 1;
    BigInt prev_pivot = 1;
    for (int k = 0; k < n - 1; ++k) {
        if (A[k][k] == 0) {
            int r = k + 1;
            while (r < n && A[r][k] == 0) r++;

            if (r == n) return 0; // whole column 0 -> singular -> det = 0

            swap(A[k], A[r]);     // swap rows
            sign = -sign;         // swaps flip sign of determinant
        }
        BigInt maxValue = 0;
        for (int i = k + 1; i < n; ++i) {
            for (int j = k + 1; j < n; ++j) {
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) / prev_pivot;
                maxValue = max(maxValue, A[i][j]);
            }
        }
        prev_pivot = A[k][k];
    }

    return sign * A[n - 1][n - 1];
}

BigInt absBigInt(const BigInt& x) {
    return x < 0 ? -x : x;
}

BigInt gcdBigInt(BigInt a, BigInt b) {
    while (b != 0) {
        BigInt r = a % b;
        a = b;
        b = r;
    }
    return a;
}

void reduceFraction(BigInt& numerator, BigInt& denominator) {
    if (denominator == 0) {
        throw std::invalid_argument("Denominator must not be zero.");
    }

    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }

    BigInt g = gcdBigInt(absBigInt(numerator), denominator);

    numerator /= g;
    denominator /= g;
}

/*
Input: Two integer matrices A and B
Output: Their Kronecker product
*/
static vector<vector<int>> kroneckerProduct(vector<vector<int>> &A, vector<vector<int>> &B) {
    vector<vector<int>> product(A.size() * B.size(), vector<int>( A[0].size() * B[0].size(),0));
    for (int i = 0;i < A.size();++i) {
        for (int j = 0;j < A[i].size();++j) {
            for (int a = 0;a < B.size();++a) {
                for (int b = 0;b < B[a].size();++b) {
                    product[B.size()*i+a][B[a].size()*j+b] = A[i][j] * B[a][b];
                }
            }
        }
    }
    return product;
}

/*
Iterates over all bicliques that are the best found solution of D^(k-1) and all bicliques of the optimal solution for D.
Computes the kroneckerproduct of every solution pair and adds it to the given dominoGraph.
In the thesis we explain why this is a valid first solution.
*/
void generateInitialSolution(DominoGraph& dominoGraph, vector<Biclique>& SolutionDKMinusOne, vector<Biclique>& SolutionDOne){
    //über alle paare von bicliquen iterieren und sie jeweils dem dominographen hinzufügen
    for (Biclique& b1 : SolutionDKMinusOne) {
        for (Biclique& b2 : SolutionDOne) {
            Biclique newBiclique = b1.kroneckerProduct(b2);
            dominoGraph.insert(newBiclique);
        }
    }
}

// Move one double towards the conservative side before decimal formatting.
// max_digits10 preserves that direction relative to the original double.
// This protects the output rounding, not the numerical accuracy of Gurobi.
static string roundedForOutput(double value, bool roundUp) {
    if (value >= GRB_INFINITY) return "inf";
    if (value <= -GRB_INFINITY) return "-inf";
    if (value != 0.0 && std::isfinite(value)) {
        value = std::nextafter(value, roundUp
            ? std::numeric_limits<double>::infinity()
            : -std::numeric_limits<double>::infinity());
    }
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return output.str();
}


/*
Computes the basis for the computed solution of Gurobi. It fills the resultK...txt file with the caluated exact values for the primal and dual solution.
The primal is a real upper bound of the aggregated LP, while the dual is based on the Gurobi outputs.
Gurobi uses floating point arithmetic, so the dual value could be bigger than the real lower bound.

The function returns a vector of all bicliques, that are part of the solution basis.
*/
vector<Biclique> evaluateSolution(
    GRBModel& model, //Master problem solved by Gurobi
    DominoGraph& dominoGraph, //DominoGraph class containing solution bicliques
    bool calculateFractions, //if this is true the exact values for x are calculated
    vector<pair<double,double>>& gaps, // the best found value as lower boudn and the best still possible value as upper bound for every maximal biclique
    double solutionValueForGaps // RMP objective before the pricing round that produced gaps
) {
    //printing current state of problem into resultFile
    const std::time_t now = std::time(nullptr);
    const std::tm localTime = *std::localtime(&now);

    std::ostringstream filename;
    filename << "resultK" << dominoGraph.getPower()
         << std::put_time(&localTime, "%Y-%m-%d_%H-%M-%S")
         << (calculateFractions ? "_final" : "")
         << ".txt";

    std::filesystem::create_directories("results");
    std::ofstream resultFile(std::filesystem::path("results") / filename.str());
    //std::ofstream resultFile("resultK"+to_string(dominoGraph.getPower())+".txt");
    resultFile << std::setprecision(10);
    resultFile << "result of optimization:" << endl;
    resultFile << "minimized value: " << roundedForOutput(model.get(GRB_DoubleAttr_ObjVal), true) << endl;
    resultFile << "solution value for gaps (rounded down): " << roundedForOutput(solutionValueForGaps, false) << endl;
    resultFile << "gaps: "<<endl;
    double worstUpperBound = 0;
    for(const auto& [lowerBound, upperBound] : gaps){
        worstUpperBound = max(worstUpperBound,upperBound);
        resultFile << "[" << roundedForOutput(lowerBound, false) <<", "<< roundedForOutput(upperBound, true) << "]"<<endl;
    }
    resultFile << "worst Upperbound: " << roundedForOutput(worstUpperBound, true) << endl;
    if (worstUpperBound < 1){
        resultFile << "no numerical lower bound calculated due to worst upperbound < 1." <<endl;
    } else if (std::isfinite(worstUpperBound) && worstUpperBound < GRB_INFINITY &&
        std::isfinite(solutionValueForGaps) && solutionValueForGaps >= 0.0) {
        const double numerator = solutionValueForGaps == 0.0 ? 0.0 : std::nextafter(
            solutionValueForGaps, -std::numeric_limits<double>::infinity());
        const double denominator = std::nextafter(worstUpperBound, std::numeric_limits<double>::infinity());
        resultFile << "numerical lower bound (rounded down): " << roundedForOutput(numerator / denominator, false) << endl;
    } else {
        resultFile << "numerical lower bound: unavailable (unknown or infinite bound)" << endl;
    }
    resultFile << "value of x:" << endl;
    int x_size = model.get(GRB_IntAttr_NumVars);
    std::unique_ptr<GRBVar[]> x(model.getVars());
    for (int i = 0;i < x_size;++i) {
        if (x[i].get(GRB_DoubleAttr_X) != 0) {
            resultFile << "x" << i << ": " << std::setprecision(10) << x[i].get(GRB_DoubleAttr_X)<<endl;
        }
    }

    vector<Biclique> bicliquesInBasis;
    //Basis Matrix:
    vector<vector<BigInt>> MBasis(dominoGraph.getNumEdgeOrbits(), vector<BigInt>(dominoGraph.getNumEdgeOrbits()));
    //indices i for which x_i is a basis variable / s_i is a slackvariable:
    vector<int> indexBasisX(dominoGraph.getNumEdgeOrbits());

    int numBasisVariables=0;
    int index = 0; //indicates which column of the MBasis matrix should be written
    resultFile << "Solution Bicliques:" << endl;

    //Filter all bicliques of the restricted master problem for basis bicliques
    for (int bicliqueIndex = 0; bicliqueIndex < x_size; ++bicliqueIndex) {
        //Gurobi-status vBasis: 0 if basis variable, -1 to -3 not basis variable
        if (x[bicliqueIndex].get(GRB_IntAttr_VBasis) == 0) {
            bicliquesInBasis.push_back(dominoGraph.getBiclique(bicliqueIndex));
            resultFile <<"biclique number "<<bicliqueIndex<<":"<< dominoGraph.getBiclique(bicliqueIndex).to_string() << endl;
            indexBasisX[index]=bicliqueIndex;
            for (int edgeorbit = 0;edgeorbit < dominoGraph.getNumEdgeOrbits();++edgeorbit) {
                MBasis[edgeorbit][index] = dominoGraph.getM(edgeorbit,bicliqueIndex);
            }
            index++;
            numBasisVariables++;
        }
    }

    //The basis could be smaller than the number of edges. 
    //Iterate over the slack variables (1 per edge orbit) and check whether it is basic by checking CBasis.
    //CBasis=0 means the i-th slackvariable (so the vector with all 0 and only a 1 in the i-th column) is basis.
    //Cbasis=-1 means it -th slackvariable is not basic.
    //add all basic slack variables
    std::unique_ptr<GRBConstr[]> constraints(model.getConstrs());
    for(int edgeOrbitIndex = 0; edgeOrbitIndex < dominoGraph.getNumEdgeOrbits(); ++edgeOrbitIndex){
        //Gurobi-status cBasis: 0 if basis variable, -1 if not
        if(constraints[edgeOrbitIndex].get(GRB_IntAttr_CBasis) == 0){
            if(index>=dominoGraph.getNumEdgeOrbits()){
                throw std::runtime_error("Basis size exceeds the number of edge orbits.");
            }
            MBasis[edgeOrbitIndex][index] = 1; // Mbasis is initialized as 0 matrix. Only need so set i-th entry to 1
            indexBasisX[index]=edgeOrbitIndex;
            ++index;
        }
    }

    //Check whether the MBasis is filled
    if(index!=dominoGraph.getNumEdgeOrbits()){
        throw std::runtime_error("Basis size does not correspond to number of edge orbits.");
    }

    //print MBasis and Domino Graph to result file
    resultFile << "MBasis: "<<MBasis.size() << endl;
    for (int i = 0;i < MBasis.size();++i) {
        for (int j = 0;j < MBasis[0].size();++j) {
            resultFile << MBasis[i][j] << " ";
        }resultFile << endl;
    }
    resultFile << "Dominograph: " << dominoGraph.to_string() << endl;

    //Calculate the exact fractions of primal solution x and dual solution y
    if (calculateFractions) {
        BigInt det = determinant_bareiss(MBasis);

        resultFile << "determinante von M: " << det << endl;

        //use Cramer's Rule (exchange every line by the sizes of the edge orbits once)
        
        //first only compute the determinants of the newly constructed matrices:
        vector<BigInt> det_I(MBasis.size()); //contains the det for the matrix where the i-th line is replaced by the edge orbit sizes
        //i is the line for which the edge orbit sizes are entered
        for (int i = 0;i < MBasis.size();++i) {
            vector<vector<BigInt>> MBasis_I(MBasis.size(), vector<BigInt>(MBasis.size()));
            for (int j = 0;j < MBasis.size();++j) {
                for (int k = 0;k < MBasis.size();++k) {
                    if (k == i) {
                        MBasis_I[j][k] = dominoGraph.getOrbitSize(j);
                    }
                    else {
                        MBasis_I[j][k] = MBasis[j][k];
                    }

                }
            }
            det_I[i] = determinant_bareiss(MBasis_I);
        }

        //then apply Cramer's Rule to compute the x values with det(M_i) / det(M).
        resultFile << "Solution values: " << endl;
        vector<pair<BigInt, BigInt>> fractions(det_I.size());
        for (int i = 0; i < det_I.size(); ++i) {
            BigInt numerator = det_I[i];
            BigInt denominator = det;
            reduceFraction(numerator, denominator);
            if(i < numBasisVariables){
                resultFile << "x_" << indexBasisX[i] << " = " << numerator << " / " << denominator <<endl;
            } else {
                resultFile << "s_" << indexBasisX[i] << " = " << numerator << " / " << denominator << endl;
            }
        }


        //------------------------------ dual values -------------------------
        //Calculate Dual solutions
        vector<int> orbitSizes(dominoGraph.getNumEdgeOrbits());
        for (int i = 0; i < dominoGraph.getNumEdgeOrbits(); ++i) {
            orbitSizes[i] = dominoGraph.getOrbitSize(i);
        }

        ExactDualCertificate dualCertificate =computeExactDual(MBasis, orbitSizes);

        resultFile << "Dual solution:" << endl;

        for (int i = 0; i < dualCertificate.integerDual.size(); ++i) {
            resultFile << "y_" << i << " = " << dualCertificate.integerDual[i]
                << " / " << dualCertificate.denominator << endl;
        }

        resultFile << "Restricted Master Problem Dual Lower Bound = " << dualCertificate.objectiveNumerator
            << " / " << dualCertificate.objectiveDenominator << endl;

        if(solveDualIntegerGap(dominoGraph,dualCertificate,0)) {
            resultFile<< "Dual lower bound verified"<<endl;
        }
    }

    resultFile.close();
    return bicliquesInBasis;
}


/*
Solves the restricted master problem as far as possible.
Returns the bicliques that form the best found solution.
*/
vector<Biclique> solveLPColumnGeneration(DominoGraph& dominoGraph, int pricingType) {
    //delete previous file contents of qubo-files
    std::ofstream file("qubo" + to_string(dominoGraph.getPower()) + ".txt", std::ios::trunc);
    file.close();
    //print size of initial solution and corresponding Matrix (submatrix of M in the aggregated LP)
    cout << "initial Solution:" <<dominoGraph.getNumEdgeOrbits() << " x " << dominoGraph.getNumBicliques() << endl;
    cout<<dominoGraph.M_to_string();
    
    //select the log environment for the restricted master problem
    GRBEnv env;
    env.set("LogFile", "mip1.log");
    GRBModel model = GRBModel(env);

    //add one continuous variable x_i with 0<=x_i<=1 per biclique (column) i of the restricted biclique to the model
    vector<GRBVar> x(dominoGraph.getNumBicliques());
    for (int i = 0; i < x.size(); ++i) {
        x[i] = model.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "x" + std::to_string(i + 1));
    }
    model.update();
    //Gurobi solves a minimization problem
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    //add constraints: for every edge orbit, the weights of the bicliques for that edge orbit should be exactly the size of the edge orbits
    //in the master's thesis we work with Mx=1 but multiplying each line with the size of the edge orbit yields
    //integer entries for M and the edge orbit sizes on the right side.
    vector<GRBConstr> constraints(dominoGraph.getNumEdgeOrbits());
    for (int edgeOrbit = 0; edgeOrbit < dominoGraph.getNumEdgeOrbits(); ++edgeOrbit) {
        GRBLinExpr expr = 0;
        for (int biclique = 0; biclique < x.size(); ++biclique) {
            expr += dominoGraph.getM(edgeOrbit,biclique) * x[biclique];
        }
        constraints[edgeOrbit] = model.addConstr(expr == dominoGraph.getOrbitSize(edgeOrbit), "row_" + std::to_string(edgeOrbit));
    }
    
    vector<pair<double,double>> gaps(dominoGraph.getMaxBicliques().size(),{-GRB_INFINITY,GRB_INFINITY});

    //solve the restricted master problem with the initial solution
    cout << "first optimization call" << endl;
    model.optimize();
    double solutionValueForGaps = model.get(GRB_DoubleAttr_ObjVal);
    evaluateSolution(model, dominoGraph, false,gaps, solutionValueForGaps);

    //previous solution stores for every maximal biclique how the best computed objective value was reached.
    vector<vector<double>> previousSolutions(dominoGraph.getMaxBicliques().size());

    //timelimit in seconds. Is increased whenever all pricing problems time out before 
    // finding violated constraints of computing that there are none.
    int timeLimit = 140;//in seconds
    bool isTimeouted = false; //true whenever at least one subproblem is timed out 
    bool everywhereOutOfMemory = true; //true whenever all subproblems reach memory limits / find no violated constraints and did not terminate early due to time limits.

    //looping stopps when in every biclique no violated constraint is found and no biclique is timed out
    //it also stopps when all subproblems terminate due to memory limits or do not find violated constraints
    while (true) {
        isTimeouted = false;
        everywhereOutOfMemory = true; 
        solutionValueForGaps = model.get(GRB_DoubleAttr_ObjVal);

        //get value of each edgeorbit from dual solution
        vector<double> dualSolution(dominoGraph.getNumEdgeOrbits());
        for (int i = 0; i < dominoGraph.getNumEdgeOrbits(); ++i) {
            dualSolution[i] = constraints[i].get(GRB_DoubleAttr_Pi);
        }
        
        int addedSubbiclique = 0; // number of different bicliques that violate constraints in any maximum biclique
        //solve subproblem for each maximal biclique
        for(int i=0;i<dominoGraph.getMaxBicliques().size();++i) {
            //print information into console for waiting users
            cout << endl;
            cout << "maximal biclique " << i << " of " << dominoGraph.getMaxBicliques().size() << endl;
            cout << "size of M: " << dominoGraph.getNumEdgeOrbits()<<" x " << dominoGraph.getNumBicliques() << endl;
            cout << "#addedSubbiclique until now: " << addedSubbiclique << endl;

            //<a,b,c,d> (maximal biclique on index i=0) gets extra treatment:
            //print interesting information in qubo file which can be read by user when pressing "f" in the start menu
            if (i == 0) {
                std::ofstream resultFile("qubo" + to_string(dominoGraph.getPower()) + ".txt", std::ios::app);
                resultFile<<printABCD(getDCollection(dominoGraph.getPower() - 1), dominoGraph.getPower(),dualSolution,dominoGraph);
                resultFile.close();
            }

            //generate pricing problem and solve it. calculateBestSubbicliques adds all bicliques 
            // that violated constraints to the dominoGraph and returns the number of added bicliques 
            // (this can be smaller than the number of violated constraints if any maximal biclique in a previous iteration allready added the biclique)
            std::unique_ptr<PricingProblem> qubo;
            if (pricingType == 2) qubo = std::make_unique<QuboSmall>();
            else if (pricingType == 3) qubo = std::make_unique<QuboLinearized>();
            else if (pricingType == 4) qubo = std::make_unique<EnumerateAll>();
            else qubo = std::make_unique<QuboEdgeorbits>();
            int numNewColumns = qubo->calculateBestSubbiclique(dominoGraph.getMaxBicliques()[i], dualSolution, dominoGraph,previousSolutions[i],timeLimit);
            cout << numNewColumns << " newColumn added."<<endl;

            //update information that is needed for termination
            isTimeouted = isTimeouted || qubo->isTimeOut();
            //whenever every biclique either does not find violated constraints or runs out of memory we do not improve by running another round of subproblems
            everywhereOutOfMemory = everywhereOutOfMemory && (qubo->reachedMemLimit() || (numNewColumns==0 && !qubo->isTimeOut()));
            cout << "adding new column" << endl;
            addedSubbiclique += numNewColumns;

            //update information about gap
            gaps[i] = qubo->getSolutionInterval();

            //add new columns to LP, by iterating over M and adding all columns from left to right, that are new.
            for(int biclique=dominoGraph.getNumBicliques() - numNewColumns;biclique <dominoGraph.getNumBicliques();++biclique) {
                GRBColumn col;
                for (int edgeorbit = 0;edgeorbit < dominoGraph.getNumEdgeOrbits();++edgeorbit) {
                    if (dominoGraph.getM(edgeorbit,biclique) != 0) {
                        col.addTerm(dominoGraph.getM(edgeorbit, biclique), constraints[edgeorbit]);
                    }
                }
                //add corresponding variable
                GRBVar newVar = model.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, col, "x" + std::to_string(x.size() + 1));
                x.push_back(newVar);
            }
        }
        //We are done, if no biclique was added and no time out occured.
        cout << "addedSubbiclique " << addedSubbiclique << endl;
        if (addedSubbiclique == 0 && !isTimeouted) break;
        //If we fail due to memory limits in every case it does not help to keep going
        else if (addedSubbiclique == 0 && everywhereOutOfMemory) break;
        //If no bicliques were added due to timeout, the time limit is increased.
        else if (addedSubbiclique == 0 && isTimeouted) {
            timeLimit *= 2;
        }

        //solve the restricted master problem.
        model.update();
        model.optimize();

        //print value to console and evaluate the new found solution of restricted master problem
        cout<<model.get(GRB_DoubleAttr_ObjVal) << endl;
        evaluateSolution(model, dominoGraph, false, gaps, solutionValueForGaps);
    }
    if(model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        throw std::runtime_error("restricted master problem is not solved until optimality");
    }
    return evaluateSolution(model, dominoGraph, true, gaps, solutionValueForGaps);
}
