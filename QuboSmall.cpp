#include "QuboSmall.h"

#include "ObjThresholdCallback.h"
#include "Printing.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

/*
add constraints and objective value
*/
void QuboSmall::setUpQubo(Biclique& maxBiclique, vector<double>& dualSolutions, DominoGraph& dominoGraph, GRBModel& qubo, vector<GRBVar>& x, int timeLimit) {
    vector<CondensedNode> top = maxBiclique.getTop();
    vector<CondensedNode> bot = maxBiclique.getBot();

    x = vector<GRBVar>(top.size() + bot.size());
    bool isABCD = true;
    for (int i = 0; i < top.size(); ++i) {
        if (top[i].getNumOccurences() > 1)isABCD = false;
        x[i] = qubo.addVar(0.0, top[i].getNumOccurences(), 0.0, GRB_INTEGER, "x_" + to_string(i));
    }
    for (int i = top.size(); i < top.size() + bot.size(); ++i) {
        if (bot[i - top.size()].getNumOccurences() > 1)isABCD = false;
        x[i] = qubo.addVar(0.0, bot[i - top.size()].getNumOccurences(), 0.0, GRB_INTEGER, "x_" + to_string(i));
    }

    GRBQuadExpr obj = 0;

    // upper triangular matrix (i < j)
    for (int i = 0; i < top.size(); ++i) {
        for (int j = top.size(); j < x.size(); ++j) {
            // i is row index and j is columnindex
            EdgeOrbit edgeOrbit = EdgeOrbit::topBottomToEdge(top[i], bot[j - top.size()]);
            double value = dualSolutions[dominoGraph.getEdgeIndex(edgeOrbit)];
            obj += value * x[i] * x[j];
        }
    }
    qubo.setObjective(obj, GRB_MAXIMIZE);
    qubo.set(GRB_DoubleParam_TimeLimit, timeLimit);
    if (!isABCD) return;
    //Code for orderings. Uncomment if needed.
    /*
    vector<int> lexKleinsterIndex; //das i-te element ist der index, der am kleinsten ist für k/2 - i einsen
    for (int anzahlEinsen = dominoGraph.getPower() / 2;anzahlEinsen >= 0;--anzahlEinsen) {
        int kleinsterIndex = 0;
        for (int i = 0;i < anzahlEinsen;++i) {
            kleinsterIndex += 1 << i;
        }
        cout << "lex kleinester Index: " << top[kleinsterIndex].to_string() << endl;
        lexKleinsterIndex.push_back(kleinsterIndex);
        //über alle orbits iterieren (nicht besonders zeiteffizient aber da k klein relativ egal...)
        for (int i = 0; i < 1 << (dominoGraph.getPower());++i) {
            int einsen = 0;
            for (int j = 0;j < (dominoGraph.getPower());++j) {
                if ((i & (1 << j)) != 0)einsen++;
            }
            if (einsen == anzahlEinsen) {
                GRBLinExpr sum = 0;
                for (int j = 0; j < lexKleinsterIndex.size(); ++j) {
                    sum += x[lexKleinsterIndex[j]];
                }
                cout << "einfügen " << top[i].to_string() << endl;
                qubo.addConstr(
                    sum >= x[i],
                    "lex_order_" + std::to_string(i) + std::to_string(anzahlEinsen)
                );
            }
        }
    }*/
    /*
    //ordnung auf knoten mit genau hälfte (abgerundet) einsen
    int haelfteEinsenIndex = 0;
    for (int i = 0;i < (dominoGraph.getPower()) / 2;++i) {
        haelfteEinsenIndex += 1 << i;
    }

    for (int i = 0; i < 1 << (dominoGraph.getPower());++i) {
        int anzEinsen = 0;
        for (int j = 0;j < (dominoGraph.getPower());++j) {
            if ((i & (1<<j)) != 0)anzEinsen++;
        }
        if (anzEinsen == (dominoGraph.getPower()) / 2){
            cout << "einfügen " << top[i].to_string()<<endl;
            qubo.addConstr(
                x[haelfteEinsenIndex] >= x[i],
                "lex_order_" + std::to_string(i)
            );
        }
    }*/
    /*
    //ordnung auf knoten mit genau 1 eins
    for (int i = 0; i < dominoGraph.getPower() - 1; ++i) {
        int idx1 = 1 << i;
        int idx2 = 1 << (i + 1);

        qubo.addConstr(
            x[idx1] >= x[idx2],
            "lex_order_one" + std::to_string(i)
        );
    }


    //ordnung auf knoten mit genau 1 null
    for (int i = 0; i < dominoGraph.getPower() - 1; ++i) {
        int idx1 = ((1 << dominoGraph.getPower()) - 1) ^ (1 << i);
        int idx2 = ((1 << dominoGraph.getPower()) - 1) ^ (1 << (i + 1));

        qubo.addConstr(
            x[idx1]+x[1] >= x[idx2], //wenn das orbit mit genau 1 eins nicht verwendet wird
            "lex_order_zero" + std::to_string(i)
        );
    }

    //ordnung auf knoten mit genau 2 einsen
    int eins1 = 0;
    int eins2 = 1;
    while (dominoGraph.getPower() >= 4 && eins2 < dominoGraph.getPower()) {
        int idx1 = (1 << eins1)+(1<<eins2);
        for (int i = eins1; i < dominoGraph.getPower() - 1; ++i) {
            for (int j = i+1; j < dominoGraph.getPower(); ++j) {
                int idx2 = (1 << i) + (1 << j);
                if (idx1 == idx2) continue;
                qubo.addConstr(
                    x[idx1]
                    +x[1]
                    +x[((1 << dominoGraph.getPower()) - 1) ^ 1]
                    >= x[idx2], //wenn x[1]=0 gibt es nichts in orbit 0...01

                    "lex_order_one_one_" + std::to_string(eins1) + "_" + std::to_string(i) + "_" + to_string(j)
                );
            }
        }
        eins1 += 2;
        eins2 += 2;
    }

    //ordnung auf knoten mit genau 2 nullen
    int null1 = 0;
    int null2 = 1;
    while (dominoGraph.getPower() >= 5 && null2 < dominoGraph.getPower()) {
        int idx1 = ((1 << dominoGraph.getPower()) - 1) ^ ((1 << null1) + (1 << null2));
        for (int i = null1; i < dominoGraph.getPower() - 1; ++i) {
            for (int j = i+1; j < dominoGraph.getPower(); ++j) {
                int idx2 = ((1 << dominoGraph.getPower()) - 1) ^ ((1 << i) + (1 << j));
                if (idx1 == idx2) continue;
                qubo.addConstr(
                    x[idx1]
                    +x[1]
                    +x[((1 << dominoGraph.getPower()) - 1) ^ 1]
                    +x[3] 
                    >= x[idx2],
                    "lex_order_zero_zero_" + std::to_string(null1) + "_" + std::to_string(i) + "_" + std::to_string(j)
                );
            }
        }
        null1 += 2;
        null2 += 2;
    }*/
}

QuboSmall::QuboSmall() {

}

QuboSmall::QuboSmall(bool useThreashold) {
    useObjThreashold = useThreashold;
}

/*
calculate best subbiclique using the approach of node compression of thesis
*/
int QuboSmall::calculateBestSubbiclique(Biclique &maxBiclique, vector<double> &dualSolutions,DominoGraph& dominoGraph, vector<double> &previousSolution, int timeLimit) {
    try {
    vector<CondensedNode> top = maxBiclique.getTop();
    vector<CondensedNode> bot = maxBiclique.getBot();
    GRBEnv environment;
    environment.set("logFile", "mip2.log");
    GRBModel qubo = GRBModel(environment);
    qubo.set(GRB_DoubleParam_SoftMemLimit, 200.0);
    qubo.set(GRB_DoubleParam_MemLimit, 205.0);
    ObjThresholdCallback cb;//this has to be outside of if. Else, undefined behavior might occur.
    if (useObjThreashold) {
        qubo.setCallback(&cb);
    }
    vector<GRBVar> x;
    setUpQubo(maxBiclique, dualSolutions,dominoGraph, qubo, x, timeLimit);
    int numTopNodes = 0;
    for (CondensedNode node : top)numTopNodes += node.getNumOccurences();

    cout << "done setting up. Starting qubo" << endl;

    if(previousSolution.size()==x.size()){
    for (int i = 0; i < x.size(); i++) {
        x[i].set(GRB_DoubleAttr_Start, previousSolution[i]);
    }
    }
    else { previousSolution = vector<double>(x.size()); }
    int numInsertions = 0;
    GRBConstr variableCount;
    bool variableCountSet = false;
    for (int numVariables = 1; numVariables <= numTopNodes;++numVariables) {
        // Commented code can be used to iterate over all sums of nodes, as described in thesis extra approach
        // remove next line and uncomment commented part if you want to try.
        numVariables = numTopNodes;
        /*if (variableCountSet) {
            qubo.remove(variableCount);
            qubo.update();
        }
        GRBLinExpr sum = 0;
        for (int i = 0; i < top.size(); ++i) {
            sum += x[i];
        }
        variableCount = qubo.addConstr(
            sum == numVariables,
            "count"
        );
        qubo.update();*/

        variableCountSet = true;
        std::ofstream log("mip2.log", std::ios::app);

        log << "\n=== START iteration "
            << maxBiclique.to_string()<<"\n"
            << " : sum(x) == "
            << numVariables
            << " ===\n";
        cout << "\n=== START iteration "
            << maxBiclique.to_string() << "\n"
            << " : sum(x) == "
            << numVariables
            << " ===\n"<<endl;
        log.flush();
        qubo.optimize();
        vector<int> solution(x.size());

        int status = qubo.get(GRB_IntAttr_Status);
        if (status == GRB_TIME_LIMIT) {
            isTimeOuted = true;
        }
        else if (status == GRB_MEM_LIMIT){
            memLimitReached = true;
        }
        try {
            upperBound = qubo.get(GRB_DoubleAttr_ObjBound);   // best possible bound
        } catch (const GRBException& e) {
            if (e.getErrorCode() != GRB_ERROR_DATA_NOT_AVAILABLE) throw;
            upperBound = GRB_INFINITY;
        }
        if(qubo.get(GRB_IntAttr_SolCount) == 0) return 0;
        double optVal = qubo.get(GRB_DoubleAttr_ObjVal);
        lowerBound = qubo.get(GRB_DoubleAttr_ObjVal);     // best feasible column
        cout << endl << endl << "optimum qubo: " << std::setprecision(10) << optVal << endl << endl;
        vector<vector<int>> columns;
        if (optVal > 1.0000001) { 
            int numSols = qubo.get(GRB_IntAttr_SolCount);
            cout << "numSols: " << numSols << endl;
            bool prevSolutionSet = false;
            for (int s = 0;s < numSols;++s) {
                qubo.set(GRB_IntParam_SolutionNumber, s);
                if (qubo.get(GRB_DoubleAttr_PoolObjVal) <= 1) {
                    continue;
                }
                cout << "qubo solution value: " << qubo.get(GRB_DoubleAttr_PoolObjVal) << endl;


                for (int i = 0; i < x.size(); ++i) {
                    solution[i] = static_cast<int>(x[i].get(GRB_DoubleAttr_Xn) + 0.5); //round up
                    if (!prevSolutionSet) previousSolution[i] = solution[i];
                }
                prevSolutionSet = true;

                vector<CondensedNode> optTop;
                vector<CondensedNode> optBottom;
                bool isABCD = true;
                for (int i = 0;i < top.size();++i) {
                    if (solution[i] >= 1) {
                        optTop.push_back(CondensedNode(top[i].getName(), ceil(solution[i]), top[i].getKroneckerPower()));
                        if (optTop[optTop.size()-1].getNumOccurences() > 1)isABCD = false;
                    }
                }
                for (int i = 0;i < bot.size();++i) {
                    int solutionIndex = top.size() + i; //shift to get the correct elements
                    if (solution[solutionIndex] >= 1) {
                        optBottom.push_back(CondensedNode(bot[i].getName(), ceil(solution[solutionIndex]), bot[i].getKroneckerPower()));
                        if (optBottom[optBottom.size()-1].getNumOccurences() > 1)isABCD = false;
                    }
                }
                Biclique opt = Biclique(optTop, optBottom);
                cout << opt.to_string() << endl;
                cout <<"isABCD: " << isABCD << endl;
                if (isABCD) {
                    std::ofstream resultFile("qubo" + to_string(dominoGraph.getPower()) + ".txt", std::ios::app);
                    resultFile << "optvalue = " << qubo.get(GRB_DoubleAttr_PoolObjVal)<<"\n";
                    resultFile << opt.to_string() << "\n";
                    resultFile << printABCD(getDCollection(dominoGraph.getPower() - 1), dominoGraph.getPower(), dualSolutions, opt, dominoGraph);
                    resultFile.close();
                }
                cout << "opt: " << opt.to_string() << endl;
                if (dominoGraph.insert(opt))numInsertions++;
            }
            cout << "column generation done for "<<numVariables<<" top variables" << endl;

        }
        else {
            // infeasible, unbounded, or interupt
            cout << "no violated constraints found for "<<numVariables<< "top variables" << endl;
        }
    }
    return numInsertions;
    }
    catch (GRBException& e) {
        if(e.getErrorCode() == GRB_ERROR_OUT_OF_MEMORY){
            memLimitReached = true;
            return 0;
        } else {
            throw;
        }
    }
}

pair<double, double> QuboSmall::getSolutionInterval() {
    return { lowerBound,upperBound };
}

const bool QuboSmall::isTimeOut() const {
    return isTimeOuted;
}


const bool QuboSmall::reachedMemLimit() const {
    return memLimitReached;
}
