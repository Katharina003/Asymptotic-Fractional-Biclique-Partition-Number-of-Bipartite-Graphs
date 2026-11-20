#include "QuboLinearized.h"
#include "ObjThresholdCallback.h"

#include <iomanip>
#include <iostream>

using namespace std;


void QuboLinearized::setUpQubo(vector<double>& dualSolutions,
        DominoGraph& dominoGraph,
        GRBModel& qubo,
        vector<GRBVar>& xTop, // left nodes variables
        vector<GRBVar>& xBot, //right nodes variables
        vector<CondensedNode>& topExtended,//left nodes
        vector<CondensedNode>& botExtended //right nodes
    ) {

        xTop = vector<GRBVar>(topExtended.size());
        xBot = vector<GRBVar>(botExtended.size());

        for (int i = 0; i < topExtended.size(); ++i) {
            xTop[i] = qubo.addVar(0.0, 1.0, 0.0, GRB_BINARY, "xTop_" + to_string(i));
        }
        for (int i = 0; i < botExtended.size(); ++i) {
            xBot[i] = qubo.addVar(0.0, 1.0, 0.0, GRB_BINARY, "xBot_" + to_string(i));
        }
        vector<vector<GRBVar>> y = vector<vector<GRBVar>>(xTop.size(),(vector<GRBVar>(xBot.size())));

        //constraints
        for (int i = 0;i < xTop.size();i++) {
            for (int j = 0;j < xBot.size();j++) {
                y[i][j] = qubo.addVar(0, 1, 0, GRB_BINARY);

                qubo.addConstr(y[i][j] <= xTop[i]);
                qubo.addConstr(y[i][j] <= xBot[j]);
                qubo.addConstr(y[i][j] >= xTop[i] + xBot[j] - 1);
            }
        }

        //objective
        GRBLinExpr obj = 0;
        for (int i = 0; i < xTop.size(); ++i) {
            for (int j = 0; j < xBot.size(); ++j) {
                EdgeOrbit edge = EdgeOrbit::topBottomToEdge(topExtended[i], botExtended[j]);
                double value = dualSolutions[dominoGraph.getEdgeIndex(edge)];
                obj += value * y[i][j];
            }
        }
        qubo.setObjective(obj, GRB_MAXIMIZE);
    }

int QuboLinearized::calculateBestSubbiclique(Biclique& maxBiclique, vector<double>& dualSolutions, DominoGraph& dominoGraph, vector<double>& previousSolution, int timeLimit) {
    lowerBound = -GRB_INFINITY;
    upperBound = GRB_INFINITY;
    isTimeOuted = false;
    memLimitReached = false;
    try {
        auto top = maxBiclique.getTop();
        auto bottom = maxBiclique.getBot();
        GRBEnv environment;
        //environment.set(GRB_IntParam_OutputFlag, 0);
        environment.set("logFile", "mip2.log");
        GRBModel qubo = GRBModel(environment);
        ObjThresholdCallback cb;
        if (useObjThreshold) qubo.setCallback(&cb);
        qubo.set(GRB_DoubleParam_TimeLimit, timeLimit);
        qubo.set(GRB_DoubleParam_SoftMemLimit, 200.0);
        qubo.set(GRB_DoubleParam_MemLimit, 205.0);
        vector<GRBVar> xTop;
        vector<GRBVar> xBot;

        //we use only binary variables so condensed nodes need to be splitted in multiple nodes
        vector<CondensedNode> topExtended;
        for (const auto& node : top) {
            for (int i = 0;i < node.getNumOccurences();++i) topExtended.push_back(CondensedNode(node.getName(), 1, node.getKroneckerPower()));
        }
        vector<CondensedNode> botExtended;
        for (const auto& node : bottom) {
            for (int i = 0;i < node.getNumOccurences();++i) botExtended.emplace_back(node.getName(), 1, node.getKroneckerPower());
        }
        setUpQubo(dualSolutions, dominoGraph, qubo, xTop,xBot, topExtended, botExtended);

        if (previousSolution.size() == xTop.size() + xBot.size()) {
            for (int i = 0; i < xTop.size(); i++) {
                xTop[i].set(GRB_DoubleAttr_Start, previousSolution[i]);
            }
            for (int i = 0; i < xBot.size(); i++) {
                xBot[i].set(GRB_DoubleAttr_Start, previousSolution[i+xTop.size()]);
            }
        }
        else { previousSolution = vector<double>(xTop.size()+xBot.size()); }

        qubo.optimize();

        vector<int> solution(xTop.size()+xBot.size());

        //check Gurobis status.
        int status = qubo.get(GRB_IntAttr_Status);
        isTimeOuted = (status == GRB_TIME_LIMIT);
        memLimitReached = (status == GRB_MEM_LIMIT);
        try {
            upperBound = qubo.get(GRB_DoubleAttr_ObjBound);
        } catch (const GRBException& e) {
            if (e.getErrorCode() != GRB_ERROR_DATA_NOT_AVAILABLE) throw;
        }
        
        //statusflags are set, if no violated constraints are found, we can return
        if (qubo.get(GRB_IntAttr_SolCount) == 0) return 0;
        lowerBound = qubo.get(GRB_DoubleAttr_ObjVal);
        double optVal = qubo.get(GRB_DoubleAttr_ObjVal);
        cout << endl << endl << "optimum qubo: " << std::setprecision(10) << optVal << endl << endl;
        int numInsertions = 0;
        if (optVal > 1.0000001) { 
            int numSols = qubo.get(GRB_IntAttr_SolCount);
            cout << "numSols: " << numSols << endl;
            bool prevSolutionSet = false;
            for (int s = 0;s < numSols;++s) {
                qubo.set(GRB_IntParam_SolutionNumber, s);//selects the s-th solution in solution pool

                //only add benefitial bicliques to restricted master problem
                if (qubo.get(GRB_DoubleAttr_PoolObjVal) <= 1) {
                    continue;
                }
                cout << "qubo solution value: " << qubo.get(GRB_DoubleAttr_PoolObjVal) << endl;
                //read solution and round it up
                for (int i = 0; i < xTop.size(); ++i) {
                    solution[i] = static_cast<int>(xTop[i].get(GRB_DoubleAttr_Xn) + 0.5); 
                    if (!prevSolutionSet) previousSolution[i] = solution[i];
                }
                for (int i = 0; i < xBot.size(); ++i) {
                    solution[i+xTop.size()] = static_cast<int>(xBot[i].get(GRB_DoubleAttr_Xn) + 0.5);
                    if (!prevSolutionSet) previousSolution[i + xTop.size()] = solution[i + xTop.size()];
                }
                prevSolutionSet = true;

                //condense the binary nodes into condensed nodes again
                vector<CondensedNode> optTop;
                vector<CondensedNode> optBottom;
                for (int i = 0; i < topExtended.size(); ++i) {
                    if (solution[i] == 1) optTop.push_back(topExtended[i]);
                }
                //solution is one vector containing both top and then bottom nodes
                for (int i = 0; i < botExtended.size(); ++i) {
                    if (solution[topExtended.size() + i] == 1) optBottom.push_back(botExtended[i]);
                }
                Biclique opt(optTop, optBottom); //condenses nodes automatically
                if (dominoGraph.insert(opt)) ++numInsertions;
            }
            cout << "column generation done" << endl;
            return numInsertions;
        }
        else {
            // infeasible, unbounded, or interrupt
            cout << "no violated constraint found" << endl;
            for (int i = 0;i < xTop.size()+ xBot.size();++i)solution[i] = -1;

            return numInsertions;
        }

    }
    catch (const GRBException& e) {
        if (e.getErrorCode() != GRB_ERROR_OUT_OF_MEMORY) throw;
        memLimitReached = true;
        return 0;
    }
}

pair<double, double> QuboLinearized::getSolutionInterval() {
    return {lowerBound, upperBound};
}

const bool QuboLinearized::isTimeOut() const {
    return isTimeOuted;
}

const bool QuboLinearized::reachedMemLimit() const {
    return memLimitReached;
}
