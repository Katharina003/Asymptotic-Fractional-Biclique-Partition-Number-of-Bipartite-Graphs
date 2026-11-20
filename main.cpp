#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <map>
#include <memory>
#include <set>
#include <algorithm>
#include <utility>
#include <tuple>
#include "gurobi_c++.h"
#include <unordered_map>
#include <queue>
#include <boost/multiprecision/cpp_int.hpp> // Arbitrary-Precision Integer
#include <numeric>   // std::gcd
#include <stdexcept> // std::invalid_argument
#include <cstdlib>   // std::abs
#include <unordered_set>
#include <regex>
#include <fstream>
#include <bitset>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

#include "CondensedNode.h"
#include "EdgeOrbit.h"
#include "CondensedFormula.h"
#include "Biclique.h"
#include "DominoGraph.h"
#include "BicliqueParser.h"
#include "IsBiclique.h"
#include "QuboSmall.h"
#include "Printing.h"
#include "QuboEdgeorbits.h"
#include "QuboLinearized.h"
#include "Solver.h"


using namespace boost::multiprecision;
using namespace std;

vector< vector<vector<string>>> DCollection;

vector<vector<string>>& getDCollection(int power){
    vector<vector<string>> D = {{"n","b","a"},{"d","c","d"},{"a","b","n"}};
    if(DCollection.size()==0){
        DCollection.push_back(D);
    }
    while(DCollection.size()<=power){
        DCollection.push_back(kroneckerProduct(DCollection[DCollection.size()-1],D,true));
    }
    return DCollection[power];
}

string readFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }
    stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}
//------------------------------------------------MAIN----------------------------------------------------------
#pragma region MAIN

int main() {
    std::ofstream("mip2.log", std::ios::trunc).close();

    print(getDCollection(2), 3, false);
    print(getDCollection(1),2,true);

    while(true){
        string c;
        cout<< "This is the program mentioned in the master's thesis 'Asymptotic Fractional Biclique Partition Number of Bipartite Graphs'."<<endl;
        cout<< "It solves the aggregated LP for Kronecker powers of the Domino graph and can print some additional information."<<endl<<endl;
        cout<< "To solve an LP, press l."<<endl;
        cout<< "To highlight bicliques, press h."<<endl;
        cout<< "To print <a,b,c,d>, press m."<<endl;
        cout<< "To check whether a condensed formula (using only mirroring symmetry) is a biclique, press b."<<endl;
        cout<< "To print a Kronecker power of D, press p."<<endl;
        cout<< "To print the current dual edge values for the k-th power of <a,b,c,d>, press f."<<endl;
        cout<< "To exit, press q."<<endl;
        cin>>c;
        if(c=="h"){
            cout<<"Please enter the Kronecker power. Then enter the edges to highlight, separated by spaces. End with a dot separated from the last edge by a space."<<endl;
            vector<string> s;
            int k;
            cin >> k;
            while(true){
                string a;
                cin>>a;
                if(a==".")break;
                else s.push_back(a);
            }

            highlight(getDCollection(k-1), k, s, false);
        }
        else if (c == "f") {
            cout << "For which Kronecker power should the file be printed?" << endl;
            int k;
            cin >> k;
            string fileName = "qubo" + to_string(k) + ".txt";
            std::ifstream file(fileName);
            if (!file.is_open()) {
                cout << "Please start solving the LP for Kronecker power " << k
                     << " first to see the latest version of " << fileName << '.' << endl;
                continue;
            }
            std::string line;

            while (std::getline(file, line)) {
                std::cout << line << '\n';
            }

            file.close();
        }
        else if (c == "m") {
            cout << "For which Kronecker power should <a,b,c,d> be printed?"<<endl;
            int k;
            cin >> k;
            printABCD(getDCollection(k-1), k);
        }
        else if(c=="b"){
            cout << "Enter the Kronecker power k, the number of left nodes in the maximal biclique, "
                    "the number of right nodes in the maximal biclique, and the number of edge orbits in the condensed formula, separated by spaces.\n"
                    "Next, enter the left node labels, followed by the right node labels, separated by spaces. "
                    "On both sides, a node label is a k-bit binary string entered as a decimal integer "
                    "(for k = 2: 00 = 0, 01 = 1, 10 = 2, 11 = 3).\n"
                    "Finally, enter the edge orbits, separated by spaces. Each orbit is a string of exactly k "
                    "letters from a, b, c, d. Each letter describes the corresponding pair of left and right bits: "
                    "a = (0,0), b = (1,0), c = (1,1), d = (0,1). "
                    "Letter order matters: ab and ba describe different edge orbits. "
                    "For k = 2, examples are aa ab ba." << endl;
            int k,numberTopNodes,numberBottomNodes, numberEdgesSumformular;
            cin>>k>>numberTopNodes>>numberBottomNodes>>numberEdgesSumformular;
            vector<int> topNodesMaxClique(numberTopNodes);
            vector<int> bottomNodesMaxClique(numberBottomNodes);
            vector<string> sumformula(numberEdgesSumformular);
            for(int i=0;i<numberTopNodes;++i){
                cin>>topNodesMaxClique[i];
            }
            for(int i=0;i<numberBottomNodes;++i){
                cin>>bottomNodesMaxClique[i];
            }
            for(int i=0;i<numberEdgesSumformular;++i){
                string s;
                cin>>s;
                sumformula[i]=s;
            }
            vector<CondensedNode> topNodes, bottomNodes;
            for (int name : topNodesMaxClique) topNodes.push_back(CondensedNode(name, 1, k));
            for (int name : bottomNodesMaxClique) bottomNodes.push_back(CondensedNode(name, 1, k));
            vector<pair<CondensedNode, CondensedNode>> formula;
            for (const string& orbit : sumformula) {
                int left = 0, right = 0;
                for (char letter : orbit) {
                    //for each letter shift binary representation by one, and change last digit to 1
                    //if the node orbit is L_1^1 for left nodes (R_1^1 for right nodes respectively)
                    left = (left << 1) | (letter == 'b' || letter == 'c');
                    right = (right << 1) | (letter == 'c' || letter == 'd');
                }
                formula.push_back({CondensedNode(left, 1, k), CondensedNode(right, 1, k)});
            }
            cout << isBiclique(Biclique(topNodes, bottomNodes), formula, k);
            
        }
        else if (c == "l") {
            cout << "Select the pricing problem: 1 = QuboEdgeorbits, 2 = QuboSmall, 3 = QuboLinearized, 4 = EnumerateAll (can take a very long time for k > 3)." << endl;
            string pricingChoice;
            if (!(cin >> pricingChoice)) break;
            if (pricingChoice != "1" && pricingChoice != "2" && pricingChoice != "3" && pricingChoice != "4") {
                cout << "Invalid pricing problem. Returning to the main menu." << endl;
                continue;
            }
            int pricingType = stoi(pricingChoice);
            cout << "Enter the Kronecker power k for which D^k should be solved. For k > 3, do not expect fast results." << endl;
            int k;
            cin >> k;
            cout << "Which Kronecker power should be used as the starting point? "
                    "Enter 1 to start from the beginning, or an integer r with 2 <= r <= k "
                    "to load the saved solution for power r from resultKr.txt (e.g., resultK2.txt for r = 2). "
                    "Computation then starts at power r + 1. If r = k, the program loads the solution "
                    "for power k and attempts to improve it." << endl;
            int answer;
            cin >> answer;

            //lösung für k = 1:
            Biclique abcd = Biclique({ CondensedNode(0,1,1), CondensedNode(1,1,1) }, { CondensedNode(0,1,1), CondensedNode(1,1,1) });
            Biclique d2 = Biclique({ CondensedNode(0,2,1) }, { CondensedNode(1,1,1) });
            Biclique ab = Biclique({ CondensedNode(0,1,1),CondensedNode(1,1,1) }, { CondensedNode(0,1,1)});
            vector<Biclique> solD1 = { abcd,d2,ab };
            vector<Biclique> solKMinusOne = solD1;
            if (answer > 1) {
                string solutionText;
                try {
                    solutionText = readFile("resultK" + std::to_string(answer) + ".txt");
                }
                catch (const runtime_error& error) {
                    cout << error.what() << " Returning to the main menu." << endl;
                    continue;
                }
                solKMinusOne = BicliqueParser::parseSolutionBicliques(solutionText);
            }
            for (Biclique b : solKMinusOne) cout << b.to_string()<<endl;
            int startvalue = min(answer + 1,k);
            for (int i = startvalue;i <= k;++i) {
                const auto start = std::chrono::steady_clock::now();
                DominoGraph DI = DominoGraph(i);
                cout << "Initialization finished." << endl;
                if(answer<k)generateInitialSolution(DI, solD1, solKMinusOne);
                else {
                    for (Biclique b : solKMinusOne)DI.insert(b);
                }
                cout << DI.to_string() << endl;
                solKMinusOne = solveLPColumnGeneration(DI, pricingType);
                const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
                cout << "Time for k = " << i << " including certificate output: " << elapsed.count() << " seconds." << endl;
            }
            cout << "Done!" << endl<<endl<<endl;
        }
        else if (c == "p") {
            cout << "Please enter the Kronecker power k for which D^k should be printed." << endl;
            int k;
            cin >> k;
            if(k<= 4)print(getDCollection(k-1),k,true);
            else { cout << "This Kronecker power is too large to display in this console."<<endl; }
        }
        else if (c == "q")break;
    }
    return 0;
}

#pragma endregion
