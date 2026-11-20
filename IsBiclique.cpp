#include "IsBiclique.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

using namespace std;

/*
check for a given maximal biclique and condensed formula (using only mirroring symmetry),
whether it is possible to construct a biclique that has the condensed formula and is a subbiclique
of the maximal biclique.

The condensed formula is given as tuples of condensed nodes, because the condensed formula class uses
Kronecker permutations. In the tuple both numberOfOcc field should represent the number of edges, in
the orbit connecting nodes of .first orbit with .second orbit in the condensed formula.
The condensed nodes DO NOT contain the number of nodes (as we do not know them and want to derive them)
*/
bool isBiclique(Biclique maxbiclique,
                vector<pair<CondensedNode,CondensedNode>> condensedFormula,//not using EdgeOrbit class, as it also uses Kronecker permutations
                int kroneckerPower){
    if(condensedFormula.empty()){
        cout<<"empty biclique is not a valid subbiclique"<<endl;
        return false;
    }
    
     //sort edge lists to compress them
    sort(condensedFormula.begin(),condensedFormula.end());

    //buckets contains all compressed right nodes, leftNodes contains the name of the left nodes the buckets correspond to
    vector<vector<CondensedNode>> buckets;
    vector<CondensedNode> leftNodes;
    CondensedNode currentLeftNode=condensedFormula[0].first;
    vector<CondensedNode> rightNodesSameLeft;
    int numEdges = 0;
    for(auto [leftNode,rightNode] : condensedFormula){
        numEdges+= rightNode.getNumOccurences();
        if(leftNode.getName()==currentLeftNode.getName()){
            rightNodesSameLeft.push_back(rightNode);
        }
        else {
            CondensedNode::condenseNodeList(rightNodesSameLeft);
            buckets.push_back(rightNodesSameLeft);
            leftNodes.push_back(currentLeftNode);
            rightNodesSameLeft.clear();
            currentLeftNode=leftNode;
            rightNodesSameLeft.push_back(rightNode);
        }
    }
    CondensedNode::condenseNodeList(rightNodesSameLeft);
    buckets.push_back(rightNodesSameLeft);
    leftNodes.push_back(currentLeftNode);


    //check whether all left nodes use the same right nodes and compute scaling factors.
    //while doing so, also check whether scaling factor applies to all listelements
    vector<pair<int,int>> scalingFactorFraction(leftNodes.size());// numerator : denominator
    scalingFactorFraction[0]={1,1};
    for(size_t i=1;i<leftNodes.size();++i) {
        int scalingNumerator =  buckets[i][0].getNumOccurences();
        int scalingDenominator = buckets[0][0].getNumOccurences();
        //left nodes have different number of right node orbits
        if(buckets[0].size()!=buckets[i].size()) {
            cout<<"Left nodes do not agree on the node orbits of the right nodes. "
                "This is not a valid subset of any maximal biclique"<<endl;
            return false;
        }

        for(size_t list0Index = 0; list0Index<buckets[0].size();++list0Index) {
            //left nodes use different node orbits
            if(buckets[0][list0Index].getName() != buckets[i][list0Index].getName()){
                cout<<"Left nodes do not agree on the node orbits of the right nodes. "
                    "This is not a valid subset of any maximal biclique"<<endl;
                return false;
            }
            //scaling factors differ between right node orbits within two left nodes
            if(1LL * scalingNumerator * buckets[0][list0Index].getNumOccurences()
                != 1LL * scalingDenominator * buckets[i][list0Index].getNumOccurences()) {
                cout<< "there is no valid scaling factor for left node "<<leftNodes[0].to_string()
                << " and left node " <<leftNodes[i].to_string() <<". This is no valid subset of any maximal biclique"<<endl;
                return false;
            }
        }
        const int divisor = std::gcd(scalingNumerator, scalingDenominator);
        scalingFactorFraction[i]={scalingNumerator / divisor,scalingDenominator / divisor};
    }

    //compute LCM
    long long denominatorLcm = scalingFactorFraction.front().second; // we allready excluded empty bicliques
    for(size_t i = 1; i < scalingFactorFraction.size(); ++i) {
        denominatorLcm = std::lcm(denominatorLcm, static_cast<long long>(scalingFactorFraction[i].second));
        //The first left count is at least this LCM, and must divide every count in bucket zero.
        if (denominatorLcm > buckets[0][0].getNumOccurences()) return false;
    }

    //rescale scaling factor to integers to get smallest number of left nodes per orbit.
    //only multiples of this can be a possible number of nodes for a left orbit
    vector<long long> leastLeftNodes(leftNodes.size());
    for(size_t i = 0; i < scalingFactorFraction.size(); ++i) {
        leastLeftNodes[i] = scalingFactorFraction[i].first * (denominatorLcm / scalingFactorFraction[i].second);
    }

    //compute a multiple of the leastLeftNode vector that works for both left and right nodes.
    //could be improved by a binary search for allowed intervals (this does however not change runtime asymptotic runningtime): 
    //When there are more left nodes than in the maximal biclique
    //we need to decrease the multiple, if there are too many right nodes, we need to increase it.
    vector<CondensedNode> maxTopNodes = maxbiclique.getTop();
    vector<CondensedNode> maxBotNodes = maxbiclique.getBot();
    for(long long multiple = 1; multiple<= numEdges;++multiple){
        size_t maxTopPointer=0;
        bool possibleMultiple=true;
        //first check whether left nodes are subset of max biclique.
        // if too many are needed, this also holds for all coming
        //multiples so we are done and there is no valid subbiclique.
        for(size_t i=0;i<leftNodes.size();++i){
            while(maxTopPointer < maxTopNodes.size() && maxTopNodes[maxTopPointer].getName()!=leftNodes[i].getName()){
                maxTopPointer++;
            }
            if(maxTopPointer==maxTopNodes.size()) {
                cout<<"condensed formula is not subset of max biclique! No valid subset for this maximal biclique"<<endl;
                return false;
            }
            //to many left nodes needed.
            if(leastLeftNodes[i] > maxTopNodes[maxTopPointer].getNumOccurences() / multiple) {
                cout<<"not a valid subformula of maximal biclique"<<endl;
                possibleMultiple=false;
                return false;
            }
        }
        //next check whether all right nodes are subsets of max biclique and integer
        size_t maxBotPointer = 0;
        for(size_t j=0;j<buckets[0].size();++j){
            // see master's thesis where this fraction comes from.
            int numerator = buckets[0][j].getNumOccurences();
            long long denominator = multiple * denominatorLcm;

            //number of right nodes is a fraction:
            if(numerator % denominator != 0) {
                possibleMultiple=false;
                break;
            }
            int numRightNodes = numerator / denominator;

            //check whether right node is subset of maximal biclique
            while(maxBotPointer < maxBotNodes.size() && maxBotNodes[maxBotPointer].getName()!=buckets[0][j].getName()){
                maxBotPointer++;
            }

            //right nodes does not appear at all in maximal biclique:
            if(maxBotPointer==maxBotNodes.size()) {
                cout<<"condensed formula is not subset of max biclique! No valid subset for this maximal biclique"<<endl;
                return false;
            }

            //number of right nodes is not subset of max biclique:
            if(maxBotNodes[maxBotPointer].getNumOccurences() < numRightNodes) {
                possibleMultiple=false;
                break;
            }
        }

        //construct biclique by applying found multiple
        if (possibleMultiple) {
            //reconstruct left nodes
            vector<CondensedNode> bicliqueTop;
            for(size_t j=0;j<leftNodes.size();++j) {
                int numOcc = static_cast<int>(multiple*leastLeftNodes[j]);
                bicliqueTop.push_back(CondensedNode(leftNodes[j].getName(), numOcc, kroneckerPower));
            }
            //reconstruct right nodes
            vector<CondensedNode> bicliqueBot;
            for(size_t j=0;j<buckets[0].size();++j) {
                int numerator = buckets[0][j].getNumOccurences();
                long long denominator = multiple * denominatorLcm;
                int numOcc = numerator / denominator;
                bicliqueBot.push_back(CondensedNode(buckets[0][j].getName(), numOcc, kroneckerPower));
            }
            Biclique biclique(bicliqueTop,bicliqueBot);
            cout <<"valid Biclique!:" <<biclique.to_string()<<endl;
            return true;
        }
    }
    return false;
}
