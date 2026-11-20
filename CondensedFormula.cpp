#include "CondensedFormula.h"

using namespace std;

/*
Looks up how many edges of edge orbit e are used in the condensed formula and 
returns this value
*/
int CondensedFormula::getNumOcc(const EdgeOrbit& e) const {
    if (edgeOrbitToNumOcc.count(e)>0) return edgeOrbitToNumOcc.at(e);
    else return 0;
}

/*
changes the condensed formula by increasing the count of edge orbit e by amount.
*/
void CondensedFormula::incrementEdgeCount(EdgeOrbit &e,int amount) {
    //inserts 0 if edge is not already part of list, before increasing by amount.
    edgeOrbitToNumOcc[e]+=amount; 
}

/*
sets the number of edges in edge orbit e to numOcc
*/
void CondensedFormula::setEdgeCount(EdgeOrbit& e, int numOcc) {
    edgeOrbitToNumOcc[e] = numOcc;
}

/*
Outputs a string of the format
(N: <orbit name n1>, #<num of edges in orbit n1>)
...
(N: <orbit name ni>, #<num of edges in orbit ni>)
without linebreaks. Only outputs edge orbits that where a some point added 
(and potentially removed at a later point) to the condensed formula.
*/
string CondensedFormula::to_string() const {
    string s;
    for (auto& [orbit, num] : edgeOrbitToNumOcc) {
        s += "(N: " + orbit.getName() + ", #" + std::to_string(num) + ") ";
    }
    return s;
}
