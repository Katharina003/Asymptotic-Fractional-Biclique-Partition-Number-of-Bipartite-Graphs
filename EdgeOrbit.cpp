#include "EdgeOrbit.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;

/*
set edge orbit to given name
*/
EdgeOrbit::EdgeOrbit(const string name) {
    orbitName = name;
    sort(orbitName.begin(),orbitName.end());
}

/*
return edge orbit name
*/
const string EdgeOrbit::getName() const {
    return orbitName;
}

/*
compute k! if it fits into int (throws errors if it does not fit)
works at least for k<12

needed for computing orbit size.
*/
int EdgeOrbit::factorial(int k) {
    if (k < 0) throw std::invalid_argument("k must be non-negative");
    int result = 1;
    for (int i = 2; i <= k; ++i) {
        if (result > std::numeric_limits<int>::max() / i) {
            throw std::overflow_error("factorial overflow");
        }
        result *= i;
    }
    return result;
}

/*
compute edge orbit size of given edge orbit.
*/
int EdgeOrbit::edgeOrbitSize(EdgeOrbit edge) {
    //evaluate which orbit is represented by edge name by counting the occurrences of a to d.
    int a = 0;
    int b = 0;
    int c = 0;
    int d = 0;
    int power = edge.getName().size();
    for (char i : edge.getName()) {
        if (i == 'a')++a;
        else if (i == 'b')++b;
        else if (i == 'c')++c;
        else if (i == 'd')++d;
        else {
            cout << "fehler, buchstabe " << i << " ist nicht erlaubt." << endl;
        }
    }

    //use formula of master's thesis applying the orbit-stabilizer theorem
    int groupSize = (1 << power) * factorial(power);
    int orbitSize = groupSize / (factorial(a) * factorial(b) * factorial(d) * (factorial(c) * (1 << c)));
    return orbitSize;
}

/*
given a left and a right node, compute for the edge connecting the two in which
orbit it lies.
As nodes are binary representations of the coordinate wise orbits they belong to,
we can apply the construction of the thesis to compute the edge orbit of the edge
*/
EdgeOrbit EdgeOrbit::topBottomToEdge(CondensedNode topNode, CondensedNode botNode) {
    string result = "";
    int power = topNode.getKroneckerPower();
    for (int i = 1; i < (1 << power);i *= 2) {
        //in the i-th coordinate the coordinates lie in L_1^1 and R_1^1
        if (((topNode.getName() & i) != 0) && ((botNode.getName() & i) != 0)) {
            result += "c";
        }
        //in the i-th coordinate the coordinates lie in L_1^1 and R_0^1
        if (((topNode.getName() & i) != 0) && ((botNode.getName() & i) == 0)) {
            result += "b";
        }
        //in the i-th coordinate the coordinates lie in L_0^1 and R_1^1
        if (((topNode.getName() & i) == 0) && ((botNode.getName() & i) != 0)) {
            result += "d";
        }
        //in the i-th coordinate the coordinates lie in L_0^1 and R_0^1
        if (((topNode.getName() & i) == 0) && ((botNode.getName() & i) == 0)) {
            result += "a";
        }
    }
    //edge orbit names are lexicographically sorted.
    sort(result.begin(), result.end());
    return EdgeOrbit(result);
}

/*
returns which edge orbit has the lexicographically smaller term.
This is good enough as the edge orbit names get sorted in the constructor and
cannot be changed anywhere else.
*/
bool EdgeOrbit::operator<(const EdgeOrbit& other) const {
    return orbitName < other.orbitName;
}

/*
returns whether two edge orbits have the same name.
This is good enough as the edge orbit names get sorted in the constructor and
cannot be changed anywhere else.
*/
bool EdgeOrbit::operator==(const EdgeOrbit& other) const {
    return orbitName == other.orbitName;
}

/*
hash hashes the only field of the class.
*/
std::size_t std::hash<EdgeOrbit>::operator()(const EdgeOrbit& e) const {
return std::hash<std::string>()(e.getName());
}
