#pragma once

#include "CondensedNode.h"

#include <cstddef>
#include <functional>
#include <string>

/*
class representing edge orbits regarding A(D^k). Only consists of a string with the orbit name
and some extra functions for functionality belonging to edge orbits.
*/
class EdgeOrbit {
private:
    std::string orbitName;//name of form a^ib^jc^ld^h where ^x means repeat the letter x times.

public:
    EdgeOrbit(const std::string name);

    const std::string getName() const;

    static int factorial(int k);


    static int edgeOrbitSize(EdgeOrbit edge);

    static EdgeOrbit topBottomToEdge(CondensedNode topNode, CondensedNode botNode);
    bool operator<(const EdgeOrbit& other) const;
    bool operator==(const EdgeOrbit& other) const;
};

namespace std {
template<>
struct hash<EdgeOrbit> {
    std::size_t operator()(const EdgeOrbit& e) const;
};
}
