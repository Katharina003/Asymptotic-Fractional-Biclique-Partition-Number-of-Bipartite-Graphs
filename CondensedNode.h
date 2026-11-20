#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

/*
Class storing node orbits regarding mirroring symmetry.
Both left and right nodes are refered to by integers whose binary representation
shows whether the i-th coordinate of the node is in the L_0^1 or L_1^1 (R_0^1 or R_1^1 respectively)
*/
class CondensedNode {
private:
    int name;
    int numOccurences;
    int kroneckerPower;

public:
    CondensedNode(int newName,
                  int newNumOccurences,
                  int newKroneckerPower);

    const int getNumOccurences() const;
    const int getName() const;
    const int getKroneckerPower() const;

    void increaseNumOccurencesBy(int amount);

    static CondensedNode unsortedProductOfNodes(
        const CondensedNode& a,
        const CondensedNode& b
    );

    static CondensedNode sortedProductOfNodes(CondensedNode a, CondensedNode b);

    static void condenseNodeList(std::vector<CondensedNode>& nodeList);

    static std::vector<CondensedNode> nodeProductUnordered(
        const std::vector<CondensedNode>& A,
        const std::vector<CondensedNode>& B
    );

    bool operator==(const CondensedNode& other) const;
    bool operator<(const CondensedNode& other) const;

    std::string to_string() const;
};

namespace std {
template<>
struct hash<CondensedNode> {
    std::size_t operator()(const CondensedNode& n) const;
};
}
