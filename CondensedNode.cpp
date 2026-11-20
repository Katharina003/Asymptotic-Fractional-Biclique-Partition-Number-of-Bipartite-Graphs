#include "CondensedNode.h"

#include <algorithm>


/*
Construct a Condensed node by providing its orbit name n (regarding only mirroring symmetry,
using an integer value whos binary representation is a coordinate wise representation whether the 
i-th coordinate lies in L_0^1 or L_1^1 (R_0^1 or R_1^1 respectively)),
the number of nodes that lie in the orbit, 
and the Kronecker power k such that n in D^k
*/
CondensedNode::CondensedNode(int newName, int newNumOccurences, int newKroneckerPower) {
    name = newName;
    numOccurences = newNumOccurences;
    kroneckerPower = newKroneckerPower;

}

/*
returns the number of nodes in the orbit that are condensed by the condensed node
*/
const int CondensedNode::getNumOccurences() const {
    return numOccurences;
}

/*
return the name of the node. The name is an integer whos bi binary representation
shows whether the i-th coordinate of the node is in the L_0^1 or L_1^1 (R_0^1 or R_1^1 respectively)
*/
const int CondensedNode::getName() const {
    return name;
}

/*
returns the kronecker power the node belongs to.
(note that this is also the length of the binary representation of the name)
*/
const int CondensedNode::getKroneckerPower() const {
    return kroneckerPower;
}

/*
increases the number of occurrences the condensed node compresses by amount.
*/
void CondensedNode::increaseNumOccurencesBy(int amount) {
    numOccurences += amount;
}

/*
produces the product of two condensed nodes that does not use swapping of Kronecker factors
the new name appends the name of node a and node b to a new name ab

The condensed nodes contains one node per possible combinations of nodes in orbit a
with nodes in orbit b therefor the num of occurences is the product of their occurences

lastly the kronecker power is always the size of the two names.
Therefore we add the two Kronecker powers together.
*/
CondensedNode CondensedNode::unsortedProductOfNodes(const CondensedNode& a, const CondensedNode& b) {
    int newName = (a.getName() << b.getKroneckerPower()) + b.getName();
    int newNumOcc = a.getNumOccurences()* b.getNumOccurences();
    int newKroneckerPower = a.getKroneckerPower() + b.getKroneckerPower();
    return CondensedNode(newName, newNumOcc, newKroneckerPower);
}

/*
produces the product of two condensed nodes that DOES use swapping of Kronecker factors.
Be careful when to use it! If you want to generate edges at any point out of this condensed node, 
this will produce wrong results.
Also mixing this and the other function can lead to desaster!

The condensed node represents a subset of L_name^KroneckerPower (or R_name^KroneckerPower respectively)
Therefore we can construct the new name as in the thesis defined orbit product
by adding the two names together. In this representation the length of the name
is no longer a binary representation of the coordinate notation of the node,
but the index of node orbit L_i^k (or R_i^k respectively)

The condensed nodes contains one node per possible combinations of nodes in orbit a
with nodes in orbit b therefor the num of occurences is the product of their occurences

lastly the kronecker power is always the size of the two names.
Therefore we add the two Kronecker powers together.
*/
CondensedNode CondensedNode::sortedProductOfNodes(CondensedNode a, CondensedNode b) {
    int newName = a.getName() + b.getName();
    int newNumOcc = a.getNumOccurences() * b.getNumOccurences();
    int newKroneckerPower = a.getKroneckerPower() + b.getKroneckerPower(); //ist + hier richtig?
    return CondensedNode(newName, newNumOcc, newKroneckerPower);
}

/*
sorts nodeList and condenses same nodes into one node.
Does all changes to the given nodeList reference
*/
void CondensedNode::condenseNodeList(std::vector<CondensedNode>& nodeList) {
    std::sort(nodeList.begin(), nodeList.end());

    //add all nodes that share a name. Keep only one representative in which the sum is saved.
    for (std::size_t i = 1; i < nodeList.size(); ++i) {
        if (nodeList[i - 1].getName() == nodeList[i].getName()) {
            nodeList[i - 1].increaseNumOccurencesBy(nodeList[i].getNumOccurences());
            nodeList.erase(nodeList.begin() + i);
            --i;
        }
    }
}

/*
    takes two sets of condensed nodes (name, num of occurences).
    Concatinates the name of the first node in front of the name of the second node, for all combinations in the sets.
    Afterwards it collects nodes with the same name into one condensed node
*/
std::vector<CondensedNode> CondensedNode::nodeProductUnordered(const std::vector<CondensedNode>& A, const std::vector<CondensedNode>& B) {
    std::vector<CondensedNode> result;
    for (const CondensedNode& a : A) {
        for (const CondensedNode& b : B) {
            result.push_back(CondensedNode::unsortedProductOfNodes(a, b));
        }
    }
    CondensedNode::condenseNodeList(result);
    return result;
}

/*
Two nodes are the same, if 
their name match 
their number of occurences match and
their kronecker power matches
*/
bool CondensedNode::operator==(const CondensedNode& other) const {
    return name == other.getName() && numOccurences == other.getNumOccurences() && kroneckerPower == other.getKroneckerPower();
}

/*
We define order on the condensed nodes as follows:
We first order by kronecker Power,
then we order by name
lastly we order by num occurences
(this way it is easy to combine nodes of the same graph that share a name by adding their occurences,
as they lie next to each other in a sorted list)
*/
bool CondensedNode::operator<(const CondensedNode& other) const {
    if (kroneckerPower != other.kroneckerPower) return kroneckerPower < other.kroneckerPower;
    if (name != other.name) return name < other.name;
    return numOccurences < other.numOccurences;
}

/*
Condensed nodes are always printed in the following format:
(N: <name>, #<number of occurences>, K: <Kronecker power>)
the integer name is resolved in a bitstring with Kronecker power many digits:
for example:
(N: 1100, #2, K: 4) for name = 12
*/
std::string CondensedNode::to_string() const {
    int binary = name;
    std::string binaryName="";
    for (int i = 0;i < kroneckerPower;++i) {
        binaryName += std::to_string(binary % 2);
        binary = binary >> 1;
    }
    std::reverse(binaryName.begin(), binaryName.end());
    std::string s = "(N: "+binaryName + ", #" + std::to_string(numOccurences)+ ", K: "+std::to_string(kroneckerPower) + ")";
    return s;
}


/*
All three fields of the class are relevant for determining whether two elements are 
the same. Thus, the hash function uses all three for hashing.
*/
std::size_t std::hash<CondensedNode>::operator()(const CondensedNode& n) const {
    std::size_t seed = 0;

    auto hashCombine = [](std::size_t& seed, std::size_t value) {
        seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };
    hashCombine(seed, std::hash<int>()(n.getName()));
    hashCombine(seed, std::hash<int>()(n.getNumOccurences()));
    hashCombine(seed, std::hash<int>()(n.getKroneckerPower()));
    return seed;
}
