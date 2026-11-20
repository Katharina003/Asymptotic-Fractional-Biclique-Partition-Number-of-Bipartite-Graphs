#pragma once

#include "Biclique.h"

#include <string>
#include <vector>

/*
Class which parses resultK<KroneckerPower>.txt files for computing initial solutions and
for solution verification.
*/
class BicliqueParser {
private:

    static CondensedNode parseNode(const std::string& s);

    static std::vector<CondensedNode> parseNodeList(const std::string& s);

public:

    static Biclique parse(const std::string& input);

    static std::vector<Biclique> parseSolutionBicliques(const std::string& input);
};
