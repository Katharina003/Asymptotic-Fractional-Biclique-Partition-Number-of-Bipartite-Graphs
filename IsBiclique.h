#pragma once

#include "Biclique.h"
#include <utility>
#include <vector>


bool isBiclique(Biclique maxbiclique,
                std::vector<std::pair<CondensedNode,CondensedNode>> condensedFormula,//not using EdgeOrbit class, as it also uses Kronecker permutations
                int kroneckerPower);
