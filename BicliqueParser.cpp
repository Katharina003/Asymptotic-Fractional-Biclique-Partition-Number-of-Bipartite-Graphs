#include "BicliqueParser.h"

#include <regex>
#include <stdexcept>

using namespace std;

/*
Parses one condensed node from a string of the format
(N: <name>, #<num occurrences>, K: <Kronecker power>)
where the name is a sequence of 0s and 1s
num occurrences is a non-negative integer and
Kronecker power is the Kronecker power of the Domino graph it belongs to.
E.g.
(N: 10, #3, K:2)

Note that the binary name string must be exactly as long as the Kronecker power specified in the node.
*/
CondensedNode BicliqueParser::parseNode(const string& s) {
    regex nodeRegex(R"(\(N\s*:\s*([01]+)\s*,\s*#(\d+)\s*,\s*K:\s*(\d+)\s*\))");
    smatch match;

    if (!regex_search(s, match, nodeRegex)) {
        throw runtime_error("Invalid node format: " + s);
    }

    string binaryName = match[1];
    int numOccurrences = stoi(match[2]);
    int kroneckerPower = stoi(match[3]);

    int name = 0;
    if (binaryName.size() != kroneckerPower) {
        throw runtime_error("Invalid node, Kronecker power "+to_string(kroneckerPower)
        +" does not match name length of node "+binaryName);
    }
    for (char c : binaryName) {
        name = (name << 1) + (c - '0');
    }
    return CondensedNode(name, numOccurrences, kroneckerPower);
}

/*
Parses a string representing a list of nodes (as they appear in Top Nodes and Bot Nodes).
Searches for all non-nested parenthesized blocks in the string and gives them to the node parser.
e.g.
(N: 1, #1, K: 1) ,;/  (N: 0, #1, K: 1)
*/
vector<CondensedNode> BicliqueParser::parseNodeList(const string& s) {
    vector<CondensedNode> result;

    regex nodeRegex(R"(\([^()]*\))");

    auto begin = sregex_iterator(s.begin(), s.end(), nodeRegex);
    auto end = sregex_iterator();

    for (auto it = begin; it != end; ++it) {
        result.push_back(parseNode(it->str()));
    }

    return result;
}

/*
Parses one biclique. Expects a string of the format:
Top Nodes: (...), ..., (...)
Bot Nodes: (...), ..., (...)
where (...), ..., (...) are lists of valid nodes, each on one line.
Commas between nodes are optional. The bottom list must end with a newline or '}'.
*/
Biclique BicliqueParser::parse(const string& input) {

    regex topRegex(R"(Top Nodes\s*:\s*([\s\S]*?)(?:\n|Bot Nodes:))");
    regex botRegex(R"(Bot Nodes\s*:\s*([\s\S]*?)(?:\n|\}))");

    smatch topMatch;
    smatch botMatch;

    if (!regex_search(input, topMatch, topRegex)) {
        throw runtime_error("Could not parse Top Nodes");
    }

    if (!regex_search(input, botMatch, botRegex)) {
        throw runtime_error("Could not parse Bot Nodes");
    }

    vector<CondensedNode> topNodes =
        parseNodeList(topMatch[1]);

    vector<CondensedNode> botNodes =
        parseNodeList(botMatch[1]);

    return Biclique(topNodes, botNodes);
}

/*
Parses the relevant information from a result file computed by this program:
It finds the part labeled "Solution Bicliques" and parses the curly-bracketed blocks as bicliques
until it hits the label "MBasis:" (or the end of the input if the label is missing).
It returns all solution bicliques as a vector.
*/
vector<Biclique> BicliqueParser::parseSolutionBicliques(const string& input) {
    //vector to store the result
    vector<Biclique> result;

    //find relevant substring
    size_t start = input.find("Solution Bicliques");
    if (start == string::npos) {
        throw runtime_error("No Solution Bicliques section found");
    }
    size_t end = input.find("MBasis:", start);
    if (end == string::npos) {
        end = input.size();
    }
    string sub = input.substr(start, end - start);

    //every biclique has non-nested curly brackets around it. Keep track of where
    //to search for the next biclique
    size_t pos = 0;
    while (true) {
        size_t open = sub.find('{', pos);
        if (open == string::npos) break; // when no more brackets are found before MBasis we can stop
        size_t close = sub.find('}', open + 1);

        //something is wrong if no closing bracket follows the opening bracket
        if (close == string::npos) {
            throw runtime_error("Unbalanced braces");
        }

        //select the string including the curly brackets and parse its biclique.
        string block = sub.substr(open, close - open + 1);
        result.push_back(parse(block));

        pos = close + 1;
    }

    return result;
}
