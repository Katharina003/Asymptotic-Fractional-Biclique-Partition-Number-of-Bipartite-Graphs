#pragma once

#include "Biclique.h"
#include "DominoGraph.h"

#include <string>
#include <tuple>
#include <utility>
#include <vector>

extern std::vector<std::vector<std::vector<std::string>>> DCollection;
std::vector<std::vector<std::string>>& getDCollection(int power);

struct RGB {
    int r, g, b;
};

std::string pairToString(std::pair<std::string, int> p);
RGB colorOfChar(char c);
RGB mixColors(const std::string& s);
std::string latexColorString(const std::string& s);
std::string printToLatex(const std::vector<std::vector<std::string>>& S);

std::vector<std::vector<std::string>> kroneckerProduct(
    std::vector<std::vector<std::string>>& D,
    std::vector<std::vector<std::string>>& P, bool dosort);

std::string colorize(const std::string& s, int potenz, bool doColorize);
std::vector<std::string> generateTernarySequence(int potenz);
std::string grayscale(const std::string& s, int potenz, bool doColorize);
void print(std::vector<std::vector<std::string>>& S, int potenz, bool doColorize);

std::string generateColorHighlight(
    const std::string& s, const std::vector<std::string>& highlights,
    const std::vector<std::tuple<int, int, int>>& colors, bool strong = false);
void highlight(std::vector<std::vector<std::string>>& S, int potenz,
               const std::vector<std::string>& highlightstrings, bool highlightMiddleOnes);

bool validRowIndex(int x);
bool validColIndex(int x);
std::string binaryLabel(int x, int k);
void printGraph(const std::vector<std::vector<std::string>>& graph,
                DominoGraph dominoGraph, std::vector<double> dual);

std::string printABCD(const std::vector<std::vector<std::string>>& S, int potenz);
std::string printABCD(const std::vector<std::vector<std::string>>& S, int potenz,
                      std::vector<double> dualSolution, DominoGraph& dominoGraph);
std::string printABCD(const std::vector<std::vector<std::string>>& S, int potenz,
                      std::vector<double> dualSolution, Biclique biclique,
                      DominoGraph& dominoGraph);
