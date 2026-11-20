#include "Printing.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

/*
A lot of code for graphical outputs.
*/

string pairToString(pair<string,int> p){
    string s = "("+p.first+", "+to_string(p.second)+")";
    return s;
}

RGB colorOfChar(char c) {
    switch(c) {
        case 'a': return {255, 0,   0};   // Red
        case 'b': return {0,   0, 255};   // Blue
        case 'c': return {0, 255,   0};   // Green
        case 'd': return {0,0,255};       // Blue
    }
    return {0,0,0};
}

RGB mixColors(const string &s) {
    int R=0, G=0, B=0;
    if (s.empty()) return {0,0,0};

    for (char c : s) {
        RGB col = colorOfChar(c);
        R += col.r;
        G += col.g;
        B += col.b;
    }
    int n = s.size();
    return {R/n, G/n, B/n};
}

string latexColorString(const string &s) {
    RGB c = mixColors(s);
    return "\\color[RGB]{" +
           to_string(c.r) + "," +
           to_string(c.g) + "," +
           to_string(c.b) + "}";
}

string printToLatex(const vector<vector<string>> &S) {
    string out;
    out += "\\begin{pmatrix}\n";

    for (int i = 0; i < S.size(); ++i) {
        for (int j = 0; j < S[i].size(); ++j) {

            if (S[i][j] == "n") {
                out += "\\cdot";
            } else {
                string val = S[i][j];
                string col = latexColorString(val);
                out += col + "{(" + val + ")}";
            }

            if (j + 1 < S[i].size()) out += " & ";
        }
        if (i + 1 < S.size()) out += " \\\\\n";
    }

    out += "\n\\end{pmatrix}";
    return out;
}

vector<vector<string>> kroneckerProduct(vector<vector<string>> &D, vector<vector<string>> &P,bool dosort){
    vector<vector<string>> product(D.size()*P.size(),vector<string> (D[0].size() * P[0].size()));
    for(int i=0;i<D.size();++i){
        for(int j=0;j<D[i].size();++j) {
            if(D[i][j]=="n"){
                for(int k=0;k<P.size();++k){
                    for(int h=0;h<P[k].size();++h) {
                        product[i*P.size()+k][j*P.size()+h] = "n";
                    }
                }
            }
            else{
                for(int k=0;k<P.size();++k){
                    for(int h=0;h<P[k].size();++h) {
                        if(P[k][h]=="n")product[i*P.size()+k][j*P.size()+h]="n";
                        else{
                            product[i*P.size()+k][j*P.size()+h] = D[i][j]+P[k][h];
                            if(dosort)sort(product[i*P.size()+k][j*P.size()+h].begin(),product[i*P.size()+k][j*P.size()+h].end());
                        }
                    }
                }
            }
        }
    }
    return product;
}

string colorize(const string &s, int potenz,bool doColorize) {
    if(!doColorize) return s;
    int r = 0, g = 0, b = 0;
    for (char ch : s) {
        if (ch == 'a') r += 250/potenz;
        else if (ch == 'b'|| ch=='d') b += 250 / potenz;
        else if (ch == 'c') g += 250/potenz;
        else if (ch == 'd') {
            r += 250 / potenz;
            b += 250 / potenz;
            g += 250 / potenz;
        }
    }
    // every color at most 255
    r = min(r, 255);
    g = min(g, 255);
    b = min(b, 255);

    // ANSI Escape for RGB-colors
    return "\033[38;2;" + to_string(r) + ";" + to_string(g) + ";" + to_string(b) + "m" + s + "\033[0m";
}

// generate all orbit names by generating ternary strings with 2 ->0
vector<string> generateTernarySequence(int potenz) {
    vector<int> digits(potenz, 0);
    long total = 1;
    for (int i = 0; i < potenz; ++i) total *= 3;

    vector<string> seq;
    seq.reserve(total);

    bool forward = true;

    for (long count = 0; count < total; ++count) {
        string s;
        for (int d : digits) {
            int out = (d == 2 ? 0 : d);
            s.push_back(char('0' + out));
        }
        seq.push_back(s);

        int pos = potenz - 1;
        while (pos >= 0) {
            if (digits[pos] < 2) {
                digits[pos]++;
                break;
            } else {
                digits[pos] = 0;
                pos--;
            }
        }

    }

    return seq;
}

string grayscale(const string &s, int potenz, bool doColorize) {
    if(!doColorize) return s;
    // count ones
    int ones = count(s.begin(), s.end(), '1');

    // Map 0 to power ones -> 0 to 255 grayscale
    int intensity = (200 * (ones+1)) / (potenz+1);
    return "\033[38;2;" + to_string(intensity) + ";" + to_string(intensity) + ";" + to_string(intensity)
           + "m" + s + "\033[0m";
}

void print(vector<vector<string>> &S, int potenz, bool doColorize) {
    vector<string> seq = generateTernarySequence(potenz);

    // node orbits (columns)
    for(int i=0;i<potenz;++i)cout<<" ";cout<< " ";
    for (auto &s : seq) cout << grayscale(s,potenz,doColorize) << " ";
    cout << "\n";

    for (int k = 0; k < S.size(); ++k) {
        cout << grayscale(seq[k],potenz,doColorize) << " ";    // node orbits (row lables)
        for (int h = 0; h < S[k].size(); ++h) {
            if (S[k][h] != "n")
                cout << colorize(S[k][h], potenz,doColorize) << " ";
            else {
                for (int i = 0; i < potenz; ++i) cout << ".";
                cout << " ";
            }
        }
        cout << "\n";
    }
}

string generateColorHighlight(const string &s,
    const vector<string> &highlights,
    const vector<tuple<int,int,int>> &colors,
    bool strong)
{
    for (size_t i = 0; i < highlights.size(); ++i) {
        if (s == highlights[i]) {
            auto [r, g, b] = colors[i];

            if (!strong) {
                r = 100;
                g = 100;
                b = 100;
            }
            return "\033[38;2;" + to_string(r) + ";" + to_string(g) + ";" + to_string(b) + "m"
                 + s + "\033[0m";
        }
    }
    return s;
}


void highlight(vector<vector<string>> &S, int potenz, const vector<string> &highlightstrings,bool highlightMiddleOnes) {
    vector<string> seq = generateTernarySequence(potenz);

   vector<tuple<int,int,int>> colors;
    vector<tuple<int,int,int>> base = {
        {255,   0,   0}, // rot
        {0,   255,   0}, // grün
        {0,     0, 255}, // blau
        {255, 255,   0}, // gelb
        {255,   0, 255}, // magenta
        {0,   255, 255}  // cyan
    };

    for (int i = 0; i < highlightstrings.size(); ++i) {
        colors.push_back(base[i % base.size()]);
    }

    // column labels
    for(int i=0;i<potenz;++i) cout << " ";
    cout << " ";

    for (size_t h = 0; h < (seq.size()); ++h) {
    bool strong = (!highlightMiddleOnes) || (highlightMiddleOnes && (h % 3 == 1));
    cout << generateColorHighlight(seq[h], highlightstrings, colors, strong) << " ";
    }
    cout << "\n";

    // rows
    for (int k = 0; k < (int)S.size(); ++k) {

        // row lables
        string rowLabel = seq[k];
        cout << generateColorHighlight(rowLabel, highlightstrings, colors) << " ";

        for (int h = 0; h < (int)(S[k].size()); ++h) {
            if (S[k][h] != "n") {
                bool strong = (!highlightMiddleOnes) || (highlightMiddleOnes&&(h%3==1)); //k%3==1 ||
                cout << generateColorHighlight(S[k][h], highlightstrings, colors,strong) << " ";
            } else {
                for (int i = 0; i < potenz; ++i) cout << ".";
                cout << " ";
            }
        }
        cout << "\n";
    }
}

bool validRowIndex(int x)
{
    // allowed: {0,1}
    while (x > 0)
    {
        int d = x % 3;

        if (d == 2)
            return false;

        x /= 3;
    }

    return true;
}

bool validColIndex(int x)
{
    // allowed: {1,2}
    while (x > 0)
    {
        int d = x % 3;

        if (d == 0)
            return false;

        x /= 3;
    }

    return true;
}

string binaryLabel(int x, int k)
{
    string s;

    for (int i = k - 1; i >= 0; --i)
    {
        s += ((x >> i) & 1) ? '1' : '0';
    }

    return s;
}



void printGraph(const std::vector<std::vector<std::string>>& graph, DominoGraph dominoGraph, vector<double> dual) {
    for (const auto& row : graph) {
        for (string s : row) {
            if (s[0] == 'n') {
                cout << "O ";
                continue;
            }
            EdgeOrbit edge = EdgeOrbit(s);
            int idx = dominoGraph.getEdgeIndex(edge);
            double v = dual[idx];
            double t = std::min(max(std::abs(v),0.8), 1.0);
            int intensity = static_cast<int>(t * 255);

            if (v >= 0) {
                // green
                std::cout << "\033[38;2;0;" + std::to_string(intensity) + ";0m" << "O" << "\033[0m ";
            }
            else {
                // red
                cout<< "\033[38;2;" + std::to_string(intensity) + ";0;0m" << "O" << "\033[0m ";
            }
        }
        std::cout << endl;
    }
}

string printABCD(const vector<vector<string>>& S, int potenz)
{
    int visibleSize = 1 << potenz; // 2^k

    //column labels
    string s;
    for (int i = 0; i < potenz; ++i) s+= " ";
    s+= " ";

    for (int i = visibleSize-1; i >= 0; --i) {
        s+= binaryLabel(i, potenz) + " ";
    }

    s+= "\n";

    //rows
    int rowLabel = 0;
    for (int r = 0; r < S.size(); ++r) {
        if (!validRowIndex(r)) continue;

        s+= binaryLabel(rowLabel, potenz) + " ";
        rowLabel++;

        for (int c = (int)((S[r].size() / 2) - 1); c < S[r].size(); ++c) {
            if (!validColIndex(c))continue;

            if (S[r][c] != "n") {
                s+= colorize(S[r][c], potenz, true) + " ";
            }
            else {
                for (int i = 0; i < potenz; ++i) s+= ".";
                s+= " ";
            }
        }
        s+= "\n";
    }
    cout <<s;
    return s;
}

string printABCD(const vector<vector<string>>& S, int potenz, vector<double> dualSolution, DominoGraph& dominoGraph) {
    string sBefore = printABCD(S, potenz);
    string s;
    int visibleSize = 1 << potenz; // 2^k
    //dual

    // column labels
    for (int i = 0; i < potenz; ++i) s+= " ";
    s+= " ";
    for (int i = visibleSize-1; i >= 0; --i) {
        s+= binaryLabel(i, potenz) + " ";
    }
    s+= "\n";

    //rows
    int rowLabel = 0;
    for (int r = 0; r < S.size(); ++r)
    {
        if (!validRowIndex(r))
            continue;

        s+= binaryLabel(rowLabel, potenz) + " ";
        rowLabel++;

        for (int c = (int)((S[r].size() / 2) - 1); c < S[r].size(); ++c) {
            if (!validColIndex(c))
                continue;

            if (S[r][c] != "n") {
                EdgeOrbit e(S[r][c]);
                double dual = dualSolution[dominoGraph.getEdgeIndex(e)];
                string value = to_string(dual);
                if (value[0] == '-')value = value.substr(1, dominoGraph.getPower());
                else value = value.substr(0, dominoGraph.getPower());
                if (dual < 0) {
                    value = "\033[38;2;240;0;0m" + value + "\033[0m";
                }
                else if (dual == 0) {
                    value = "\033[38;2;240;240;240m" + value + "\033[0m";
                }
                else {
                    value = "\033[38;2;0;240;0m" + value + "\033[0m";
                }
                s+=value + " ";
            }
            else {
                for (int i = 0; i < potenz; ++i)
                    s+= ".";
                s+= " ";
            }
        }
        s+= "\n";
    }
    cout << s;
    return sBefore + s;
}

string printABCD(const vector<vector<string>>& S, int potenz, vector<double> dualSolution, Biclique biclique, DominoGraph& dominoGraph) {
    string sBefore = printABCD(S, potenz);
    string s;
    int visibleSize = 1 << potenz; // 2^k
    vector<CondensedNode> bot = biclique.getBot();
    vector<CondensedNode> top = biclique.getTop();
    //dual
    //columnlabels
    for (int i = 0; i < potenz; ++i)
        s += " ";
    s += " ";
    for (int i = visibleSize - 1; i >= 0; --i) {
        if (find(bot.begin(), bot.end(), CondensedNode(i, 1, dominoGraph.getPower())) != bot.end()) {
            s += "\033[38;2;0;0;240m" + binaryLabel(i, potenz) + "\033[0m" + " ";
        }
        else s += binaryLabel(i, potenz) + " ";
    }
    s += "\n";

    //rows
    int rowLabel = 0;
    for (int r = 0; r < S.size(); ++r)
    {
        if (!validRowIndex(r))
            continue;
        if (find(top.begin(), top.end(), CondensedNode(rowLabel, 1, dominoGraph.getPower())) != top.end()) {
            s += "\033[38;2;0;0;240m"+binaryLabel(rowLabel, potenz)+ "\033[0m" + " ";
        }
        else s += binaryLabel(rowLabel, potenz) + " ";


        int colLabel = visibleSize;
        for (int c = (int)((S[r].size() / 2) - 1); c < S[r].size(); ++c) {
            if (!validColIndex(c))
                continue;
            colLabel--;
            if (S[r][c] != "n") {
                EdgeOrbit e(S[r][c]);
                double dual = dualSolution[dominoGraph.getEdgeIndex(e)];
                if (find(top.begin(), top.end(), CondensedNode(rowLabel, 1, dominoGraph.getPower())) == top.end()
                    || find(bot.begin(), bot.end(), CondensedNode(colLabel, 1, dominoGraph.getPower())) == bot.end()) {
                    dual = 0;
                }
                string value = to_string(dual);
                if (value[0] == '-')value = value.substr(1, dominoGraph.getPower());
                else value = value.substr(0, dominoGraph.getPower());
                if (dual < 0) {
                    value = "\033[38;2;240;0;0m" + value + "\033[0m";
                }
                else if (dual == 0) {
                    value = "\033[38;2;240;240;240m" + value + "\033[0m";
                }
                else {
                    value = "\033[38;2;0;240;0m" + value + "\033[0m";
                }
                s += value + " ";
            }
            else {
                for (int i = 0; i < potenz; ++i)
                    s += ".";
                s += " ";
            }
        }
        s += "\n";
        rowLabel++;
    }
    cout << s;
    return sBefore + s;
}
