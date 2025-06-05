#ifndef SUBSEQHASH3_HPP
#define SUBSEQHASH3_HPP

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace std;

#define POS_INF numeric_limits<int>::max()
#define NEG_INF numeric_limits<int>::min()

struct BaseDPCell {
    int fracOmega;
    string fracSeed;
    BaseDPCell* prevSubproblem;
};

struct PivotDPCell {
    int omega;
    string seed;
};

struct PiCell {
    int* psi;
    PivotDPCell* seedData;
    int optimalA;
    int optimalB;
};

class SubseqHash3 {
    int k;
    int d;
    map<char, int> alphabet;

    int* tableAF;
    int* tableBF1;
    int* tableBF2;
    int* tableCF;

    int* tableAR;
    int* tableBR1;
    int* tableBR2;
    int* tableCR;

    int* tableAP;
    int* tableBP1;
    int* tableBP2;
    int* tableBP3;
    int* tableCP;

    int returnPivotTableIndex(int, int, int, int);

    void solveForwardDP(string, int, BaseDPCell*, BaseDPCell*);
    void solveReverseDP(string, int, BaseDPCell*, BaseDPCell*);

public:
    SubseqHash3();
    SubseqHash3(int, int, map<char, int>);
    ~SubseqHash3();

    void generateTables();
    void loadTables();

    int getK() const;
    int getD() const;
    map<char, int> getAlphabet() const;

    void solvePivotDP(string, int);
};

#endif