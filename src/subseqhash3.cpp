#include<algorithm>
#include<chrono>
#include<fstream>
#include<iostream>
#include<limits>
#include<map>
#include<random>
#include<string>
#include<vector>

using namespace std;

#define POS_INF numeric_limits<int>::max()
#define NEG_INF numeric_limits<int>::min()

class SubseqHash3 {
    int k;
    int d;
    map<char, int> alphabet;

    int* tableAF;
    uint8_t* tableBF1;
    uint8_t* tableBF2;
    int* tableCF;

    int* tableAR;
    uint8_t* tableBR1;
    uint8_t* tableBR2;
    int* tableCR;

    int* tableAP;
    uint8_t* tableBP1;
    uint8_t* tableBP2;
    uint8_t* tableBP3;
    int* tableCP;

    int returnPivotTableIndex(int, int, int, int);
    void generateTables();
    void loadTables();

public:
    SubseqHash3();
    SubseqHash3(int, int, map<char, int>);
    ~SubseqHash3();
};

int SubseqHash3::returnPivotTableIndex(int i, int j, int sigmaI, int sigmaJ) {
    int pivotTableIndex = 0;

    for(int iIndex = 0; iIndex < i; iIndex++) {
        pivotTableIndex += (this->k - (iIndex + 1)) * this->alphabet.size() * this->alphabet.size();
    }

    pivotTableIndex += (j - (i + 1)) * this->alphabet.size() * this->alphabet.size();
    pivotTableIndex += sigmaI * this->alphabet.size() + sigmaJ;

    return pivotTableIndex;
}

void SubseqHash3::generateTables() {
    random_device seedSource;
    mt19937 randomNumberEngine(seedSource());
    uniform_int_distribution<int> uniformIntegerDistribution(1 << 5, 1 << 10);

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                this->tableAF[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = uniformIntegerDistribution(randomNumberEngine);
                this->tableAR[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = uniformIntegerDistribution(randomNumberEngine);
            }
        }
    }

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    this->tableAP[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)] = uniformIntegerDistribution(randomNumberEngine);
                }
            }
        }
    }

    vector<uint8_t> pairValues, tripletValues;

    for(uint8_t i = 0; i < 4; i++) {
        pairValues.push_back(i);
    }

    if(this->alphabet.size() > 4) {
        for(uint8_t i = 4; i < this->alphabet.size(); i++) {
            pairValues.push_back(i % 4);
        }
    }

    for(uint8_t i = 0; i < 8; i++) {
        tripletValues.push_back(i);
    }

    if(this->alphabet.size() > 8) {
        for(uint8_t i = 8; i < this->alphabet.size(); i++) {
            tripletValues.push_back(i % 8);
        }
    }

    unsigned currentTimeBasedSeed;

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            currentTimeBasedSeed = unsigned(chrono::system_clock::now().time_since_epoch().count());
            shuffle(pairValues.begin(), pairValues.end(), default_random_engine(currentTimeBasedSeed));

            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                this->tableBF1[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = (uint8_t) (pairValues[sigma] >> 1) % 2;
                this->tableBF2[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = (uint8_t) (pairValues[sigma] >> 0) % 2;
            }

            currentTimeBasedSeed = unsigned(chrono::system_clock::now().time_since_epoch().count());
            shuffle(pairValues.begin(), pairValues.end(), default_random_engine(currentTimeBasedSeed));

            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                this->tableBR1[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = (uint8_t) (pairValues[sigma] >> 1) % 2;
                this->tableBR2[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = (uint8_t) (pairValues[sigma] >> 0) % 2;
            }
        }
    }

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                currentTimeBasedSeed = unsigned(chrono::system_clock::now().time_since_epoch().count());
                shuffle(tripletValues.begin(), tripletValues.end(), default_random_engine(currentTimeBasedSeed));

                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    this->tableBP1[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)] = (tripletValues[sigmaJ] >> 2) % 2;
                    this->tableBP2[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)] = (tripletValues[sigmaJ] >> 1) % 2;
                    this->tableBP3[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)] = (tripletValues[sigmaJ] >> 0) % 2;
                }
            }
        }
    }

    vector<int> dValues;

    for(int i = 0; i < this->d; i++) {
        dValues.push_back(i);
    }

    if(this->alphabet.size() > d) {
        for(int i = d; i < this->alphabet.size(); i++) {
            dValues.push_back(i % d);
        }
    }

    for(int u = 0; u < this->k; u++) {
        currentTimeBasedSeed = unsigned(chrono::system_clock::now().time_since_epoch().count());
        shuffle(dValues.begin(), dValues.end(), default_random_engine(currentTimeBasedSeed));

        for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
            this->tableCF[u * this->alphabet.size() + sigma] = dValues[sigma];
        }

        currentTimeBasedSeed = unsigned(chrono::system_clock::now().time_since_epoch().count());
        shuffle(dValues.begin(), dValues.end(), default_random_engine(currentTimeBasedSeed));

        for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
            this->tableCR[u * this->alphabet.size() + sigma] = dValues[sigma];
        }
    }

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                currentTimeBasedSeed = unsigned(chrono::system_clock::now().time_since_epoch().count());
                shuffle(dValues.begin(), dValues.end(), default_random_engine(currentTimeBasedSeed));

                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    this->tableCP[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)] = dValues[sigmaJ];
                }
            }
        }
    }
}

void SubseqHash3::loadTables() {

}

SubseqHash3::SubseqHash3() : SubseqHash3(24, 11, {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}}) {
    // Parameterized constructor is called from default constructor with default arguments
}

SubseqHash3::SubseqHash3(int k, int d, map<char, int> alphabet) {
    this->k = k;
    this->d = d;
    this->alphabet = alphabet;

    this->tableAF = new int[k * d * alphabet.size()];
    this->tableBF1 = new uint8_t[k * d * alphabet.size()];
    this->tableBF2 = new uint8_t[k * d * alphabet.size()];
    this->tableCF = new int[k * alphabet.size()];

    this->tableAR = new int[k * d * alphabet.size()];
    this->tableBR1 = new uint8_t[k * d * alphabet.size()];
    this->tableBR2 = new uint8_t[k * d * alphabet.size()];
    this->tableCR = new int[k * alphabet.size()];

    this->tableAP = new int[(k * (k - 1) * alphabet.size() * alphabet.size()) / 2];
    this->tableBP1 = new uint8_t[(k * (k - 1) * alphabet.size() * alphabet.size()) / 2];
    this->tableBP2 = new uint8_t[(k * (k - 1) * alphabet.size() * alphabet.size()) / 2];
    this->tableBP3 = new uint8_t[(k * (k - 1) * alphabet.size() * alphabet.size()) / 2];
    this->tableCP = new int[(k * (k - 1) * alphabet.size() * alphabet.size()) / 2];

    this->generateTables();
}

SubseqHash3::~SubseqHash3() {
    delete[] this->tableAF;
    delete[] this->tableBF1;
    delete[] this->tableBF2;
    delete[] this->tableCF;

    this->tableAF = nullptr;
    this->tableBF1 = nullptr;
    this->tableBF2 = nullptr;
    this->tableCF = nullptr;

    delete[] this->tableAR;
    delete[] this->tableBR1;
    delete[] this->tableBR2;
    delete[] this->tableCR;

    this->tableAR = nullptr;
    this->tableBR1 = nullptr;
    this->tableBR2 = nullptr;
    this->tableCR = nullptr;

    delete[] this->tableAP;
    delete[] this->tableBP1;
    delete[] this->tableBP2;
    delete[] this->tableBP3;
    delete[] this->tableCP;

    this->tableAP = nullptr;
    this->tableBP1 = nullptr;
    this->tableBP2 = nullptr;
    this->tableBP3 = nullptr;
    this->tableCP = nullptr;
}