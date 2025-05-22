#include <algorithm>
#include <chrono>
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

    void solveForwardDP(string, int, int*, int*);
    void solveReverseDP(string, int, int*, int*);

public:
    SubseqHash3();
    SubseqHash3(int, int, map<char, int>);
    ~SubseqHash3();

    void solvePivotDP(string, int);
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

    unsigned int currentTimeBasedSeed;

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            currentTimeBasedSeed = chrono::system_clock::now().time_since_epoch().count();
            shuffle(pairValues.begin(), pairValues.end(), default_random_engine(currentTimeBasedSeed));

            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                this->tableBF1[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = (uint8_t) (pairValues[sigma] >> 1) % 2;
                this->tableBF2[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = (uint8_t) (pairValues[sigma] >> 0) % 2;
            }

            currentTimeBasedSeed = chrono::system_clock::now().time_since_epoch().count();
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
                currentTimeBasedSeed = chrono::system_clock::now().time_since_epoch().count();
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
        currentTimeBasedSeed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(dValues.begin(), dValues.end(), default_random_engine(currentTimeBasedSeed));

        for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
            this->tableCF[u * this->alphabet.size() + sigma] = dValues[sigma];
        }

        currentTimeBasedSeed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(dValues.begin(), dValues.end(), default_random_engine(currentTimeBasedSeed));

        for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
            this->tableCR[u * this->alphabet.size() + sigma] = dValues[sigma];
        }
    }

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                currentTimeBasedSeed = chrono::system_clock::now().time_since_epoch().count();
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

void SubseqHash3::solveForwardDP(string sequence, int windowLength, int* dpFmin, int* dpFmax) {
    int N = sequence.length(), n = windowLength;

    // Base case initialization: Forward[w][s][0][0] = 0 else Forward[w][s][u][v] = NaN
    for(int w = 0; w < N; w++) {
        for(int s = 0; s < n + 1; s++) {
            for(int u = 0; u < this->k + 1; u++) {
                for(int v = 0; v < this->d; v++) {
                    dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v] = (u == 0 && v == 0) ? 0 : POS_INF;
                    dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v] = (u == 0 && v == 0) ? 0 : NEG_INF;
                }
            }
        }
    }

    // DP table population: Fmin and Fmax
    for(int w = 0; w < N; w++) {
        for(int s = 1; s <= min(n, N - w); s++) {
            for(int u = 1; u <= min(this->k, N - w); u++) {
                for(int v = 0; v < this->d; v++) {
                    int dpFoption1, dpFoption2, vF = (v - this->tableCF[(u - 1) * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] + this->d) % this->d;
                    
                    // Solving Fmin[w][s][u][v]
                    dpFoption1 = dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v];

                    dpFoption2 = this->tableAF[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] * this->tableBF2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]];
                    dpFoption2 += (this->tableBF1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] == 1) ? dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF] : -dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF];

                    dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v] = min(dpFoption1, dpFoption2);

                    // Solving Fmax[w][s][u][v]
                    dpFoption1 = dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v];

                    dpFoption2 = this->tableAF[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] * this->tableBF2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]];
                    dpFoption2 += (this->tableBF1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] == 1) ? dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF] : -dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF];

                    dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v] = max(dpFoption1, dpFoption2);
                }
            }
        }
    }
}

void SubseqHash3::solveReverseDP(string sequence, int windowLength, int* dpRmin, int* dpRmax) {
    int N = sequence.length(), n = windowLength;

    // Base case initialization: Reverse[w][s][0][0] = 0 else Reverse[w][s][u][v] = NaN
    for(int w = 0; w < N; w++) {
        for(int s = 0; s < n + 1; s++) {
            for(int u = 0; u < this->k + 1; u++) {
                for(int v = 0; v < this->d; v++) {
                    dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v] = (u == 0 && v == 0) ? 0 : POS_INF;
                    dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v] = (u == 0 && v == 0) ? 0 : NEG_INF;
                }
            }
        }
    }

    // DP table population: Rmin and Rmax
    for(int w = N - 1; w >= 0; w--) {
        for(int s = 1; s <= min(n, N - w); s++) {
            for(int u = 1; u <= min(this->k, N - w); u++) {
                for(int v = 0; v < this->d; v++) {
                    int dpRoption1, dpRoption2, vR = (v - this->tableCR[(u - 1) * this->alphabet.size() + this->alphabet[sequence[w]]] + this->d) % this->d;

                    // Solving Rmin[w][s][u][v]
                    if(w == N - 1) {
                        // Special case: Rmin[N][1][1][v] = Rmin[N+1][0][1][v] -> NaN (1 length subseq from 0 length substr) or +Rmin/-Rmax[N+1][0][0][vR] -> Reverse[N+1][0][0][vR]=0 if vR=0 else NaN (0 length subseq with psi=0 from 0 length substr)
                        dpRoption1 = POS_INF;
                        dpRoption2 = (vR == 0) ? this->tableAR[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] * this->tableBR2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] : POS_INF;
                    } else {
                        dpRoption1 = dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v];

                        dpRoption2 = this->tableAR[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] * this->tableBR2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]];
                        dpRoption2 += (this->tableBR1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] == 1) ? dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR] : -dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR];
                    }

                    dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v] = min(dpRoption1, dpRoption2);

                    // Solving Rmax[w][s][u][v]
                    if(w == N - 1) {
                        // Special case: Rmax[N][1][1][v] = Rmax[N+1][0][1][v] -> NaN (1 length subseq from 0 length substr) or +Rmax/-Rmin[N+1][0][0][vR] -> Reverse[N+1][0][0][vR]=0 if vR=0 else NaN (0 length subseq with psi=0 from 0 length substr)
                        dpRoption1 = NEG_INF;
                        dpRoption2 = (vR == 0) ? this->tableAR[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] * this->tableBR2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] : NEG_INF;
                    } else {
                        dpRoption1 = dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v];

                        dpRoption2 = this->tableAR[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] * this->tableBR2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]];
                        dpRoption2 += (this->tableBR1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] == 1) ? dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR] : -dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR];
                    }

                    dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v] = max(dpRoption1, dpRoption2);
                }
            }
        }
    }
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

    delete[] this->tableAR;
    delete[] this->tableBR1;
    delete[] this->tableBR2;
    delete[] this->tableCR;

    delete[] this->tableAP;
    delete[] this->tableBP1;
    delete[] this->tableBP2;
    delete[] this->tableBP3;
    delete[] this->tableCP;
}

void SubseqHash3::solvePivotDP(string sequence, int windowLength) {
    int N = sequence.length(), n = windowLength;

    int *dpFmin, *dpFmax, *dpRmin, *dpRmax;

    // Fmin, Fmax, Rmin, and Rmax have a dimension of [N][n+1][k+1][d]
    dpFmin = new int[N * (n + 1) * (this->k + 1) * this->d];
    dpFmax = new int[N * (n + 1) * (this->k + 1) * this->d];
    dpRmin = new int[N * (n + 1) * (this->k + 1) * this->d];
    dpRmax = new int[N * (n + 1) * (this->k + 1) * this->d];

    solveForwardDP(sequence, windowLength, dpFmin, dpFmax);
    solveReverseDP(sequence, windowLength, dpRmin, dpRmax);

    // Psi and Omega have a dimension of [N - n + 1][n][n][k][k]
    int psi[N - n + 1][n][n][this->k][this->k], omega[N - n + 1][n][n][this->k][this->k];

    // Psi and Omega initialization with d and NaN, respectively
    for(int w = 0; w < N - n + 1; w++) {
        for(int a = 0; a < n; a++) {
            for(int b = 0; b < n; b++) {
                for(int i = 0; i < this->k; i++) {
                    for(int j = 0; j < this->k; j++) {
                        psi[w][a][b][i][j] = this->d;
                        omega[w][a][b][i][j] = NEG_INF;
                    }
                }
            }
        }
    }

    // Psi table population
    for(int w = 0; w < N - n + 1; w++) {
        for(int a = 1; a < n; a++) {
            for(int b = a + 1; b <= n; b++) {
                for(int i = 1; i < this->k; i++) {
                    for(int j = i + 1; j <= this->k; j++) {
                        // Finding minimum psi from all possible combinations of v1, v2, and v3 that satisfy constraints
                        for(int v1 = 0; v1 < this->d; v1++) {
                            for(int v2 = 0; v2 < this->d; v2++) {
                                for(int v3 = 0; v3 < this->d; v3++) {
                                    if(dpRmin[w * (n + 1) * (this->k + 1) * this->d + (a - 1) * (this->k + 1) * this->d + (i - 1) * this->d + v1] < POS_INF && dpFmin[(w + a) * (n + 1) * (this->k + 1) * this->d + (b - a - 1) * (this->k + 1) * this->d + (j - i - 1) * this->d + v2] < POS_INF && dpRmin[(w + b) * (n + 1) * (this->k + 1) * this->d + (n - b) * (this->k + 1) * this->d + (k - j) * this->d + v3] < POS_INF) {
                                        psi[w][a - 1][b - 1][i - 1][j - 1] = min(psi[w][a - 1][b - 1][i - 1][j - 1], (this->tableCP[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] + v1 + v2 + v3) % this->d);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Omega table popultation
    for(int w = 0; w < N - n + 1; w++) {
        for(int a = 1; a < n; a++) {
            for(int b = a + 1; b <= n; b++) {
                for(int i = 1; i < this->k; i++) {
                    for(int j = i + 1; j <= this->k; j++) {
                        // Finding maximum omega from all possible combinations of v1 and v2
                        for(int v1 = 0; v1 < this->d; v1++) {
                            for(int v2 = 0; v2 < this->d; v2++) {
                                int v3, reverse1, forward2, reverse3;

                                v3 = (psi[w][a - 1][b - 1][i - 1][j - 1] - this->tableCP[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] - v1 - v2 + this->d) % this->d;

                                reverse1 = (this->tableBP1[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] == 1) ? dpRmax[w * (n + 1) * (this->k + 1) * this->d + (a - 1) * (this->k + 1) * this->d + (i - 1) * this->d + v1] : -dpRmin[w * (n + 1) * (this->k + 1) * this->d + (a - 1) * (this->k + 1) * this->d + (i - 1) * this->d + v1];
                                forward2 = (this->tableBP2[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] == 1) ? dpFmax[(w + a) * (n + 1) * (this->k + 1) * this->d + (b - a - 1) * (this->k + 1) * this->d + (j - i - 1) * this->d + v2] : -dpFmin[(w + a) * (n + 1) * (this->k + 1) * this->d + (b - a - 1) * (this->k + 1) * this->d + (j - i - 1) * this->d + v2];
                                reverse3 = (this->tableBP3[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] == 1) ? dpRmax[(w + b) * (n + 1) * (this->k + 1) * this->d + (n - b) * (this->k + 1) * this->d + (k - j) * this->d + v3] : -dpRmin[(w + b) * (n + 1) * (this->k + 1) * this->d + (n - b) * (this->k + 1) * this->d + (k - j) * this->d + v3];

                                omega[w][a - 1][b - 1][i - 1][j - 1] = max(omega[w][a - 1][b - 1][i - 1][j - 1], this->tableAP[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] + reverse1 + forward2 + reverse3);
                            }
                        }
                    }
                }
            }
        }
    }

    delete[] dpFmin;
    delete[] dpFmax;
    delete[] dpRmin;
    delete[] dpRmax;
}