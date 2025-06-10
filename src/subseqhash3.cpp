#include "subseqhash3.hpp"

int SubseqHash3::returnPivotTableIndex(int i, int j, int sigmaI, int sigmaJ) {
    int pivotTableIndex = 0;

    for(int iIndex = 0; iIndex < i; iIndex++) {
        pivotTableIndex += (this->k - (iIndex + 1)) * this->alphabet.size() * this->alphabet.size();
    }

    pivotTableIndex += (j - (i + 1)) * this->alphabet.size() * this->alphabet.size();
    pivotTableIndex += sigmaI * this->alphabet.size() + sigmaJ;

    return pivotTableIndex;
}

void SubseqHash3::solveForwardDP(string sequence, int windowLength, BaseDPCell* dpFmin, BaseDPCell* dpFmax) {
    int N = sequence.length(), n = windowLength;

    // Base case initialization: Forward[w][s][0][0] = 0 (as long as [s] within window starting from [w]) else Forward[w][s][u][v] = NaN
    for(int w = 0; w < N; w++) {
        for(int s = 0; s < n + 1; s++) {
            for(int u = 0; u < this->k + 1; u++) {
                for(int v = 0; v < this->d; v++) {
                    if(s <= N - w && u == 0 && v == 0) {
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = 0;
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("");
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;

                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = 0;
                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("");
                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;
                    } else {
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = POS_INF;
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("X");
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;

                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = NEG_INF;
                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("X");
                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;
                    }
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
                    dpFoption1 = dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v].fracOmega;

                    if(this->tableBF1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] == 1) {
                        dpFoption2 = (dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracOmega < POS_INF) ? (this->tableAF[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] * this->tableBF2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] + dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracOmega) : POS_INF;
                    } else {
                        dpFoption2 = (dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracOmega > NEG_INF) ? (this->tableAF[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] * this->tableBF2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] + -dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracOmega) : POS_INF;
                    }

                    if(dpFoption1 < dpFoption2) {
                        // Do not add X[w + s - 1] to seed
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpFoption1;
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v].fracSeed + string("");
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = &dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v];
                    } else {
                        // Add X[w + s - 1] to seed
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpFoption2;
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = ((this->tableBF1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] == 1) ? dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracSeed : dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracSeed) + sequence[w + s - 1];
                        dpFmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = (this->tableBF1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] == 1) ? &dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF] : &dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF];
                    }

                    // Solving Fmax[w][s][u][v]
                    dpFoption1 = dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v].fracOmega;

                    if(this->tableBF1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] == 1) {
                        dpFoption2 = (dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracOmega > NEG_INF) ? (this->tableAF[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] * this->tableBF2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] + dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracOmega) : NEG_INF;
                    } else {
                        dpFoption2 = (dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracOmega < POS_INF) ? (this->tableAF[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] * this->tableBF2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] + -dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracOmega) : NEG_INF;
                    }

                    if(dpFoption1 > dpFoption2) {
                        // Do not add X[w + s - 1] to seed
                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpFoption1;
                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v].fracSeed + string("");
                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = &dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v];
                    } else {
                        // Add X[w + s - 1] to seed
                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpFoption2;
                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = ((this->tableBF1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] == 1) ? dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracSeed : dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF].fracSeed) + sequence[w + s - 1];
                        dpFmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = (this->tableBF1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w + s - 1]]] == 1) ? &dpFmax[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF] : &dpFmin[w * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vF];
                    }
                }
            }
        }
    }
}

void SubseqHash3::solveReverseDP(string sequence, int windowLength, BaseDPCell* dpRmin, BaseDPCell* dpRmax) {
    int N = sequence.length(), n = windowLength;

    // Base case initialization: Reverse[w][s][0][0] = 0 (as long as [s] within window starting from [w]) else Reverse[w][s][u][v] = NaN
    for(int w = 0; w < N; w++) {
        for(int s = 0; s < n + 1; s++) {
            for(int u = 0; u < this->k + 1; u++) {
                for(int v = 0; v < this->d; v++) {
                    if(s <= N - w && u == 0 && v == 0) {
                        dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = 0;
                        dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("");
                        dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;

                        dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = 0;
                        dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("");
                        dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;
                    } else {
                        dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = POS_INF;
                        dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("X");
                        dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;

                        dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = NEG_INF;
                        dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("X");
                        dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;
                    }
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

                        if(dpRoption1 < dpRoption2) {
                            // Do not add X[w] to seed
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpRoption1;
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("X");
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;
                        } else {
                            // Add X[w] to seed but condition applies (vR must be equal to 0, in other words, 0-length subseq from previous subproblem must have psi-value of 0)
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpRoption2;
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = (vR == 0) ? sequence[w] + string("") : sequence[w] + string("X");
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;
                        }
                    } else {
                        dpRoption1 = dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v].fracOmega;

                        if(this->tableBR1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] == 1) {
                            dpRoption2 = (dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracOmega < POS_INF) ? (this->tableAR[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] * this->tableBR2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] + dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracOmega) : POS_INF;
                        } else {
                            dpRoption2 = (dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracOmega > NEG_INF) ? (this->tableAR[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] * this->tableBR2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] + -dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracOmega) : POS_INF;
                        }

                        if(dpRoption1 < dpRoption2) {
                            // Do not add X[w] to seed
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpRoption1;
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("") + dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v].fracSeed;
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = &dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v];
                        } else {
                            // Add X[w] to seed
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpRoption2;
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = sequence[w] + ((this->tableBR1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] == 1) ? dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracSeed : dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracSeed);
                            dpRmin[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = (this->tableBR1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] == 1) ? &dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR] : &dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR];
                        }
                    }

                    // Solving Rmax[w][s][u][v]
                    if(w == N - 1) {
                        // Special case: Rmax[N][1][1][v] = Rmax[N+1][0][1][v] -> NaN (1 length subseq from 0 length substr) or +Rmax/-Rmin[N+1][0][0][vR] -> Reverse[N+1][0][0][vR]=0 if vR=0 else NaN (0 length subseq with psi=0 from 0 length substr)
                        dpRoption1 = NEG_INF;
                        dpRoption2 = (vR == 0) ? this->tableAR[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] * this->tableBR2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] : NEG_INF;

                        if(dpRoption1 > dpRoption2) {
                            // Do not add X[w] to seed
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpRoption1;
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("X");
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;
                        } else {
                            // Add X[w] to seed but condition applies (vR must be equal to 0, in other words, 0-length subseq from previous subproblem must have psi-value of 0)
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpRoption2;
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = (vR == 0) ? sequence[w] + string("") : sequence[w] + string("X");
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = nullptr;
                        }
                    } else {
                        dpRoption1 = dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v].fracOmega;

                        if(this->tableBR1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] == 1) {
                            dpRoption2 = (dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracOmega > NEG_INF) ? (this->tableAR[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] * this->tableBR2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] + dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracOmega) : NEG_INF;
                        } else {
                            dpRoption2 = (dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracOmega < POS_INF) ? (this->tableAR[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] * this->tableBR2[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] + -dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracOmega) : NEG_INF;
                        }

                        if(dpRoption1 > dpRoption2) {
                            // Do not add X[w] to seed
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpRoption1;
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = string("") + dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v].fracSeed;
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = &dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + u * this->d + v];
                        } else {
                            // Add X[w] to seed
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracOmega = dpRoption2;
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].fracSeed = sequence[w] + ((this->tableBR1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] == 1) ? dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracSeed : dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR].fracSeed);
                            dpRmax[w * (n + 1) * (this->k + 1) * this->d + s * (this->k + 1) * this->d + u * this->d + v].prevSubproblem = (this->tableBR1[(u - 1) * this->d * this->alphabet.size() + v * this->alphabet.size() + this->alphabet[sequence[w]]] == 1) ? &dpRmax[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR] : &dpRmin[(w + 1) * (n + 1) * (this->k + 1) * this->d + (s - 1) * (this->k + 1) * this->d + (u - 1) * this->d + vR];
                        }
                    }
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
    this->tableBF1 = new int[k * d * alphabet.size()];
    this->tableBF2 = new int[k * d * alphabet.size()];
    this->tableCF = new int[k * alphabet.size()];

    this->tableAR = new int[k * d * alphabet.size()];
    this->tableBR1 = new int[k * d * alphabet.size()];
    this->tableBR2 = new int[k * d * alphabet.size()];
    this->tableCR = new int[k * alphabet.size()];

    this->tableAP = new int[(k * (k - 1) * alphabet.size() * alphabet.size()) / 2];
    this->tableBP1 = new int[(k * (k - 1) * alphabet.size() * alphabet.size()) / 2];
    this->tableBP2 = new int[(k * (k - 1) * alphabet.size() * alphabet.size()) / 2];
    this->tableBP3 = new int[(k * (k - 1) * alphabet.size() * alphabet.size()) / 2];
    this->tableCP = new int[(k * (k - 1) * alphabet.size() * alphabet.size()) / 2];
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

    vector<int> pairValues, tripletValues;

    for(int i = 0; i < 4; i++) {
        pairValues.push_back(i);
    }

    if(this->alphabet.size() > 4) {
        for(int i = 4; i < this->alphabet.size(); i++) {
            pairValues.push_back(i % 4);
        }
    }

    for(int i = 0; i < 8; i++) {
        tripletValues.push_back(i);
    }

    if(this->alphabet.size() > 8) {
        for(int i = 8; i < this->alphabet.size(); i++) {
            tripletValues.push_back(i % 8);
        }
    }

    unsigned int currentTimeBasedSeed;

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            currentTimeBasedSeed = chrono::system_clock::now().time_since_epoch().count();
            shuffle(pairValues.begin(), pairValues.end(), default_random_engine(currentTimeBasedSeed));

            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                this->tableBF1[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = (pairValues[sigma] >> 1) % 2;
                this->tableBF2[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = (pairValues[sigma] >> 0) % 2;
            }

            currentTimeBasedSeed = chrono::system_clock::now().time_since_epoch().count();
            shuffle(pairValues.begin(), pairValues.end(), default_random_engine(currentTimeBasedSeed));

            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                this->tableBR1[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = (pairValues[sigma] >> 1) % 2;
                this->tableBR2[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] = (pairValues[sigma] >> 0) % 2;
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

    if(this->alphabet.size() > this->d) {
        for(int i = this->d; i < this->alphabet.size(); i++) {
            dValues.push_back(i % this->d);
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

    string tablesDirectoryPath = string("..") + filesystem::path::preferred_separator + string("tables");
    filesystem::path tablesDirectory(tablesDirectoryPath);

    if(!filesystem::exists(tablesDirectory)) {
        if(!filesystem::create_directories(tablesDirectory)) {
            exit(EXIT_FAILURE);
        }
    }

    string tablesFilePath = string("..") + filesystem::path::preferred_separator + string("tables") + filesystem::path::preferred_separator + string("tables_k") + to_string(this->k) + string("_d") + to_string(this->d) + string("_sigma") + to_string(this->alphabet.size());
    ofstream tablesFile(tablesFilePath);

    if(!tablesFile.is_open()) {
        exit(EXIT_FAILURE);
    }

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile << this->tableAF[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] << " ";
            }

            tablesFile << endl;
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile << this->tableAR[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] << " ";
            }

            tablesFile << endl;
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    tablesFile << this->tableAP[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)] << " ";
                }

                tablesFile << endl;
            }

            tablesFile << endl;
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile << this->tableBF1[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] << " ";
            }

            tablesFile << endl;
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile << this->tableBF2[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] << " ";
            }

            tablesFile << endl;
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile << this->tableBR1[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] << " ";
            }

            tablesFile << endl;
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile << this->tableBR2[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma] << " ";
            }

            tablesFile << endl;
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    tablesFile << this->tableBP1[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)] << " ";
                }

                tablesFile << endl;
            }

            tablesFile << endl;
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    tablesFile << this->tableBP2[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)] << " ";
                }

                tablesFile << endl;
            }

            tablesFile << endl;
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    tablesFile << this->tableBP3[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)] << " ";
                }

                tablesFile << endl;
            }

            tablesFile << endl;
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int u = 0; u < this->k; u++) {
        for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
            tablesFile << this->tableCF[u * this->alphabet.size() + sigma] << " ";
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int u = 0; u < this->k; u++) {
        for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
            tablesFile << this->tableCR[u * this->alphabet.size() + sigma] << " ";
        }

        tablesFile << endl;
    }

    tablesFile << endl;

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    tablesFile << this->tableCP[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)] << " ";
                }

                tablesFile << endl;
            }

            tablesFile << endl;
        }

        tablesFile << endl;
    }

    tablesFile << endl;
    tablesFile.close();
}

void SubseqHash3::loadTables() {
    string tablesFilePath = string("..") + filesystem::path::preferred_separator + string("tables") + filesystem::path::preferred_separator + string("tables_k") + to_string(this->k) + string("_d") + to_string(this->d) + string("_sigma") + to_string(this->alphabet.size());
    ifstream tablesFile;

    tablesFile.open(tablesFilePath);

    if(!tablesFile.is_open()) {
        exit(EXIT_FAILURE);
    }

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile >> this->tableAF[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma];
            }
        }
    }

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile >> this->tableAR[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma];
            }
        }
    }

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    tablesFile >> this->tableAP[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)];
                }
            }
        }
    }

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile >> this->tableBF1[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma];
            }
        }
    }

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile >> this->tableBF2[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma];
            }
        }
    }

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile >> this->tableBR1[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma];
            }
        }
    }

    for(int u = 0; u < this->k; u++) {
        for(int v = 0; v < this->d; v++) {
            for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
                tablesFile >> this->tableBR2[u * this->d * this->alphabet.size() + v * this->alphabet.size() + sigma];
            }
        }
    }

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    tablesFile >> this->tableBP1[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)];
                }
            }
        }
    }

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    tablesFile >> this->tableBP2[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)];
                }
            }
        }
    }

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    tablesFile >> this->tableBP3[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)];
                }
            }
        }
    }

    for(int u = 0; u < this->k; u++) {
        for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
            tablesFile >> this->tableCF[u * this->alphabet.size() + sigma];
        }
    }

    for(int u = 0; u < this->k; u++) {
        for(int sigma = 0; sigma < this->alphabet.size(); sigma++) {
            tablesFile >> this->tableCR[u * this->alphabet.size() + sigma];
        }
    }

    for(int i = 0; i < this->k - 1; i++) {
        for(int j = i + 1; j < this->k; j++) {
            for(int sigmaI = 0; sigmaI < this->alphabet.size(); sigmaI++) {
                for(int sigmaJ = 0; sigmaJ < this->alphabet.size(); sigmaJ++) {
                    tablesFile >> this->tableCP[this->returnPivotTableIndex(i, j, sigmaI, sigmaJ)];
                }
            }
        }
    }

    tablesFile.close();
}

int SubseqHash3::getK() const {
    return this->k;
}

int SubseqHash3::getD() const {
    return this->d;
}

map<char, int> SubseqHash3::getAlphabet() const {
    return this->alphabet;
}

vector<PiCell> SubseqHash3::solvePivotDP(string sequence, int windowLength) {
    int N = sequence.length(), n = windowLength;

    BaseDPCell *dpFmin, *dpFmax, *dpRmin, *dpRmax;

    // Fmin, Fmax, Rmin, and Rmax have a dimension of [N][n+1][k+1][d]
    dpFmin = new BaseDPCell[N * (n + 1) * (this->k + 1) * this->d];
    dpFmax = new BaseDPCell[N * (n + 1) * (this->k + 1) * this->d];
    dpRmin = new BaseDPCell[N * (n + 1) * (this->k + 1) * this->d];
    dpRmax = new BaseDPCell[N * (n + 1) * (this->k + 1) * this->d];

    solveForwardDP(sequence, windowLength, dpFmin, dpFmax);
    solveReverseDP(sequence, windowLength, dpRmin, dpRmax);

    // Psi and Omega have a dimension of [N - n + 1][n][n][k][k]
    // [future task] We may want to switch to dynamically allocated arrays with non-trivial access to avoid storing unused array cells
    int psi[N - n + 1][n][n][this->k][this->k];
    PivotDPCell omega[N - n + 1][n][n][this->k][this->k];
    
    // Psi and Omega initialization with d and {NaN, empty_string}, respectively
    for(int w = 0; w < N - n + 1; w++) {
        for(int a = 0; a < n; a++) {
            for(int b = 0; b < n; b++) {
                for(int i = 0; i < this->k; i++) {
                    for(int j = 0; j < this->k; j++) {
                        psi[w][a][b][i][j] = this->d;

                        omega[w][a][b][i][j].omega = NEG_INF;
                        omega[w][a][b][i][j].seed = string("");
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
                                // We need to explicitly handle the case when w + b = N, that is, there is an empty substr after second pivot position, or in other words, last character is second pivot
                                if(w + b < N) {
                                    for(int v3 = 0; v3 < this->d; v3++) {
                                        if(dpRmin[w * (n + 1) * (this->k + 1) * this->d + (a - 1) * (this->k + 1) * this->d + (i - 1) * this->d + v1].fracOmega < POS_INF && dpFmin[(w + a) * (n + 1) * (this->k + 1) * this->d + (b - a - 1) * (this->k + 1) * this->d + (j - i - 1) * this->d + v2].fracOmega < POS_INF && dpRmin[(w + b) * (n + 1) * (this->k + 1) * this->d + (n - b) * (this->k + 1) * this->d + (this->k - j) * this->d + v3].fracOmega < POS_INF) {
                                            psi[w][a - 1][b - 1][i - 1][j - 1] = min(psi[w][a - 1][b - 1][i - 1][j - 1], (this->tableCP[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] + v1 + v2 + v3) % this->d);
                                        }
                                    }
                                } else {
                                    if(dpRmin[w * (n + 1) * (this->k + 1) * this->d + (a - 1) * (this->k + 1) * this->d + (i - 1) * this->d + v1].fracOmega < POS_INF && dpFmin[(w + a) * (n + 1) * (this->k + 1) * this->d + (b - a - 1) * (this->k + 1) * this->d + (j - i - 1) * this->d + v2].fracOmega < POS_INF && j == this->k) {
                                        // v3 = 0 only when subseq length is also 0, or in other words, only a psi value of 0 is acceptable for a 0-length subseq
                                        psi[w][a - 1][b - 1][i - 1][j - 1] = min(psi[w][a - 1][b - 1][i - 1][j - 1], (this->tableCP[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] + v1 + v2) % this->d);
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
                                int v3, reverseFracOmega1, forwardFracOmega2, reverseFracOmega3;
                                string reverseFracSeed1, forwardFracSeed2, reverseFracSeed3;

                                v3 = (psi[w][a - 1][b - 1][i - 1][j - 1] - this->tableCP[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] - v1 - v2 + 3 * this->d) % this->d;

                                if(this->tableBP1[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] == 1) {
                                    reverseFracOmega1 = dpRmax[w * (n + 1) * (this->k + 1) * this->d + (a - 1) * (this->k + 1) * this->d + (i - 1) * this->d + v1].fracOmega;
                                    reverseFracSeed1 = dpRmax[w * (n + 1) * (this->k + 1) * this->d + (a - 1) * (this->k + 1) * this->d + (i - 1) * this->d + v1].fracSeed;
                                } else {
                                    reverseFracOmega1 = -dpRmin[w * (n + 1) * (this->k + 1) * this->d + (a - 1) * (this->k + 1) * this->d + (i - 1) * this->d + v1].fracOmega;
                                    reverseFracSeed1 = dpRmin[w * (n + 1) * (this->k + 1) * this->d + (a - 1) * (this->k + 1) * this->d + (i - 1) * this->d + v1].fracSeed;
                                }

                                if(this->tableBP2[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] == 1) {
                                    forwardFracOmega2 = dpFmax[(w + a) * (n + 1) * (this->k + 1) * this->d + (b - a - 1) * (this->k + 1) * this->d + (j - i - 1) * this->d + v2].fracOmega;
                                    forwardFracSeed2 = dpFmax[(w + a) * (n + 1) * (this->k + 1) * this->d + (b - a - 1) * (this->k + 1) * this->d + (j - i - 1) * this->d + v2].fracSeed;
                                } else {
                                    forwardFracOmega2 = -dpFmin[(w + a) * (n + 1) * (this->k + 1) * this->d + (b - a - 1) * (this->k + 1) * this->d + (j - i - 1) * this->d + v2].fracOmega;
                                    forwardFracSeed2 = dpFmin[(w + a) * (n + 1) * (this->k + 1) * this->d + (b - a - 1) * (this->k + 1) * this->d + (j - i - 1) * this->d + v2].fracSeed;
                                }

                                // We need to explicitly handle the case when w + b = N, that is, there is an empty substr after second pivot position, or in other words, last character is second pivot
                                if(this->tableBP3[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] == 1) {
                                    reverseFracOmega3 = (w + b < N) ? dpRmax[(w + b) * (n + 1) * (this->k + 1) * this->d + (n - b) * (this->k + 1) * this->d + (this->k - j) * this->d + v3].fracOmega : ((j == this->k) ? 0 : NEG_INF);
                                    reverseFracSeed3 = (w + b < N) ? dpRmax[(w + b) * (n + 1) * (this->k + 1) * this->d + (n - b) * (this->k + 1) * this->d + (this->k - j) * this->d + v3].fracSeed : ((j == this->k) ? string("") : string("X"));
                                } else {
                                    reverseFracOmega3 = (w + b < N) ? -dpRmin[(w + b) * (n + 1) * (this->k + 1) * this->d + (n - b) * (this->k + 1) * this->d + (this->k - j) * this->d + v3].fracOmega : ((j == this->k) ? 0 : NEG_INF);
                                    reverseFracSeed3 = (w + b < N) ? dpRmin[(w + b) * (n + 1) * (this->k + 1) * this->d + (n - b) * (this->k + 1) * this->d + (this->k - j) * this->d + v3].fracSeed : ((j == this->k) ? string("") : string("X"));
                                }

                                if(reverseFracOmega1 > -POS_INF && forwardFracOmega2 > -POS_INF && reverseFracOmega3 > -POS_INF) {
                                    if(omega[w][a - 1][b - 1][i - 1][j - 1].omega < this->tableAP[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] + reverseFracOmega1 + forwardFracOmega2 + reverseFracOmega3) {
                                        omega[w][a - 1][b - 1][i - 1][j - 1].omega = this->tableAP[this->returnPivotTableIndex(i - 1, j - 1, this->alphabet[sequence[w + a - 1]], this->alphabet[sequence[w + b - 1]])] + reverseFracOmega1 + forwardFracOmega2 + reverseFracOmega3;
                                        omega[w][a - 1][b - 1][i - 1][j - 1].seed = reverseFracSeed1 + sequence[w + a - 1] + forwardFracSeed2 + sequence[w + b - 1] + reverseFracSeed3;
                                    }
                                }
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

    // Pi has a dimension of [N - n + 1][k][k]
    // [future task] We may want to switch to a dynamically allocated array with non-trivial access to avoid storing unused array cells
    PiCell pi[N - n + 1][this->k][this->k];

    for(int w = 0; w < N - n + 1; w++) {
        for(int i = 0; i < this->k - 1; i++) {
            for(int j = i + 1; j < this->k; j++) {
                // Finding optimal seed from all possible combinations of pivot positions a and b
                for(int a = 0; a < n - 1; a++) {
                    for(int b = a + 1; b < n; b++) {
                        if(a == 0 && b == 1) {
                            pi[w][i][j].psi = &psi[w][a][b][i][j];
                            pi[w][i][j].seedData = &omega[w][a][b][i][j];
                            pi[w][i][j].windowStartPosition = w + 1;
                            pi[w][i][j].optimalA = a + 1;
                            pi[w][i][j].optimalB = b + 1;
                            pi[w][i][j].pivotI = i + 1;
                            pi[w][i][j].pivotJ = j + 1;
                        } else {
                            if(*pi[w][i][j].psi > psi[w][a][b][i][j] || (*pi[w][i][j].psi == psi[w][a][b][i][j] && pi[w][i][j].seedData->omega < omega[w][a][b][i][j].omega)) {
                                pi[w][i][j].psi = &psi[w][a][b][i][j];
                                pi[w][i][j].seedData = &omega[w][a][b][i][j];
                                pi[w][i][j].windowStartPosition = w + 1;
                                pi[w][i][j].optimalA = a + 1;
                                pi[w][i][j].optimalB = b + 1;
                                pi[w][i][j].pivotI = i + 1;
                                pi[w][i][j].pivotJ = j + 1;
                            }
                        }
                    }
                }
            }
        }
    }

    vector<PiCell> seeds;

    for(int w = 0; w < N - n + 1; w++) {
        for(int i = 0; i < this->k - 1; i++) {
            for(int j = i + 1; j < this->k; j++) {
                seeds.push_back(pi[w][i][j]);
            }
        }
    }
    
    return seeds;
}