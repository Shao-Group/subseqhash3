#include "subseqhash3.hpp"

bool comparePiCellByOmega(const PiCell& piCell1, const PiCell& piCell2) {
    return piCell1.omega > piCell2.omega;
}

int main(int argc, char** argv) {
    int k = stoi(argv[1]), d = stoi(argv[2]);

    SubseqHash3 subseqHash3(k, d, {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}});

    // subseqHash3.generateTables();
    subseqHash3.loadTables();

    string sequencePairsFilePath = argv[3], sequence1, sequence2;

    ifstream sequencePairsFile;

    sequencePairsFile.open(sequencePairsFilePath);

    if(!sequencePairsFile.is_open()) {
        exit(EXIT_FAILURE);
    }

    bool bDoSingleComparison = true, bVerbose = false, bDoHigherOmegaComparison = true;
    int numHashCollisions = 0, numTotalSequencePairs = 0, numTopPicks = 1;

    while(sequencePairsFile >> sequence1) {
        sequencePairsFile >> sequence2;

        vector<PiCell> sequence1Seeds = subseqHash3.solvePivotDP(sequence1, sequence1.length()), sequence2Seeds = subseqHash3.solvePivotDP(sequence2, sequence2.length());

        if(sequence1Seeds.size() != sequence2Seeds.size()) {
            exit(EXIT_FAILURE);
        }

        if(bVerbose) {
            cout << "Sequence1: " << sequence1 << "\nSequence2: " << sequence2 << "\n\nwindow_start (w), first_pivot (i), second_pivot (j), first_pivot_win_pos (a), second_pivot_win_pos (b), psi, omega, seed\n" << endl;

            for(int i = 0; i < sequence1Seeds.size(); i++) {
                cout << sequence1Seeds[i].windowStartPosition << ", " << sequence1Seeds[i].pivotI << ", " << sequence1Seeds[i].pivotJ << ", " << sequence1Seeds[i].optimalA << ", " << sequence1Seeds[i].optimalB << ", " << sequence1Seeds[i].psi << ", " << sequence1Seeds[i].omega << ", " << sequence1Seeds[i].seed << "\n" << sequence2Seeds[i].windowStartPosition << ", " << sequence2Seeds[i].pivotI << ", " << sequence2Seeds[i].pivotJ << ", " << sequence2Seeds[i].optimalA << ", " << sequence2Seeds[i].optimalB << ", " << sequence2Seeds[i].psi << ", " << sequence2Seeds[i].omega << ", " << sequence2Seeds[i].seed << "\n" << endl;
            }

            cout << "-\n" << endl;
        }

        if(bDoHigherOmegaComparison) {
            // Omega score-based top k seeds vs omega score-based top k seeds all vs all comparison
            sort(sequence1Seeds.begin(), sequence1Seeds.end(), comparePiCellByOmega);
            sort(sequence2Seeds.begin(), sequence2Seeds.end(), comparePiCellByOmega);

            bool bHashCollision = false;

            for(int i = 0; i < numTopPicks; i++) {
                if(sequence1Seeds[i].seed.empty()) {
                    exit(EXIT_FAILURE);
                }

                for(int j = 0; j < numTopPicks; j++) {
                    if(sequence2Seeds[j].seed.empty()) {
                        exit(EXIT_FAILURE);
                    }

                    if(sequence1Seeds[i].seed == sequence2Seeds[j].seed) {
                        numHashCollisions++;
                        bHashCollision = true;
                        break;
                    }
                }

                if(bHashCollision) {
                    break;
                }
            }
        } else {
            bool bHashCollision = false;

            for(int i = 0; i < sequence1Seeds.size(); i++) {
                if(bDoSingleComparison) {
                    // Specific pivot (i, j) vs same pivot (i, j) seeds comparison (following hash collision experiment in SubseqHash2)
                    if(sequence1Seeds[i].seed.empty() || sequence2Seeds[i].seed.empty()) {
                        exit(EXIT_FAILURE);
                    }

                    if(sequence1Seeds[i].seed == sequence2Seeds[i].seed) {
                        numHashCollisions++;
                        break;
                    }
                } else {
                    // All pivots vs all pivots seeds comparison
                    if(sequence1Seeds[i].seed.empty()) {
                        exit(EXIT_FAILURE);
                    }

                    for(int j = 0; j < sequence2Seeds.size(); j++) {
                        if(sequence2Seeds[j].seed.empty()) {
                            exit(EXIT_FAILURE);
                        }

                        if(sequence1Seeds[i].seed == sequence2Seeds[j].seed) {
                            numHashCollisions++;
                            bHashCollision = true;
                            break;
                        }
                    }

                    if(bHashCollision) {
                        break;
                    }
                }
            }
        }

        numTotalSequencePairs++;
    }
    
    sequencePairsFile.close();

    cout << "Hash collision probability = " << (float) numHashCollisions / numTotalSequencePairs << endl;

    return 0;
}