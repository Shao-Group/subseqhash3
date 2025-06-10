#include "subseqhash3.hpp"

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

    int numHashCollisions = 0, numTotalSequencePairs = 0;

    while(sequencePairsFile >> sequence1) {
        sequencePairsFile >> sequence2;

        vector<PiCell> sequence1Seeds = subseqHash3.solvePivotDP(sequence1, sequence1.length()), sequence2Seeds = subseqHash3.solvePivotDP(sequence2, sequence2.length());

        for(int i = 0; i < sequence1Seeds.size(); i++) {
            // [future task] We should look into why we are getting empty seeds with psi and omega values of d and NEG_INF respectively
            if(!sequence1Seeds[i].seedData->seed.empty() && !sequence2Seeds[i].seedData->seed.empty()) {
                if(sequence1Seeds[i].seedData->seed == sequence2Seeds[i].seedData->seed) {
                    numHashCollisions++;
                    break;
                }
            }
        }

        numTotalSequencePairs++;
    }

    sequencePairsFile.close();

    cout << "Hash collision probability = " << numHashCollisions / numTotalSequencePairs << endl;
    
    return 0;
}