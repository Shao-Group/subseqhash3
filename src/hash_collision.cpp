#include "subseqhash3.hpp"

int main(int argc, char** argv) {
    SubseqHash3 subseqhash3(5, 7, {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}});

    // subseqhash3.generateTables();
    subseqhash3.loadTables();

    subseqhash3.solvePivotDP("ACGTACGTACGTACGT", 10);
    
    return 0;
}