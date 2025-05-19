#include<fstream>
#include<iostream>
#include<map>
#include<random>
#include<string>

using namespace std;

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

    void generateTables();
    void loadTables();

public:
    SubseqHash3();
    SubseqHash3(int, int, map<char, int>);
    ~SubseqHash3();
};

void SubseqHash3::generateTables() {

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
    this->loadTables();
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