#ifndef VFCALC_H
#define VFCALC_H

#include <unordered_map>
#include <map>
#include <vector>

using namespace std;

// board with potential unsolved tiles, 0 = voltorb, 1/2/3 themselves, 4 = unflipped
// you really need to make this presentable
class vfBoard {
public:
    uint8_t tiles[5][5];

    vfBoard(); // TEST REMOVE THIS
};

// Complete vfBoard, indicating a possible board for a specific round. Weight = relative frequency of this board vs others
class vfBoardC : public vfBoard {
public:
    double weight;

    vfBoardC(): vfBoard(), weight(0) {}
};

// In-game, possibly incomplete vfBoard, for usage during calculations
class vfBoardI : public vfBoard {
public:
    size_t hash;

    vfBoardI(): vfBoard(), hash(0) {}
    void manualHash();
    void setTileAndHash(int x, int y, int v);
    bool operator== (const vfBoardI& other) const;
    bool operator< (const vfBoardI& other) const;
};

template<>
struct std::hash<vfBoardI> {
    size_t operator()(const vfBoardI &v) const;
};

struct posData {
    double winChance;
    double gameTime;
    double totalWeight;
};

class VFCalc {
private:
    int level;
    int rowData[2][5][3];       // first axis: [0]: row, [1]: col
                                // 3: score, voltorbs, plusness
                                // Consider "plusness" of a row to be (number of 2s) + 2*(number of 3s) for any valid board
                                // A 4/1 row is +0, 5/1 is +1, 5/2 is +2, 6/3 is +4, 10/1 is +6, etc.

    vector<vfBoardC> possBoards;
    vector<int> possProfiles;

    unordered_map<vfBoardI, posData> T;
    // map<vfBoardI, posData> T;

    // debug info
    int wc_ex, wc_te;

    void getPossProfiles();
    void getAllBoards();
    void getSubBoards(vfBoardI& vt, vector<int>& subBoards);
    double bruteForceWinChances();
    double bruteForceWinChances(vfBoardI& vt);
    double bruteForceWinChances(vfBoardI& curPos, vector<uint16_t>& sub, int plusData[2][5], double tw);
public:
    VFCalc(int level, int rd[2][5][2]);
    void exec();
    bool doable(vfBoardI& vt);
    double doabilityEst(double d, double I, double T);
    void getInDepthData(vfBoardI& vt, double wc[5][5], double gt[5][5]);
    void getHeuristics(vfBoardI& vt, double wh[5][5], double vr[5][5]);
    void display(vfBoardI& vt, bool useInDepthData);
    int getNumBoards();
};

#endif
