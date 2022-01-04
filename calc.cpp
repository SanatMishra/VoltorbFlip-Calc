#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <unordered_map>
#include <map>
#include <vector>
#include <limits>

// https://graphics.stanford.edu/~seander/bithacks.html
// iterates over combinations of s bits with x of them set to 1.
#define ITERPERM(s, x, i, t) for(i = (1 << x) - 1; i < 1 << s; t = (i | (i - 1)) + 1, i = (i == 0 ? 1 << s : t | ( ((t & -t)/(i & -i) >> 1) - 1)))
#define MASK(a, b) ( ((1 << b) - 1) ^ ((1 << a) - 1) )

using namespace std;

const double dbl_inf = numeric_limits<double>::infinity();

// https://graphics.stanford.edu/~seander/bithacks.html
int cnt(int x) {
    x = (x & 0x55555555) + (x >> 1 & 0x55555555);
    x = (x & 0x33333333) + (x >> 2 & 0x33333333);
    return ((x + (x >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
}

bool dblEq(double a, double b) {
    return fabs(a - b) < 0.000001;
}

size_t pow5table[5][5] = {
    {1, 5, 25, 125, 625},
    {3125, 15625, 78125, 390625, 1953125},
    {9765625, 48828125, 244140625, 1220703125, 6103515625},
    {30517578125, 152587890625, 762939453125, 3814697265625, 19073486328125},
    {95367431640625, 476837158203125, 2384185791015625, 11920928955078125, 59604644775390625}
};

// board with potential unsolved tiles, 0 = voltorb, 1/2/3 themselves, 4 = unflipped
// you really need to make this presentable
class vfBoard {
public:
    uint8_t tiles[5][5];

    vfBoard() {
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                tiles[i][j] = 4;
    }
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

    void manualHash() {
        hash = 0;
        for (int i = 4; i >= 0; i--)
            for (int j = 4; j >= 0; j--)
                hash = 5*hash + (tiles[i][j] & 3);
    }

    // hashing all 25 tiles becomes expensive, but we can recalculate the hash
    // with a few operations if only one tile is changed (as is usually the case).
    // split into set/unset if this is still significant
    void setTileAndHash(int x, int y, int v) {
        hash += pow5table[x][y]*((v & 3) - (tiles[x][y] & 3));
        tiles[x][y] = v;
    }

    bool operator== (const vfBoardI& other) const {
        return hash == other.hash;
    }

    bool operator< (const vfBoardI& other) const {
        return hash < other.hash;
    }
};

template<>
struct std::hash<vfBoardI> {
    size_t operator()(const vfBoardI &v) const {
        return v.hash;
    }
};

struct posData {
    double winChance;
    double gameTime;
    double totalWeight;
};

class VFCalc {

public: // let these remain public until test code requiring them is gone
    int level;
    int rowData[2][5][3];       // first axis: [0]: row, [1]: col
                                // 3: score, voltorbs, plusness
                                // Consider "plusness" of a row to be (number of 2s) + 2*(number of 3s) for any valid board
                                // A 4/1 row is +0, 5/1 is +1, 5/2 is +2, 6/3 is +4, 10/1 is +6, etc.

    static int boardProfiles[8][5][3];

    vector<vfBoardC> possBoards;
    vector<int> possProfiles;
    
    unordered_map<vfBoardI, posData> T;
    // map<vfBoardI, posData> T;

// public:
    VFCalc(int level, int rd[2][5][2]) {
        this->level = level;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 5; j++) {
                rowData[i][j][0] = rd[i][j][0];
                rowData[i][j][1] = rd[i][j][1];
                rowData[i][j][2] = rd[i][j][1] + rd[i][j][0] - 5;
            }
        }

        getPossibleBoards();
        cout << possBoards.size() << " possible boards" << endl;
    }

    void getPossibleBoards() {

        int orbs = 0;
        for (int i = 0; i < 5; i++) {
            orbs += rowData[0][i][1];
        }
        int totalPlusness = 0;
        for (int i = 0; i < 5; i++) {
            totalPlusness += rowData[0][i][2];
        }
        for (int i = 0; i < 5; i++) {
            if (boardProfiles[level - 1][i][0] + 2*boardProfiles[level - 1][i][1] == totalPlusness
                  && boardProfiles[level - 1][i][2] == orbs) {
                possProfiles.push_back(i);
            }
        }

        int pbIndex = 0; // tracks size of possBoards before next loop iteration begins, so we know which boards to assign weights to

        for(int p : possProfiles) {
            vfBoardC cur;
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++)
                    cur.tiles[i][j] = 1;

            int cnt3 = boardProfiles[level - 1][p][1];
            int cnt2 = boardProfiles[level - 1][p][0];
            int cntv = boardProfiles[level - 1][p][2];
            int cnt1 = 25 - cnt3 - cnt2 - cntv;

            int rd2[2][5], v2[2][5];
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 5; j++) {
                    rd2[i][j] = rowData[i][j][2];
                    v2[i][j] = rowData[i][j][1];
                }
            }

            int sp = 0; // total space for 3s
            for (int i = 0; i < 5; i++) sp += rd2[0][i]/2;

            int ind[6] = {0}; // running sum of sp
            for (int j = 1; j <= 5; j++) {
                ind[j] = ind[j - 1] + rd2[0][j - 1]/2;
            }

            int k[5]; //3s
            int m[5];
            vector<int> nzi[5];

            int r[5]; // 2s
            vector<int> oi[5];

            int s[5]; // voltorbs
            vector<int> pi[5];

            int i;
            int tp; // used for permutations

            int variants = 0; // number of boards for this profile, used for weights at the end

            // insert tiles in the order 3, 2, v, from most to least restrictive.
            // Probably unnecessarily complex but I thought this would have significant runtime
            // For 3s we must account for boards that do not fill up all potential spaces with 3s
            // There might be a boardProfile with 4 3s, even though the board itself can fit 5.
            // We use k instead of rd2, and iterate over possible arrangements of k over the rows
            ITERPERM(sp, cnt3, i, tp) {
                for (int j = 0; j < 5; j++) {
                    k[j] = cnt(i & MASK(ind[j], ind[j + 1]));
                }

                // Is there really no way to simplify this
                nzi[0].clear();
                for (int n = 0; n < 5; n++)
                    if (rd2[1][n]/2 > 0)
                        nzi[0].push_back(n);
                ITERPERM(nzi[0].size(), k[0], m[0], tp) {
                    for (int n = 0; n < nzi[0].size(); n++) {
                        if (m[0] & (1 << n)) {
                            cur.tiles[nzi[0][n]][0] = 3;
                            rd2[0][0] -= 2;
                            rd2[1][nzi[0][n]] -= 2;
                        }
                    }

                    nzi[1].clear();
                    for (int n = 0; n < 5; n++)
                        if (rd2[1][n]/2 > 0)
                            nzi[1].push_back(n);
                    ITERPERM(nzi[1].size(), k[1], m[1], tp) {
                        for (int n = 0; n < nzi[1].size(); n++) {
                            if (m[1] & (1 << n)) {
                                cur.tiles[nzi[1][n]][1] = 3;
                                rd2[0][1] -= 2;
                                rd2[1][nzi[1][n]] -= 2;
                            }
                        }

                        nzi[2].clear();
                        for (int n = 0; n < 5; n++)
                            if (rd2[1][n]/2 > 0)
                                nzi[2].push_back(n);
                        ITERPERM(nzi[2].size(), k[2], m[2], tp) {
                            for (int n = 0; n < nzi[2].size(); n++) {
                                if (m[2] & (1 << n)) {
                                    cur.tiles[nzi[2][n]][2] = 3;
                                    rd2[0][2] -= 2;
                                    rd2[1][nzi[2][n]] -= 2;
                                }
                            }

                            nzi[3].clear();
                            for (int n = 0; n < 5; n++)
                                if (rd2[1][n]/2 > 0)
                                    nzi[3].push_back(n);
                            ITERPERM(nzi[3].size(), k[3], m[3], tp) {
                                for (int n = 0; n < nzi[3].size(); n++) {
                                    if (m[3] & (1 << n)) {
                                        cur.tiles[nzi[3][n]][3] = 3;
                                        rd2[0][3] -= 2;
                                        rd2[1][nzi[3][n]] -= 2;
                                    }
                                }

                                nzi[4].clear();
                                for (int n = 0; n < 5; n++)
                                    if (rd2[1][n]/2 > 0)
                                        nzi[4].push_back(n);
                                ITERPERM(nzi[4].size(), k[4], m[4], tp) {
                                    for (int n = 0; n < nzi[4].size(); n++) {
                                        if (m[4] & (1 << n)) {
                                            cur.tiles[nzi[4][n]][4] = 3;
                                            rd2[0][4] -= 2;
                                            rd2[1][nzi[4][n]] -= 2;
                                        }
                                    }

                                    // 3s inserted, insert 2s
                                    oi[0].clear();
                                    for (int n = 0; n < 5; n++)
                                        if (cur.tiles[n][0] == 1 && rd2[1][n] > 0)
                                            oi[0].push_back(n);
                                    if (oi[0].size() < rd2[0][0])
                                        goto L1;
                                    ITERPERM(oi[0].size(), rd2[0][0], r[0], tp) {
                                        for (int n = 0; n < oi[0].size(); n++) {
                                            if (r[0] & (1 << n)) {
                                                cur.tiles[oi[0][n]][0] = 2;
                                                rd2[0][0] -= 1;
                                                rd2[1][oi[0][n]] -= 1;
                                            }
                                        }

                                        oi[1].clear();
                                        for (int n = 0; n < 5; n++)
                                            if (cur.tiles[n][1] == 1 && rd2[1][n] > 0)
                                                oi[1].push_back(n);
                                        if (oi[1].size() < rd2[0][1])
                                            goto L2;
                                        ITERPERM(oi[1].size(), rd2[0][1], r[1], tp) {
                                            for (int n = 0; n < oi[1].size(); n++) {
                                                if (r[1] & (1 << n)) {
                                                    cur.tiles[oi[1][n]][1] = 2;
                                                    rd2[0][1] -= 1;
                                                    rd2[1][oi[1][n]] -= 1;
                                                }
                                            }

                                            oi[2].clear();
                                            for (int n = 0; n < 5; n++)
                                                if (cur.tiles[n][2] == 1 && rd2[1][n] > 0)
                                                    oi[2].push_back(n);
                                            if (oi[2].size() < rd2[0][2])
                                                goto L3;
                                            ITERPERM(oi[2].size(), rd2[0][2], r[2], tp) {
                                                for (int n = 0; n < oi[2].size(); n++) {
                                                    if (r[2] & (1 << n)) {
                                                        cur.tiles[oi[2][n]][2] = 2;
                                                        rd2[0][2] -= 1;
                                                        rd2[1][oi[2][n]] -= 1;
                                                    }
                                                }

                                                oi[3].clear();
                                                for (int n = 0; n < 5; n++)
                                                    if (cur.tiles[n][3] == 1 && rd2[1][n] > 0)
                                                        oi[3].push_back(n);
                                                if (oi[3].size() < rd2[0][3])
                                                    goto L4;
                                                ITERPERM(oi[3].size(), rd2[0][3], r[3], tp) {
                                                    for (int n = 0; n < oi[3].size(); n++) {
                                                        if (r[3] & (1 << n)) {
                                                            cur.tiles[oi[3][n]][3] = 2;
                                                            rd2[0][3] -= 1;
                                                            rd2[1][oi[3][n]] -= 1;
                                                        }
                                                    }

                                                    oi[4].clear();
                                                    for (int n = 0; n < 5; n++)
                                                        if (cur.tiles[n][4] == 1 && rd2[1][n] > 0)
                                                            oi[4].push_back(n);
                                                    if (oi[4].size() < rd2[0][4])
                                                        goto L5;
                                                    ITERPERM(oi[4].size(), rd2[0][4], r[4], tp) {
                                                        for (int n = 0; n < oi[4].size(); n++) {
                                                            if (r[4] & (1 << n)) {
                                                                cur.tiles[oi[4][n]][4] = 2;
                                                                rd2[0][4] -= 1;
                                                                rd2[1][oi[4][n]] -= 1;
                                                            }
                                                        }

                                                        // 2s inserted, insert voltorbs
                                                        pi[0].clear();
                                                        for (int n = 0; n < 5; n++)
                                                            if (cur.tiles[n][0] == 1 && v2[1][n] > 0)
                                                                pi[0].push_back(n);
                                                        if (pi[0].size() < v2[0][0])
                                                            goto L6;
                                                        ITERPERM(pi[0].size(), v2[0][0], s[0], tp) {
                                                            for (int n = 0; n < pi[0].size(); n++) {
                                                                if (s[0] & (1 << n)) {
                                                                    cur.tiles[pi[0][n]][0] = 0;
                                                                    v2[0][0] -= 1;
                                                                    v2[1][pi[0][n]] -= 1;
                                                                }
                                                            }

                                                            pi[1].clear();
                                                            for (int n = 0; n < 5; n++)
                                                                if (cur.tiles[n][1] == 1 && v2[1][n] > 0)
                                                                    pi[1].push_back(n);
                                                            if (pi[1].size() < v2[0][1])
                                                                goto L7;
                                                            ITERPERM(pi[1].size(), v2[0][1], s[1], tp) {
                                                                for (int n = 0; n < pi[1].size(); n++) {
                                                                    if (s[1] & (1 << n)) {
                                                                        cur.tiles[pi[1][n]][1] = 0;
                                                                        v2[0][1] -= 1;
                                                                        v2[1][pi[1][n]] -= 1;
                                                                    }
                                                                }

                                                                pi[2].clear();
                                                                for (int n = 0; n < 5; n++)
                                                                    if (cur.tiles[n][2] == 1 && v2[1][n] > 0)
                                                                        pi[2].push_back(n);
                                                                if (pi[2].size() < v2[0][2])
                                                                    goto L8;
                                                                ITERPERM(pi[2].size(), v2[0][2], s[2], tp) {
                                                                    for (int n = 0; n < pi[2].size(); n++) {
                                                                        if (s[2] & (1 << n)) {
                                                                            cur.tiles[pi[2][n]][2] = 0;
                                                                            v2[0][2] -= 1;
                                                                            v2[1][pi[2][n]] -= 1;
                                                                        }
                                                                    }

                                                                    pi[3].clear();
                                                                    for (int n = 0; n < 5; n++)
                                                                        if (cur.tiles[n][3] == 1 && v2[1][n] > 0)
                                                                            pi[3].push_back(n);
                                                                    if (pi[3].size() < v2[0][3])
                                                                        goto L9;
                                                                    ITERPERM(pi[3].size(), v2[0][3], s[3], tp) {
                                                                        for (int n = 0; n < pi[3].size(); n++) {
                                                                            if (s[3] & (1 << n)) {
                                                                                cur.tiles[pi[3][n]][3] = 0;
                                                                                v2[0][3] -= 1;
                                                                                v2[1][pi[3][n]] -= 1;
                                                                            }
                                                                        }

                                                                        pi[4].clear();
                                                                        for (int n = 0; n < 5; n++)
                                                                            if (cur.tiles[n][4] == 1 && v2[1][n] > 0)
                                                                                pi[4].push_back(n);
                                                                        if (pi[4].size() < v2[0][4])
                                                                            goto L10;
                                                                        ITERPERM(pi[4].size(), v2[0][4], s[4], tp) {
                                                                            for (int n = 0; n < pi[4].size(); n++) {
                                                                                if (s[4] & (1 << n)) {
                                                                                    cur.tiles[pi[4][n]][4] = 0;
                                                                                    v2[0][4] -= 1;
                                                                                    v2[1][pi[4][n]] -= 1;
                                                                                }
                                                                            }

                                                                            // we have a board
                                                                            possBoards.push_back(cur);
                                                                            variants++;

                                                                            for (int n = 0; n < pi[4].size(); n++) {
                                                                                if (s[4] & (1 << n)) {
                                                                                    cur.tiles[pi[4][n]][4] = 1;
                                                                                    v2[0][4] += 1;
                                                                                    v2[1][pi[4][n]] += 1;
                                                                                }
                                                                            }
                                                                        }
                                                                        L10:
                                                                        for (int n = 0; n < pi[3].size(); n++) {
                                                                            if (s[3] & (1 << n)) {
                                                                                cur.tiles[pi[3][n]][3] = 1;
                                                                                v2[0][3] += 1;
                                                                                v2[1][pi[3][n]] += 1;
                                                                            }
                                                                        }
                                                                    }
                                                                    L9:
                                                                    for (int n = 0; n < pi[2].size(); n++) {
                                                                        if (s[2] & (1 << n)) {
                                                                            cur.tiles[pi[2][n]][2] = 1;
                                                                            v2[0][2] += 1;
                                                                            v2[1][pi[2][n]] += 1;
                                                                        }
                                                                    }
                                                                }
                                                                L8:
                                                                for (int n = 0; n < pi[1].size(); n++) {
                                                                    if (s[1] & (1 << n)) {
                                                                        cur.tiles[pi[1][n]][1] = 1;
                                                                        v2[0][1] += 1;
                                                                        v2[1][pi[1][n]] += 1;
                                                                    }
                                                                }
                                                            }
                                                            L7:
                                                            for (int n = 0; n < pi[0].size(); n++) {
                                                                if (s[0] & (1 << n)) {
                                                                    cur.tiles[pi[0][n]][0] = 1;
                                                                    v2[0][0] += 1;
                                                                    v2[1][pi[0][n]] += 1;
                                                                }
                                                            }
                                                        }
                                                        L6:
                                                        for (int n = 0; n < oi[4].size(); n++) {
                                                            if (r[4] & (1 << n)) {
                                                                cur.tiles[oi[4][n]][4] = 1;
                                                                rd2[0][4] += 1;
                                                                rd2[1][oi[4][n]] += 1;
                                                            }
                                                        }
                                                    }
                                                    L5:
                                                    for (int n = 0; n < oi[3].size(); n++) {
                                                        if (r[3] & (1 << n)) {
                                                            cur.tiles[oi[3][n]][3] = 1;
                                                            rd2[0][3] += 1;
                                                            rd2[1][oi[3][n]] += 1;
                                                        }
                                                    }
                                                }
                                                L4:
                                                for (int n = 0; n < oi[2].size(); n++) {
                                                    if (r[2] & (1 << n)) {
                                                        cur.tiles[oi[2][n]][2] = 1;
                                                        rd2[0][2] += 1;
                                                        rd2[1][oi[2][n]] += 1;
                                                    }
                                                }
                                            }
                                            L3:
                                            for (int n = 0; n < oi[1].size(); n++) {
                                                if (r[1] & (1 << n)) {
                                                    cur.tiles[oi[1][n]][1] = 1;
                                                    rd2[0][1] += 1;
                                                    rd2[1][oi[1][n]] += 1;
                                                }
                                            }
                                        }
                                        L2:
                                        for (int n = 0; n < oi[0].size(); n++) {
                                            if (r[0] & (1 << n)) {
                                                cur.tiles[oi[0][n]][0] = 1;
                                                rd2[0][0] += 1;
                                                rd2[1][oi[0][n]] += 1;
                                            }
                                        }
                                    }

                                    L1:
                                    for (int n = 0; n < nzi[4].size(); n++) {
                                        if (m[4] & (1 << n)) {
                                            cur.tiles[nzi[4][n]][4] = 1;
                                            rd2[0][4] += 2;
                                            rd2[1][nzi[4][n]] += 2;
                                        }
                                    }
                                }
                                for (int n = 0; n < nzi[3].size(); n++) {
                                    if (m[3] & (1 << n)) {
                                        cur.tiles[nzi[3][n]][3] = 1;
                                        rd2[0][3] += 2;
                                        rd2[1][nzi[3][n]] += 2;
                                    }
                                }
                            }
                            for (int n = 0; n < nzi[2].size(); n++) {
                                if (m[2] & (1 << n)) {
                                    cur.tiles[nzi[2][n]][2] = 1;
                                    rd2[0][2] += 2;
                                    rd2[1][nzi[2][n]] += 2;
                                }
                            }
                        }
                        for (int n = 0; n < nzi[1].size(); n++) {
                            if (m[1] & (1 << n)) {
                                cur.tiles[nzi[1][n]][1] = 1;
                                rd2[0][1] += 2;
                                rd2[1][nzi[1][n]] += 2;
                            }
                        }
                    }
                    for (int n = 0; n < nzi[0].size(); n++) {
                        if (m[0] & (1 << n)) {
                            cur.tiles[nzi[0][n]][0] = 1;
                            rd2[0][0] += 2;
                            rd2[1][nzi[0][n]] += 2;
                        }
                    }
                }
            }

            for (int ind = pbIndex; ind < pbIndex + variants; ind++)
                possBoards[ind].weight = 1.0/variants;
            pbIndex += variants;

            // cout << p << ": " << variants << endl;

        }
    }

    void getSubBoards(vfBoardI& vt, vector<int>& subBoards) {
        for (int i = 0; i < possBoards.size(); i++) {
            for (int j = 0; j < 5; j++)
                for (int k = 0; k < 5; k++)
                    if (vt.tiles[j][k] != 4 && possBoards[i].tiles[j][k] != vt.tiles[j][k])
                        goto gsb_ENDLOOP;
            subBoards.push_back(i);
            // cout << endl;
            // for (int j = 0; j < 5; j++) {
            //     for (int k = 0; k < 5; k++) {
            //         cout << (int)(possBoards[i].tiles[k][j]) << " ";
            //     }
            //     cout << endl;
            // }
            // cout << endl;

            gsb_ENDLOOP: ;
        }
    }

    void display(vfBoardI& vt, bool useWinData) {
        double wh[5][5], vr[5][5], wc[5][5], gt[5][5];
        int x, y;
        getHeuristics(vt, wh, vr);
        getInDepthData(vt, wc, gt);

        cout << endl;
        for (int j = 0; j < 5; j++) {
            for (int i = 0; i < 5; i++) {
                cout << "  " << (int)(vt.tiles[i][j]) << "   ";
            }
            cout << "| " << rowData[0][j][0] << " " << rowData[0][j][1] << " ";

            for (int i = 0; i < 5; i++) {
                if (useWinData) {
                    if (wc[i][j] == -1)
                        printf("     ?     |");
                    else
                        printf("%.3f %.3f|", wc[i][j], wh[i][j]);
                } else {
                    printf("%.3f %.3f|", wh[i][j], vr[i][j]);
                }
            }
            cout << endl;
        }

        cout << "------------------------------+" << endl;
        for (int i = 0; i < 5; i++)
            cout << "  " << rowData[1][i][0] << "   ";
        cout << endl;
        for (int i = 0; i < 5; i++)
            cout << "  " << rowData[1][i][1] << "   ";
        cout << endl << endl;

        if (useWinData) {
            x = 0; y = 0;
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 5; j++) {
                    if (wc[i][j] > wc[x][y] || dblEq(wc[i][j], wc[x][y]) &&
                        (gt[i][j] < gt[x][y] || dblEq(gt[i][j], gt[x][y]) &&
                        (wh[i][j] > wh[x][y] || dblEq(wh[i][j], wh[x][y]) &&
                        vr[i][j] < vr[x][y]))) {
                        x = i; y = j;
                    }
                }
            }
            cout << x << " " << y << " " << wc[x][y] << " " << wh[x][y] << endl;
        } else {
            x = 0; y = 0;
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 5; j++) {
                    if (wh[i][j] > wh[x][y] || dblEq(wh[i][j], wh[x][y]) && vr[i][j] < vr[x][y]) {
                        x = i; y = j;
                    }
                }
            }
            cout << x << " " << y << " " << wh[x][y] << " " << vr[x][y] << endl;
        }
    }

    bool doable(vfBoardI& vt) {
        vector<int> sub;
        getSubBoards(vt, sub);

        if (sub.size() >= 15000) return 0; // arbitrary decision

        bool poss[5][5][3] = {{{0}}};

        double d = 0;
        double I = 0;
        double T = 0;

        for (int p : possProfiles) {
            d += boardProfiles[level - 1][p][0] + boardProfiles[level - 1][p][1];
        }
        d /= possProfiles.size();

        int rd2[2][5];
        for (int i = 0; i < 5; i++) {
            rd2[0][i] = rowData[0][i][2];
            rd2[1][i] = rowData[1][i][2];
        }
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                if (vt.tiles[i][j] == 2 || vt.tiles[i][j] == 3) {
                    rd2[0][j] -= vt.tiles[i][j] - 1;
                    rd2[1][i] -= vt.tiles[i][j] - 1;
                    d--;
                }
            }
        }

        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                if (vt.tiles[i][j] == 4 && rd2[0][j] != 0 && rd2[1][i] != 0) {
                    T++;
                    for (auto b : possBoards) {
                        if (b.tiles[i][j] > 0) {
                            if (!poss[i][j][b.tiles[i][j] - 1])
                                I++;
                            poss[i][j][b.tiles[i][j] - 1] = 1;
                        }
                    }
                }
            }
        }

        cout << I << " " << T << " " << d << " " << doabilityEst(d, I, T) << endl;

        return doabilityEst(d, I, T) < pow(10, 29);//doabilityEst(8.5, 45, 20); why does 10^29 work at all
        // return (tgamma(d + 1)*pow(I/d, d) < 542097412.116);
    }

    double doabilityEst(double d, double I, double T) {
        return d*pow(I, (I/T - 1)*(d + 1));
    }

    bool doable_old(vfBoardI& vt) {
        vector<int> sub;
        getSubBoards(vt, sub);
        if (sub.size() >= 15000) return 0; // arbitrary decision

        bool poss[5] = {0};
        int orbs = 0;
        for (int i = 0; i < 5; i++) {
            orbs += rowData[0][i][1];
        }
        int totalPlusness = 0;
        for (int i = 0; i < 5; i++) {
            totalPlusness += rowData[0][i][2];
        }

        int found2 = 0, found3 = 0;

        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                if (vt.tiles[i][j] == 2)
                    found2++;
                if (vt.tiles[i][j] == 3)
                    found3++;
            }
        }

        int worstDepth = 0;

        for (int i = 0; i < 5; i++) {
            if (boardProfiles[level - 1][i][0] + 2*boardProfiles[level - 1][i][1] == totalPlusness
                  && boardProfiles[level - 1][i][2] == orbs
                  && boardProfiles[level - 1][i][0] >= found2 && boardProfiles[level - 1][i][1] >= found3) {

                if (boardProfiles[level - 1][i][0] + boardProfiles[level - 1][i][1] - found2 - found3 > worstDepth) {
                    worstDepth = boardProfiles[level - 1][i][0] + boardProfiles[level - 1][i][1] - found2 - found3;
                }
            }
        }

        int validTiles = 0;

        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                if (rowData[0][j][2] != 0 && rowData[0][j][2] != 0)
                    validTiles++;
            }
        }

        cout << "v = " << validTiles << ", d = " << worstDepth << endl;

        // rough estimate of complexity is (d!*(v/d)^d)^2
        // d = depth, assumed proportional to total 2s and 3s. To keep it simple say its 1.5x that, it may depend on v/1 layout
        // v = amount of valid tiles (not 4/1 etc)
        // so we check the value of d!*(v/d)^d

        double mv = 0.7;
        double md = 1.5;

        double d = md*worstDepth;
        double v = mv*validTiles;
        double x = tgamma(d + 1)*pow(v/d, d);

        cout << tgamma(d + 1) << " " << pow(v/d, d) << " " << x << endl;

        return x < pow(10, 13); // ?
    }

    void exec() {
        vfBoardI vt;

        double wh[5][5], vr[5][5];
        int x, y, v;

        while (!doable(vt)) {
            display(vt, 0);
            cin >> x >> y >> v;
            vt.tiles[x][y] = v;
        }

        bool h = 1;
        cout << "Brute-forcing now..." << endl;
        vt.manualHash();
        while(T.find(vt) == T.end() || T[vt].gameTime != 0) {
            if (T.find(vt) == T.end()) {
                cout << "-- NEW POSITION --" << endl;
                bruteForceWinChances(vt);
            }
            display(vt, 1);
            cin >> x >> y >> v;
            vt.setTileAndHash(x, y, v);
        }

    }

    // debug info
    int wc_ex, wc_te;

    double bruteForceWinChances() {
        vfBoardI vt;
        return bruteForceWinChances(vt);
    }

    double bruteForceWinChances(vfBoardI& vt) {
        vector<int> subi;
        getSubBoards(vt, subi);

        vector<uint16_t> sub;
        double tw = 0;
        for (auto x : subi) {
            sub.push_back(x);
            tw += possBoards[x].weight;
        }

        int plusData[2][5];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 5; j++)
                plusData[i][j] = rowData[i][j][2];

        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                if (vt.tiles[i][j] == 2 || vt.tiles[i][j] == 3) {
                    plusData[0][j] -= vt.tiles[i][j] - 1;
                    plusData[1][i] -= vt.tiles[i][j] - 1;
                }
            }
        }

        vt.manualHash();

        wc_ex = 0;
        wc_te = 0;
        double ret = bruteForceWinChances(vt, sub, plusData, tw);
        cout << "E " << wc_ex << endl;
        cout << "T " << wc_te << endl;
        return ret;
    }

    // todo: take 0 voltorb rows immediately
    //       get win chance for all tiles
    double bruteForceWinChances(vfBoardI& curPos, vector<uint16_t>& sub, int plusData[2][5], double tw) {
        wc_ex++;

        if (T.find(curPos) != T.end()) return T[curPos].winChance;
        // do we need this?
        // gprof later, maybe initial solved check thats faster than the loop
        // if (plusData[0][0] == 0 && plusData[0][1] == 0 && plusData[0][2] == 0 && plusData[0][3] == 0 && plusData[0][4] == 0) {
        //     T[curPos] = {5, 5, 1, 0, (int)(sub.size())};
        //     return 1;
        // }

        double ret = 0, retd = 25;
        bool ran = 0;
        for (int x = 0; x < 5; x++) {
            for (int y = 0; y < 5; y++) {
                if (curPos.tiles[x][y] != 4 || plusData[0][y] == 0 || plusData[1][x] == 0)
                    continue;
                ran = 1;

                double ret2 = 0, retd2 = 0, w[3] = {0};
                vector<uint16_t> subn[3];

                // skip count if curPos [xy]=i is solved --> get count from child calls?
                // does that require knowing whether solved or not before entering loop
                for (int i = 0; i < sub.size(); i++) {
                    if (possBoards[sub[i]].tiles[x][y] > 0) {
                        subn[possBoards[sub[i]].tiles[x][y] - 1].push_back(sub[i]);
                        w[possBoards[sub[i]].tiles[x][y] - 1] += possBoards[sub[i]].weight;
                    }
                }
                for (int i = 0; i < 3; i++) {
                    if (subn[i].size() > 0) {
                        // maybe conglomerate this @ some point
                        plusData[0][y] -= i + 1 - (curPos.tiles[x][y] == 4 ? 1 : curPos.tiles[x][y]);
                        plusData[1][x] -= i + 1 - (curPos.tiles[x][y] == 4 ? 1 : curPos.tiles[x][y]);
                        // plusData[0][y] += (curPos.tiles[x][y] - 1) % 3 - i % 3;
                        // plusData[1][x] += (curPos.tiles[x][y] - 1) % 3 - i % 3;
                        curPos.setTileAndHash(x, y, i + 1);
                        ret2 += w[i]*bruteForceWinChances(curPos, subn[i], plusData, w[i]);
                        // T[curPos] should always exist
                        retd2 += w[i]*T[curPos].gameTime;
                    }
                }
                ret2 /= tw;
                retd2 = 1 + retd2/tw;

                plusData[0][y] += (curPos.tiles[x][y] == 4 ? 1 : curPos.tiles[x][y]) - 1;
                plusData[1][x] += (curPos.tiles[x][y] == 4 ? 1 : curPos.tiles[x][y]) - 1;
                curPos.setTileAndHash(x, y, 4);
                if (ret2 > ret || dblEq(ret2, ret) && retd2 < retd) {
                    ret = ret2;
                    retd = retd2;
                }
            }
        }

        wc_te++;
        T[curPos] = {ran ? ret : 1, ran ? retd : 0, tw};
        return (ran ? ret : 1);
    }

    void getInDepthData(vfBoardI& vt, double wc[5][5], double gt[5][5]) {
        if (T.find(vt) == T.end()) {
            for (int x = 0; x < 5; x++) {
                for (int y = 0; y < 5; y++) {
                    wc[x][y] = gt[x][y] = -1;
                }
            }
            return;
        }
        for (int x = 0; x < 5; x++) {
            for (int y = 0; y < 5; y++) {
                if (vt.tiles[x][y] != 4) {
                    wc[x][y] = gt[x][y] = -1;
                    continue;
                }

                wc[x][y] = gt[x][y] = 0;
                int temp = vt.tiles[x][y];
                for (int i = 1; i <= 3; i++) {
                    vt.setTileAndHash(x, y, i);
                    if (T.find(vt) != T.end()) {
                        wc[x][y] += T[vt].totalWeight*T[vt].winChance;
                        gt[x][y] += T[vt].totalWeight*T[vt].gameTime;
                    }
                }
                vt.setTileAndHash(x, y, temp);

                if (wc[x][y] == 0) {
                    wc[x][y] = -1; gt[x][y] = -1;
                    continue;
                }

                wc[x][y] /= T[vt].totalWeight;
                gt[x][y] = 1 + gt[x][y]/T[vt].totalWeight;
            }
        }

    }

    // maximize ln(p_v)/s, p_v is chance of voltorb, s is plusness of tile if not voltorb (2*p_3 + p_2)/(1 - p_v)
    void getHeuristics(vfBoardI& vt, double wh[5][5], double vr[5][5]) {
        vector<int> subBoards;
        double tw = 0;
        double den[5][5];
        for (int j = 0; j < 5; j++) {
            for (int k = 0; k < 5; k++) {
                wh[j][k] = 0;
                den[j][k] = 0;
                vr[j][k] = 0;
            }
        }

        getSubBoards(vt, subBoards);
        for (auto x : subBoards)
            tw += possBoards[x].weight;
        // cout << subBoards.size() << endl;

        if (subBoards.size() == 0 || tw == 0) {
            cout << "?????" << endl;
            return;
        }

        for (int i = 0; i < subBoards.size(); i++) {
            for (int j = 0; j < 5; j++) {
                for (int k = 0; k < 5; k++) {
                    int m = possBoards[subBoards[i]].tiles[j][k];
                    double w = possBoards[subBoards[i]].weight;
                    if (m == 2) den[j][k] += w;
                    if (m == 3) den[j][k] += 2*w;
                    if (m == 0) vr[j][k] += w;
                }
            }
            // cout << endl;
            // for (int k = 0; k < 5; k++) {
            //     for (int j = 0; j < 5; j++) {
            //         cout << den[j][k] << "," << pv[j][k] << " ";
            //     }
            //     cout << endl;
            // }
            // cout << endl;
        }

        for (int j = 0; j < 5; j++) {
            for (int k = 0; k < 5; k++) {
                den[j][k] /= tw;
                vr[j][k] /= tw;
                // cout << j << " " << k << " " << "den: " << den[j][k] << " pv: " << pv[j][k] << " " << mx << " " << my << " " << best << endl;
                if (vt.tiles[j][k] != 4 || den[j][k] == 0 || dblEq(vr[j][k], tw)) {
                    wh[j][k] = -dbl_inf;
                    continue;
                }
                wh[j][k] = (1 - vr[j][k])*log(1 - vr[j][k])/den[j][k];
            }
        }
        // cout << mx << " " << my << " " << best << endl;


    }

};

int VFCalc::boardProfiles[8][5][3] = {
    {{3, 1, 6}, // 2s, 3s, voltorbs
     {0, 3, 6},
     {5, 0, 6},
     {2, 2, 6},
     {4, 1, 6}
    },
    {{1, 3, 7},
     {6, 0, 7},
     {3, 2, 7},
     {0, 4, 7},
     {5, 1, 7}
    },
    {{2, 3, 8},
     {7, 0, 8},
     {4, 2, 8},
     {1, 4, 8},
     {6, 1, 8}
    },
    {{3, 3, 8},
     {0, 5, 8},
     {8, 0, 10},
     {5, 2, 10},
     {2, 4, 10}
    },
    {{7, 1, 10},
     {4, 3, 10},
     {1, 5, 10},
     {9, 0, 10},
     {6, 2, 10}
    },
    {{3, 4, 10},
     {0, 6, 10},
     {8, 1, 10},
     {5, 3, 10},
     {2, 5, 10}
    },
    {{7, 2, 10},
     {4, 4, 10},
     {1, 6, 13},
     {9, 1, 13},
     {6, 3, 10}
    },
    {{0, 7, 10},
     {8, 2, 10},
     {5, 4, 10},
     {2, 6, 10},
     {7, 3, 10}
    }
};

void runBoards() {
    unsigned seed = time(0); // 16412511
    cout << seed << endl;
    srand(seed);

    int rd[2][5][2];
    int lvl = 8, profn = 0, prof[3], board[5][5];
    int bcnt = 0;
    while (1) {
        // cout << "Level " << lvl << ", profile " << profn << endl;
        prof[0] = VFCalc::boardProfiles[lvl - 1][profn][0];
        prof[1] = VFCalc::boardProfiles[lvl - 1][profn][1];
        prof[2] = VFCalc::boardProfiles[lvl - 1][profn][2];
        int i = 0;

        for (int j = 0; j < prof[0]; j++)
            *((int*)board + (i++)) = 2;
        for (int j = 0; j < prof[1]; j++)
            *((int*)board + (i++)) = 3;
        for (int j = 0; j < prof[2]; j++)
            *((int*)board + (i++)) = 0;
        for (int j = i; j < 25; j++)
            *((int*)board + (i++)) = 1;

        for (int j = 0; j < 25; j++) {
            int k = rand() % 25;
            int temp = *((int*)board + j);
            *((int*)board + j) = *((int*)board + k);
            *((int*)board + k) = temp;
        }

        for (int j = 0; j < 5; j++) {
            rd[0][j][0] = rd[0][j][1] = 0;
            rd[1][j][0] = rd[1][j][1] = 0;
        }

        for (int j = 0; j < 5; j++) {
            for (int k = 0; k < 5; k++) {
                // cout << board[j][k] << " ";
                rd[1][k][0] += board[j][k];
                rd[0][j][0] += board[j][k];
                rd[0][j][1] += (board[j][k] == 0);
                rd[1][k][1] += (board[j][k] == 0);
            }
            // cout << endl;
        }

        // for (int i = 0; i < 2; i++) {
        //     for (int j = 0; j < 5; j++)
        //         cout << rd[i][j][0] + rd[i][j][1] - 5 << " ";
        //     cout << endl;
        // }
        // for (int i = 0; i < 2; i++)
        //     for (int j = 0; j < 5; j++)
        //         cout << rd[i][j][0] << " " << rd[i][j][1] << " ";
        // cout << endl;

        VFCalc calc(lvl, rd);
        bcnt = calc.possBoards.size();

        vfBoardI vt;

        bool d = calc.doable(vt);
        // cout << d << endl;
        if (!d) {
            // calc.bruteForceWinChances();
            cout << "Level " << lvl << endl;
            calc.exec();
            cout << endl;
        }

        if (bcnt >= 65536) break;
        profn++;
        if (profn >= 5) {
            profn = 0;
            lvl++;
            if (lvl >= 9)
                lvl = 5;
        }
    }

    // calc.exec();

    for (int j = 0; j < 5; j++) {
        for (int k = 0; k < 5; k++) {
            cout << board[k][j] << " ";
        }
        cout << endl;
    }

    cout << lvl << " " << profn << endl;
}

void testHash() {
    vfBoardI vt;
    cout << vt.hash << endl;

    vt.setTileAndHash(2, 3, 3);
    cout << vt.hash << endl;
    vt.manualHash();
    cout << vt.hash << endl;

    vt.setTileAndHash(1, 4, 1);
    cout << vt.hash << endl;
    vt.manualHash();
    cout << vt.hash << endl;

    vt.setTileAndHash(3, 0, 2);
    cout << vt.hash << endl;
    vt.manualHash();
    cout << vt.hash << endl;

    vt.setTileAndHash(1, 4, 4);
    cout << vt.hash << endl;
    vt.manualHash();
    cout << vt.hash << endl;
}

int main() {

    // testHash();
    // return 0;
    runBoards();
    return 0;

    int level = 8;
    int rd[2][5][2] = { { {7, 2}, {6, 1}, {3, 3}, {6, 2}, {6, 2} },
                        { {8, 1}, {3, 3}, {8, 1}, {3, 3}, {6, 2} } };

    cout << "Level: ";
    cin >> level;
    cout << endl;
    cout << "Row data: " << endl;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 2; k++)
                cin >> rd[i][j][k];
    cout << endl;
    VFCalc calc(level, rd);

    calc.exec();

    // cout << calc.bruteForceWinChances() << endl;

    // cout << calc.possBoards.size() << endl;
    //
    // for (auto x : calc.possBoards) {
    //     cout << "--------" << endl;
    //     for (int j = 0; j < 5; j++) {
    //         for (int i = 0; i < 5; i++) {
    //             cout << (int)(x.tiles[i][j]) << " ";
    //         }
    //         cout << endl;
    //     }
    // }



}
