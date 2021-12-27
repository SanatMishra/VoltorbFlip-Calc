#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <unordered_map>
#include <vector>

// https://graphics.stanford.edu/~seander/bithacks.html
// iterates over combinations of s bits with x of them set to 1. May be target for speedup
#define ITERPERM(s, x, i, t) for(i = (1 << x) - 1; i < 1 << s; t = (i | (i - 1)) + 1, i = (i == 0 ? 1 << s : t | ( ((t & -t)/(i & -i) >> 1) - 1)))
#define MASK(a, b) ( ((1 << b) - 1) ^ ((1 << a) - 1) )

using namespace std;

// https://graphics.stanford.edu/~seander/bithacks.html
int cnt(int x) {
    x = (x & 0x55555555) + (x >> 1 & 0x55555555);
    x = (x & 0x33333333) + (x >> 2 & 0x33333333);
    return ((x + (x >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
}

bool dblEq(double a, double b) {
    return fabs(a - b) < 0.000001;
}

// board with potential unsolved tiles, does not include row data (should it?)
// does not include plusness (might need to, new class?)
struct vfBoard {
    uint8_t tiles[5][5];

    bool operator== (const vfBoard& other) const {
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                if (tiles[i][j] != other.tiles[i][j])
                    return 0;
            }
        }
        return 1;
    }
};

template<>
struct std::hash<vfBoard> {
    size_t operator()(const vfBoard &v) const {
        size_t ret = 0;
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                ret = 5*ret + v.tiles[i][j];
        return ret;
    }
};

struct posData {
    int bmX, bmY;
    double winChance;
    double gameTime;
    int possBoards;
};

class VFCalc {

public:
    int level;
    int rowData[2][5][3];       // first axis: [0]: row, [1]: col
                                // 3: score, voltorbs, plusness
                                // Consider "plusness" of a row to be (number of 2s) + 2*(number of 3s) for any valid board
                                // A 4/1 row is +0, 5/1 is +1, 5/2 is +2, 6/3 is +4, 10/1 is +6, etc.

    static int boardProfiles[8][5][3];

    vector<vfBoard> possBoards;
    unordered_map<vfBoard, posData> T;

// public:
    VFCalc() {
        // todo: make plusness data inaccessible, auto-updated in set tile func
        level = 1;
    }

    void newBoard(int level, int rd[2][5][2]) {
        // todo: make plusness data inaccessible, auto-updated in setTile func
        // move some shit 2 constructor when you know what 2 do
        this->level = level;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 5; j++) {
                rowData[i][j][0] = rd[i][j][0];
                rowData[i][j][1] = rd[i][j][1];
                rowData[i][j][2] = rd[i][j][1] + rd[i][j][0] - 5;
            }
        }
        possBoards.clear();
        T.clear();
        getPossibleBoards();
        cout << possBoards.size() << " possible boards" << endl;
    }

    // not necessary for now
    bool solved(vfBoard& board) {
        int h = 0;
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                if (board.tiles[i][j] == 2 || board.tiles[i][j] == 3)
                    h += board.tiles[i][j] - 1;
            }
        }
        return h == 0;
        // return h == totalPlus;
    }

    void getPossibleBoards() {
        possBoards.clear();

        bool poss[5] = {0};
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
                poss[i] = 1;
            }
        }
        for(int p = 0; p < 5; p++) {
            if (!poss[p])
                continue;

            vfBoard cur;
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++)
                    cur.tiles[i][j] = 1;

            // Some of this function was written with 3, 2, 1 notated as 2, 1, 0 (contribution to row plusness). Apologies in advance
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

            // insert tiles in the order 3, 2, v, from most to least restrictive.
            // Probably unnecessarily complex but I thought this would have significant runtime
            // For 3s we must account for boards that do not fill up all potential spaces with 3s
            // There might be a boardProfile with 4 3s, even though the board itself can fit 5.
            // We use k instead of the rd2, and iterate over possible arrangements of k over the rows
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

        }
    }

    void exec() {
        cout << "k" << endl;
        winChance();
        cout << "k" << endl;

        vfBoard vt;
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                vt.tiles[i][j] = 4;
            }
        }

        int x, y, v;
        double wc, vr, gt;
        while(1) {
            cout << endl;
            for (int j = 0; j < 5; j++) {
                for (int i = 0; i < 5; i++) {
                    cout << "  " << (int)(vt.tiles[i][j]) << "   ";
                }
                cout << "| " << rowData[0][j][0] << " " << rowData[0][j][1] << " ";

                for (int i = 0; i < 5; i++) {
                    winDataFromPosition(vt, i, j, wc, vr, gt);
                    if (wc != -1) {
                        printf("%.3f %.3f %.3f|", wc, vr, gt);
                    } else {
                        cout << "        ?        |";
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

            cout << T[vt].bmX << " " << T[vt].bmY << " " << T[vt].winChance << endl;

            cin >> x >> y >> v;
            // x = T[vt].bmX;
            // y = T[vt].bmY;
            // cin >> v;

            vt.tiles[x][y] = v;
            if (T.find(vt) == T.end()) {
                cout << "-- NEW POSITION --" << endl;
                winChance(vt, 0);
            }
        }
    }

    // debug info
    int wc_ex, wc_te;

    double winChance() {
        vfBoard vt;
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                vt.tiles[i][j] = 4;
            }
        }

        return winChance(vt, 0);
    }

    double winChance(vfBoard& vt, bool debug) {
        vector<int> sub;
        for (int i = 0; i < possBoards.size(); i++)
            sub.push_back(i);

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

        wc_ex = 0;
        wc_te = 0;
        double ret = winChance(vt, sub, plusData, debug);
        cout << "E " << wc_ex << endl;
        cout << "T " << wc_te << endl;
        return ret;
    }

    // todo: take 0 voltorb rows immediately
    //       get win chance for all tiles
    double winChance(vfBoard& curPos, vector<int>& sub, int plusData[2][5], bool debug) {
        if (debug) {
            cout << endl << "wc " << sub.size() << endl;
            for (int j = 0; j < 5; j++) {
                for (int i = 0; i < 5; i++) {
                    cout << (int)(curPos.tiles[i][j]) << " ";
                }
                cout << "| " << plusData[0][j];
                cout << endl;
            }

            cout << "----------+" << endl;
            for (int i = 0; i < 5; i++)
                cout << plusData[1][i] << " ";
            cout << endl << endl;

        }
        wc_ex++;

        if (T.find(curPos) != T.end()) return T[curPos].winChance;
        // do we need this?
        // gprof later, maybe initial solved check thats faster than the loop
        // if (solved(curPos)) {
        //     T[curPos] = {5, 5, 1};
        //     return 1;
        // }

        double ret = 0, retd = 25;
        vfBoard workingPos = curPos;
        int bmX = 5, bmY = 5;
        for (int x = 0; x < 5; x++) {
            for (int y = 0; y < 5; y++) {
                if (curPos.tiles[x][y] != 4 || plusData[0][y] == 0 || plusData[1][x] == 0)
                    continue;

                double ret2 = 0, retd2 = 0;
                vector<int> subn[3];

                // skip count if curPos [xy]=i is solved --> get count from child calls?
                // does that require knowing whether solved or not before entering loop
                for (int i = 0; i < sub.size(); i++) {
                    if (possBoards[sub[i]].tiles[x][y] > 0) {
                        subn[possBoards[sub[i]].tiles[x][y] - 1].push_back(sub[i]);
                    }
                }
                for (int i = 0; i < 3; i++) {
                    if (subn[i].size() > 0) {
                        // maybe conglomerate this @ some point
                        plusData[0][y] -= (i + 1) - (curPos.tiles[x][y] == 4 ? 1 : curPos.tiles[x][y]);
                        plusData[1][x] -= (i + 1) - (curPos.tiles[x][y] == 4 ? 1 : curPos.tiles[x][y]);
                        curPos.tiles[x][y] = i + 1;
                        ret2 += subn[i].size()*winChance(curPos, subn[i], plusData, debug);
                        // T[curPos] should always exist
                        retd2 += subn[i].size()*T[curPos].gameTime;
                    }
                }
                ret2 /= sub.size();
                retd2 = 1 + retd2/sub.size();

                plusData[0][y] += (curPos.tiles[x][y] == 4 ? 1 : curPos.tiles[x][y]) - 1;
                plusData[1][x] += (curPos.tiles[x][y] == 4 ? 1 : curPos.tiles[x][y]) - 1;
                curPos.tiles[x][y] = 4;
                if (ret2 > ret || dblEq(ret2, ret) && retd2 < retd) {
                    bmX = x;
                    bmY = y;
                    ret = ret2;
                    retd = retd2;
                }
            }
        }

        wc_te++;
        T[curPos] = {bmX, bmY, bmX == 5 ? 1 : ret, bmX == 5 ? 0 : retd, (int)(sub.size())};
        return bmX == 5 ? 1 : ret;
    }

    // returns -1 if position not in table, (hopefully) meaning dead tile
void winDataFromPosition(vfBoard& vt, int x, int y, double& wc, double& vr, double& gt) {
        if (vt.tiles[x][y] != 4 || T.find(vt) == T.end()) {
            wc = vr = gt = -1;
            return;
        }

        wc = vr = gt = 0;
        int temp = vt.tiles[x][y];
        for (int i = 1; i <= 3; i++) {
            vt.tiles[x][y] = i;
            if (T.find(vt) != T.end()) {
                wc += T[vt].possBoards*T[vt].winChance;
                vr += T[vt].possBoards;
                gt += T[vt].possBoards*T[vt].gameTime;
            }
        }
        vt.tiles[x][y] = temp;

        if (wc == 0) {
            wc = -1; vr = -1; gt = -1;
            return;
        }

        wc /= T[vt].possBoards;
        vr = 1 - vr/T[vt].possBoards;
        gt = 1 + gt/T[vt].possBoards;

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

int main() {
    srand(time(0));

    VFCalc calc;
    int rd[2][5][2];
    int lvl = 6, profn = 0, prof[3], board[5][5];
    int bcnt = 0;
    while (1) {
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
                rd[0][j][0] += board[j][k];
                rd[1][k][0] += board[j][k];
                rd[0][j][1] += (board[j][k] == 0);
                rd[1][k][1] += (board[j][k] == 0);
            }
        }

        calc.newBoard(lvl, rd);
        bcnt = calc.possBoards.size();
        calc.winChance();

        if (bcnt >= 5000) break;
        profn++;
        if (profn >= 5) {
            profn = 0;
            lvl++;
            if (lvl >= 9)
                lvl = 6;
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

    // VFCalc calc;
    // int level = 8;
    // int rd[2][5][2] = { { {7, 2}, {6, 1}, {3, 3}, {6, 2}, {6, 2} },
    //                     { {8, 1}, {3, 3}, {8, 1}, {3, 3}, {6, 2} } };
    //
    //
    //
    // // cout << "Level: ";
    // // cin >> level;
    // // cout << endl;
    // // cout << "Row data: " << endl;
    // // for (int i = 0; i < 2; i++)
    // //     for (int j = 0; j < 5; j++)
    // //         for (int k = 0; k < 2; k++)
    // //             cin >> rd[i][j][k];
    // // cout << endl;
    // calc.newBoard(level, rd);
    //
    // calc.exec();
    //
    // // cout << calc.winChance() << endl;
    //
    // // cout << calc.possBoards.size() << endl;
    // //
    // // for (auto x : calc.possBoards) {
    // //     cout << "--------" << endl;
    // //     for (int j = 0; j < 5; j++) {
    // //         for (int i = 0; i < 5; i++) {
    // //             cout << (int)(x.tiles[i][j]) << " ";
    // //         }
    // //         cout << endl;
    // //     }
    // // }
    //
    // // while (true) {
    // //     (menu)
    // //     for (int lvl = 1; lvl <= 8; lvl++) {
    // //         input for lvl and rd
    // //         calc.newBoard(lvl, rd);
    // //         calc.display();
    // //         calc.getWinChanceData(); // make sure update tile display
    // //
    // //         int x, y, v;
    // //         cin >> x >> y >> v;
    // //         if (v == 0) break;
    // //
    // //         calc.enterTile(x, y, v);
    // //     }
    // // }
    // // // ...


}
