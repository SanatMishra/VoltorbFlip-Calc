#include "VFCalc.h"
#include "util.h"
#include <limits>
#include <iostream>
#include <cmath>

vfBoard::vfBoard() {
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            tiles[i][j] = 4;
}

void vfBoardI::manualHash() {
    hash = 0;
    for (int i = 4; i >= 0; i--)
        for (int j = 4; j >= 0; j--)
            hash = 5*hash + (tiles[i][j] & 3);
}

// Hashing all 25 tiles becomes expensive, but we can recalculate the hash
// with a few operations if only one tile is changed (as is usually the case).
// We also ignore the case of the tile being 0, as it is never set to 0 during the calculation.
void vfBoardI::setTileAndHash(int x, int y, int v) {
    hash += pow5table[x][y]*((v & 3) - (tiles[x][y] & 3));
    tiles[x][y] = v;
}

bool vfBoardI::operator== (const vfBoardI& other) const {
    return hash == other.hash;
}

bool vfBoardI::operator< (const vfBoardI& other) const {
    return hash < other.hash;
}

size_t std::hash<vfBoardI>::operator()(const vfBoardI &v) const {
    return v.hash;
}

void VFCalc::getSubBoards(vfBoardI& vt, vector<int>& subBoards) {
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

double VFCalc::bruteForceWinChances() {
        vfBoardI vt;
        return bruteForceWinChances(vt);
}

double VFCalc::bruteForceWinChances(vfBoardI& vt) {
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

double VFCalc::bruteForceWinChances(vfBoardI& curPos, vector<uint16_t>& sub, int plusData[2][5], double tw) {
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
                        plusData[0][y] += (curPos.tiles[x][y] - 1) % 3 - i % 3;
                        plusData[1][x] += (curPos.tiles[x][y] - 1) % 3 - i % 3;
                        curPos.setTileAndHash(x, y, i + 1);
                        ret2 += w[i]*bruteForceWinChances(curPos, subn[i], plusData, w[i]);
                        // T[curPos] should always exist
                        retd2 += w[i]*T[curPos].gameTime;
                    }
                }
                ret2 /= tw;
                retd2 = 1 + retd2/tw;

                plusData[0][y] += (curPos.tiles[x][y] - 1) % 3;
                plusData[1][x] += (curPos.tiles[x][y] - 1) % 3;
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

VFCalc::VFCalc(int level, int rd[2][5][2]) {
    this->level = level;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 5; j++) {
            rowData[i][j][0] = rd[i][j][0];
            rowData[i][j][1] = rd[i][j][1];
            rowData[i][j][2] = rd[i][j][1] + rd[i][j][0] - 5;
        }
    }

    getPossProfiles();
    getAllBoards();
    cout << possBoards.size() << " possible boards" << endl;
}

void VFCalc::exec() {
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

bool VFCalc::doable(vfBoardI& vt) {
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

    return doabilityEst(d, I, T) < doabilityEst(8.5, 45, 20); // pow(10, 29);
    // return (tgamma(d + 1)*pow(I/d, d) < 542097412.116);
}

double VFCalc::doabilityEst(double d, double I, double T) {
    return d*pow(I, (I/T - 1)*(d + 1));
}

void VFCalc::getInDepthData(vfBoardI& vt, double wc[5][5], double gt[5][5]) {
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

// Our goal is to find a combination of tiles that lowers the board plusness to 0 with the highest win chance possible.
// We can treat each tile as an independent choice that has p_v chance to fail and makes (2*p_3 + p_2)/(1 - p_v) progress
// towards 0 plusness on success. We need to maximize success chance, i.e. product of (1 - p_v), or sum of ln(1 - p_v),
// for all valid combinations of tiles picked. Estimating that the tiles are independent, we treat this like a knapsack
// problem and just take the highest ratio of ln(1 - p_v)/( (2*p_3 + p_2)/(1 - p_v) ) = (1 - p_v)ln(1 - p_v)/(2*p_3 + p_2)
void VFCalc::getHeuristics(vfBoardI& vt, double wh[5][5], double vr[5][5]) {
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

void VFCalc::display(vfBoardI& vt, bool useInDepthData) {
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
            if (useInDepthData) {
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

    if (useInDepthData) {
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

int VFCalc::getNumBoards() {
        return possBoards.size();
    }
