#include "util.h"
#include "VFCalc.h"

void VFCalc::getPossProfiles() {
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
}

void VFCalc::getAllBoards() {
    int pbIndex = 0; // tracks size of possBoards before next loop iteration begins, so we know which boards to assign weights to

    for(int p : possProfiles) {
        // cout << p << endl;
        vfBoardC cur;
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                cur.tiles[i][j] = 1;

        int cnt3 = boardProfiles[level - 1][p][1];

        int rd2[2][5], v2[2][5];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 5; j++) {
                rd2[i][j] = rowData[i][j][2];
                v2[i][j] = rowData[i][j][1];
            }
        }

        int k[5]; //3s
        int m[5];
        vector<int> nzi[5];

        int r[5]; // 2s
        vector<int> oi[5];

        int s[5]; // voltorbs
        vector<int> pi[5];

        int tp; // used for permutations

        int variants = 0; // number of boards for this profile, used for weights at the end

        // insert tiles in the order 3, 2, v, from most to least restrictive.
        // Probably unnecessarily complex but this used to have significant runtime
        // For 3s we must account for boards that do not fill up all potential spaces with 3s
        // There might be a boardProfile with 4 3s, even though the board itself can fit 5.
        // k stores the amount of 3s in each row; we iterate over possible arrangements of k over the rows

        for (k[0] = 0; k[0] <= min(rd2[0][0]/2, cnt3); k[0]++) {
        for (k[1] = 0; k[1] <= min(rd2[0][1]/2, cnt3 - k[0]); k[1]++) {
        for (k[2] = 0; k[2] <= min(rd2[0][2]/2, cnt3 - k[0] - k[1]); k[2]++) {
        for (k[3] = 0; k[3] <= min(rd2[0][3]/2, cnt3 - k[0] - k[1] - k[2]); k[3]++) {
            k[4] = cnt3 - k[0] - k[1] - k[2] - k[3];
            if (k[4] > rd2[0][4]/2)
                continue;
            // cout << k[0] << k[1] << k[2] << k[3] << k[4] << endl;
        // ITERPERM(sp, cnt3, i, tp) {
        //     #define MASK(a, b) ( ((1 << b) - 1) ^ ((1 << a) - 1) )
        //     for (int j = 0; j < 5; j++) {
        //         k[j] = cnt(i & MASK(ind[j], ind[j + 1]));
        //     }

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

                                                                        // cout << endl;
                                                                        // for (int x = 0; x < 5; x++) {
                                                                        //     for (int y = 0; y < 5; y++) {
                                                                        //         cout << (int)(cur.tiles[x][y]) << " ";
                                                                        //     }
                                                                        //     cout << endl;
                                                                        // }
                                                                        // cout << endl;

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
        }}}}

        for (int ind = pbIndex; ind < pbIndex + variants; ind++)
            possBoards[ind].weight = 1.0/variants;
        pbIndex += variants;

        // cout << p << ": " << variants << endl;

    }
}
