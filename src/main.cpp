#include "VFCalc.h"
#include "util.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>

void runBoards() {
    unsigned seed = time(0); // 16412511
    cout << seed << endl;
    srand(seed);

    int rd[2][5][2];
    int lvl = 8, profn = 0, prof[3], board[5][5];
    int bcnt = 0;
    while (1) {
        cout << "Level " << lvl << ", profile " << profn << endl;
        prof[0] = boardProfiles[lvl - 1][profn][0];
        prof[1] = boardProfiles[lvl - 1][profn][1];
        prof[2] = boardProfiles[lvl - 1][profn][2];
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
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 5; j++)
                cout << rd[i][j][0] << " " << rd[i][j][1] << " ";
        cout << endl;

        VFCalc calc(lvl, rd);
        bcnt = calc.getNumBoards();

        vfBoardI vt;

        bool d = calc.doable(vt);
        d = 1;
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

int main() {

    // testHash();
    // return 0;
    // runBoards();
    // return 0;

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
