#ifndef UTIL_H
#define UTIL_H

#include <limits>

using namespace std;

// https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
// iterates over combinations of s bits with x of them set to 1.
#define ITERPERM(s, x, i, t) for(i = (1 << x) - 1; i < 1 << s; t = (i | (i - 1)) + 1, i = (i == 0 ? 1 << s : t | ( ((t & -t)/(i & -i) >> 1) - 1)))

inline constexpr double dbl_inf = numeric_limits<double>::infinity();

inline constexpr size_t pow5table[5][5] = {
    {1, 5, 25, 125, 625},
    {3125, 15625, 78125, 390625, 1953125},
    {9765625, 48828125, 244140625, 1220703125, 6103515625},
    {30517578125, 152587890625, 762939453125, 3814697265625, 19073486328125},
    {95367431640625, 476837158203125, 2384185791015625, 11920928955078125, 59604644775390625}
};

inline constexpr int boardProfiles[8][5][3] = {
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

int cnt(int x);
bool dblEq(double a, double b);

#endif
