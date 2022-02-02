#include "util.h"
#include <cmath>

// https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetKernighan
int cnt(int x) {
    int c;
    for (c = 0; x; c++)
        x &= x - 1;
    return c;
}

bool dblEq(double a, double b) {
    return fabs(a - b) < 0.000001;
}
