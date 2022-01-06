#include "util.h"
#include <cmath>

// https://graphics.stanford.edu/~seander/bithacks.html PUT SECTION HERE
int cnt(int x) {
    x = (x & 0x55555555) + (x >> 1 & 0x55555555);
    x = (x & 0x33333333) + (x >> 2 & 0x33333333);
    return ((x + (x >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
}

bool dblEq(double a, double b) {
    return fabs(a - b) < 0.000001;
}
