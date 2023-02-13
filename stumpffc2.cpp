#ifndef stumpff_c_2
#define stumpff_c_2
#include <cmath>

double stumpffc2(double z) {
    if (z > 0) {
        return (1 - cos(sqrt(z))) / z;
    } else if (z == 0) {
        return 0.5;
    } else {
        return (cosh(sqrt(-1 * z)) - 1) / (-1 * z);
    }
}

#endif

