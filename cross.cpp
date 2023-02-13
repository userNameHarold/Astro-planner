#ifndef cross_h
#define cross_h

#include <array>

std::array<double, 3> cross(const std::array<double, 3> &v1, const std::array<double, 3> &v2) {
    std::array<double, 3> result;
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return result;
}

#endif
