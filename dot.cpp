#ifndef dot_h
#define dot_h

#include <array>

double dot(const std::array<double,3>& vec1, const std::array<double,3>& vec2) {
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
}

#endif
