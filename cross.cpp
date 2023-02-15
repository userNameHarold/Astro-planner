#ifndef cross_h
#define cross_h

#include <vector>

std::vector<double> cross(const std::vector<double> *v1, const std::vector<double> *v2) {
    std::vector<double> result;
    result[0] = (*v1)[1] * (*v2)[2] - (*v1)[2] * (*v2)[1];
    result[1] = (*v1)[2] * (*v2)[0] - (*v1)[0] * (*v2)[2];
    result[2] = (*v1)[0] * (*v2)[1] - (*v1)[1] * (*v2)[0];
    return result;
}

#endif
