#ifndef dot_h
#define dot_h

#include <vector>

/* void dot(const double *matrix1, const double *matrix2, double *result, int rows1, int cols1, int cols2) {
    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols2; ++j) {
            double dotProduct = 0.0;
            for (int k = 0; k < cols1; ++k) {
                dotProduct += matrix1[i*cols1 + k] * matrix2[k*cols2 + j];
            }
            result[i*cols2 + j] = dotProduct;
        }
    }
} */


double dot(const std::vector<double> *vec1, const std::vector<double> *vec2) {
    return (*vec1)[0] * (*vec2)[0] + (*vec1)[1] * (*vec2)[1] + (*vec1)[2] * (*vec2)[2];
} 

#endif
