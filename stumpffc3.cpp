#ifndef stumpff_c_3
#define stumpff_c_3
#include <cmath>

double stumpffc3(double z)
{
    if (z > 0)
    {
        return (sqrt(z) - sin(sqrt(z))) / pow(sqrt(z), 3);
    }
    else if (z == 0)
    {
        return 1.0 / 6.0;
    }
    else
    {
        return (sinh(sqrt(-z)) - sqrt(-z)) / pow(sqrt(-z), 3);
    }
}

#endif

