#include <cstdlib>

// Functions for allocating/deallocating memory

double * create_vector (size_t a)
{
    double * v = new double[a];
    return v;
}

void free_vector (double * v)
{
    delete [] v;
}
