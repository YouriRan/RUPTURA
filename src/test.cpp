#include <vector>

#include "test.h"

Test::Test(size_t size_vect):
    vect(static_cast<double>(size_vect))
{}
Test::Test(std::vector<double> _vect) : vect(_vect) {}

std::vector<double> Test::getVector() {
    return vect;
}
