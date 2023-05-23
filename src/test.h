#ifndef TEST_H
#define TEST_H

#include <vector>

struct Test
{
    public:

        Test(size_t size_vect);
        Test(std::vector<double> _vect);

        std::vector<double> getVector();
    
    private:
        std::vector<double> vect;
};

#endif