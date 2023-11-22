#include <iostream>
#include "domain.hpp"

// to set, for the discretization of the domain, the N: -n number_of_elements

using namespace std;

int main(int argc, char** argv)
{
    DOMAIN domain(argc, argv);
    domain.printN();
    

    return 0;
}