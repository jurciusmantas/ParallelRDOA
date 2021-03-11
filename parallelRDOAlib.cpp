#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "parallelRDOAlib.h"

using std::cout;
using std::endl;

int rouletteWheel(double* probabilities, int count)
{
    double rndNumber = rand() / (double) RAND_MAX;
    double offset = 0.0;
    int pick = 0;

    for (int i = 0; i < count; i++)
    {
        offset += probabilities[i];
        if (rndNumber < offset)
        {
            pick = i;
            break;
        }
    }

    return pick;
}