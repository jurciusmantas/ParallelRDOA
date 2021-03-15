#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <float.h>
#include "populationLib.h"

#pragma region Helpers

// O(n) speed for checking if two arrays contain the same items
bool containsSameItems(int* arr1, int* arr2, int size)
{
    int xor1 = arr1[0];
    int xor2 = arr2[0];

    for (int i = 1; i < size; i++)
    {
        xor1 ^= arr1[i];
        xor2 ^= arr2[i];
    }

    int xorBoth = xor1 ^ xor2;

    return xorBoth == 0;
}

#pragma endregion Helpers

void initPopulation(populationItem** population, int populationSize, int numX)
{
    std::cout << "initPopulation START" << std::endl;

    (*population) = new populationItem[populationSize];
    for (int i = 0; i < populationSize; i++)
    {
        (*population)[i].solution = -1;
        (*population)[i].locations = new int[numX];
    }

    std::cout << "initPopulation END" << std::endl;
}

populationItem search(populationItem* population, int populationSize, int* X, int numX)
{
    populationItem result = { -1.0 /* solution value of emtpy struct */ };
    for (int i = 0; i < populationSize; i++)
    {
        if (population[i].solution == -1)
            break;

        if (containsSameItems(X, population[i].locations, numX))
        {
            result = population[i];
            break;
        }
    }

    return result;
}

void insert(populationItem* population, int* locations, int numX, double solution, int* itemsInPopulation, int populationSize)
{
    if ((*itemsInPopulation) < populationSize)
    {
        for (int j = 0; j < numX; j++)
            population[(*itemsInPopulation)].locations[j] = locations[j];
        population[(*itemsInPopulation)].solution = solution;

        (*itemsInPopulation)++;
    }

    else
    {
        double worst = DBL_MAX;
        int worstIndex = -1;
        for (int i = 0; i < populationSize; i++)
        {
            if (worst > population[i].solution)
            {
                worst = population[i].solution;
                worstIndex = i;
            }
        }

        if (worstIndex > 0 && solution > worst)
        {
            for (int j = 0; j < numX; j++)
                population[worstIndex].locations[j] = locations[j];
            population[worstIndex].solution = solution;
        }
    }
}
