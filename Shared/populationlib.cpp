#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <float.h>
#include "populationlib.h"
#include "RDOAlib.h"

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

        if (worstIndex > -1 && solution > worst)
        {
            for (int j = 0; j < numX; j++)
                population[worstIndex].locations[j] = locations[j];
            population[worstIndex].solution = solution;
        }
    }
}
