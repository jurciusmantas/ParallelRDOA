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

void insert(populationItem* population, int* locations, int numX, double solution, int* itemsInPopulation, int populationSize, int* popRanks = NULL)
{
    if ((*itemsInPopulation) < populationSize)
    {
        for (int j = 0; j < numX; j++)
            population[(*itemsInPopulation)].locations[j] = locations[j];
        population[(*itemsInPopulation)].solution = solution;

        if (popRanks != NULL)
            updatePopulationRanks(popRanks, numX, population[(*itemsInPopulation)], { -1 });

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
            if (popRanks != NULL)
                updatePopulationRanks(popRanks, numX, { solution, locations }, population[worstIndex]);
            
            for (int j = 0; j < numX; j++)
                population[worstIndex].locations[j] = locations[j];
            population[worstIndex].solution = solution;
        }
    }
}

void updatePopulationRanks(int* ranks, int numX, populationItem inserted, populationItem removed)
{
    for (int i = 0; i < numX; i++)
    {
        if (inserted.solution != -1 && inserted.locations != NULL)
            ranks[inserted.locations[i]]++;

        if (removed.solution != -1 && removed.locations != NULL)
            ranks[removed.locations[i]]--;
    }
}
