#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "../Shared/populationlib.h"

using namespace std;

#define popSize 3
#define numX 3

/* Population variables */
populationItem* population;
int itemsInPopulation = 0;

/* Test data */
double u0 = 100, 
       u1 = 50,
       u2 = 75,
       u3 = 51;
int x0[] = { 0, 1, 2 },
    x1[] = { 0, 1, 3 },
    x2[] = { 0, 2, 4 },
    x3[] = { 1, 2, 5 };

bool testInsert();
bool testSearch();

int main() {
    cout << "Population tests: " << endl;

    initPopulation(&population, popSize, numX);

    bool testInsertResult = testInsert();
    if (!testInsertResult)
    {
        cout << "testInsert - FALSE" << endl;
        abort();
    }

    cout << "----------" << endl;

    bool testSearchResut = testSearch();
    if (!testSearchResut)
    {
        cout << "testSearch - FALSE" << endl;
        abort();
    }

    cout << endl;
    cout << "Tests OK" << endl;
    system("pause");
}

bool testInsert(){
    insert(population, x0, numX, u0, &itemsInPopulation, popSize);
    if (population[0].locations[0] == x0[0] &&
        population[0].locations[1] == x0[1] &&
        population[0].locations[2] == x0[2] &&
        population[0].solution == u0 &&
        itemsInPopulation == 1 &&
        population[1].solution == -1)
        cout << "testInsert - 1 - OK" << endl;
    else
    {
        cout << "testInsert - 1 - FALSE" << endl;
        return false;
    }

    insert(population, x1, numX, u1, &itemsInPopulation, popSize);
    if (x0[0] == 0 && 
        x0[1] == 1 && 
        x0[2] == 2 &&
        population[0].locations[0] == x0[0] &&
        population[0].locations[1] == x0[1] &&
        population[0].locations[2] == x0[2] &&
        population[0].solution == u0 &&
        population[1].locations[0] == x1[0] &&
        population[1].locations[1] == x1[1] &&
        population[1].locations[2] == x1[2] &&
        population[1].solution == u1 &&
        itemsInPopulation == 2 &&
        population[2].solution == -1)
        cout << "testInsert - 2 - OK" << endl;
    else
    {
        cout << "testInsert - 2 - FALSE" << endl;
        return false;
    }

    insert(population, x2, numX, u2, &itemsInPopulation, popSize);
    if (population[0].locations[0] == x0[0] &&
        population[0].locations[1] == x0[1] &&
        population[0].locations[2] == x0[2] &&
        population[0].solution == u0 &&
        population[1].locations[0] == x1[0] &&
        population[1].locations[1] == x1[1] &&
        population[1].locations[2] == x1[2] &&
        population[1].solution == u1 &&
        population[2].locations[0] == x2[0] &&
        population[2].locations[1] == x2[1] &&
        population[2].locations[2] == x2[2] &&
        population[2].solution == u2 &&
        itemsInPopulation == 3)
        cout << "testInsert - 3 - OK" << endl;
    else
    {
        cout << "testInsert - 3 - FALSE" << endl;
        return false;
    }

    /* Insert - successfull (changes x1 u1) */
    insert(population, x3, numX, u3, &itemsInPopulation, popSize);
    if (population[0].locations[0] == x0[0] &&
        population[0].locations[1] == x0[1] &&
        population[0].locations[2] == x0[2] &&
        population[0].solution == u0 &&
        population[1].locations[0] == x3[0] &&
        population[1].locations[1] == x3[1] &&
        population[1].locations[2] == x3[2] &&
        population[1].solution == u3 &&
        population[2].locations[0] == x2[0] &&
        population[2].locations[1] == x2[1] &&
        population[2].locations[2] == x2[2] &&
        population[2].solution == u2 &&
        itemsInPopulation == 3)
        cout << "testInsert - 4 - OK" << endl;
    else
    {
        cout << "testInsert - 4 - FALSE" << endl;
        return false;
    }

    /* Insert - NOT successfull (insert x1 u1 back) */
    insert(population, x1, numX, u1, &itemsInPopulation, popSize);
    if (population[0].locations[0] == x0[0] &&
        population[0].locations[1] == x0[1] &&
        population[0].locations[2] == x0[2] &&
        population[0].solution == u0 &&
        population[1].locations[0] == x3[0] &&
        population[1].locations[1] == x3[1] &&
        population[1].locations[2] == x3[2] &&
        population[1].solution == u3 &&
        population[2].locations[0] == x2[0] &&
        population[2].locations[1] == x2[1] &&
        population[2].locations[2] == x2[2] &&
        population[2].solution == u2 &&
        itemsInPopulation == 3)
        cout << "testInsert - 5 - OK" << endl;
    else
    {
        cout << "testInsert - 5 - FALSE" << endl;
        return false;
    }

    return true;
}

bool testSearch(){
    populationItem x0Search = search(population, popSize, x0, numX);
    if (x0Search.solution == u0 &&
        x0[0] == 0 && 
        x0[1] == 1 && 
        x0[2] == 2 &&
        x0Search.locations[0] == 0 && 
        x0Search.locations[1] == 1 && 
        x0Search.locations[2] == 2)
        cout << "testSearch - 1 - OK" << endl;
    else
    {
        cout << "testSearch - 1 - FALSE" << endl;
        return false;
    }

    populationItem x1Search = search(population, popSize, x1, numX);
    if (x1Search.solution == -1.0 &&
        x1Search.locations == NULL &&
        x1[0] == 0 &&
        x1[1] == 1 && 
        x1[2] == 3)
        cout << "testSearch - 2 - OK" << endl;
    else
    {
        cout << "testSearch - 2 - FALSE" << endl;
        return false;
    }

    populationItem x2Search = search(population, popSize, x2, numX);
    if (x2Search.solution > -1.0)
        cout << "testSearch - 3 - OK" << endl;
    else
    {
        cout << "testSearch - 3 - FALSE" << endl;
        return false;
    }

    populationItem x3Search = search(population, popSize, x3, numX);
    if (x3Search.solution > -1.0)
        cout << "testSearch - 4 - OK" << endl;
    else
    {
        cout << "testSearch - 4 - FALSE" << endl;
        return false;
    }

    return true;
}