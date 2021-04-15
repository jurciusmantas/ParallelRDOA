#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sstream>
#include "../Shared/RDOAlib.h"
#include "../Shared/populationLib.h"

/*  
    GenSolution
    0 - RDOA, 1 - RDOA-D
    EvalSolution
    0 - BinaryRule, 1 - PartialyBinaryRule 
*/
#define GEN_SOLUTION 0
#define EVAL_SOLUTION 0
#define ITERS 100
#define POP_SIZE 50

using namespace std;

/* Configuration */
int numDP   = 100;      // Vietoviu skaicius (demand points, max 10000)
int numPF   = 5;          // Esanciu objektu skaicius (preexisting facilities)
int numF    = 3;          // Esanciu imoniu skaicius (firms)
int numCL   = 25;         // Kandidatu naujiems objektams skaicius (candidate locations)
int numX    = 3;          // Nauju objektu skaicius

double **demandPoints, **distances;
int *X, *bestX, *ranks;

// Population variables
populationItem* population;
int itemsInPopulation = 0;
int timesPopulationSaved = 0;

int main() {
    double ts_start = getTime();

    //New seed on every run
    srand((unsigned)time(0));

	loadDemandPoints(numDP, &demandPoints);
	calculateDistances(numDP, &distances, demandPoints);

    initPopulation(&population, POP_SIZE, numX);
	
    X = new int[numX];
	bestX = new int[numX];
	double bestU = -1;
    
    //Init ranks
    ranks = new int[numCL];
    for (int i = 0; i < numCL; i++)
        ranks[i] = 1;

    randomSolution(numCL, numX, X);
    double u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, EVAL_SOLUTION);
    bestU = u;
    for (int i = 0; i < numX; i++) 
        bestX[i] = X[i];
    insert(population, X, numX, u, &itemsInPopulation, POP_SIZE);
	
    for (int iters = 0; iters < ITERS; iters++) {
        printf("iteration - %d \n", iters);
        
        generateSolution(numX, numCL, X, bestX, ranks, distances, GEN_SOLUTION);

        //Search for solution in population
        populationItem popItem = search(population, POP_SIZE, X, numX);
        if (popItem.solution > -1.0)
        {
            /* Generated solution was found in population */
            u = popItem.solution;
            timesPopulationSaved++;
        }
        else
        {
            /* Generated solution was not found in population */
            u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, EVAL_SOLUTION);
            insert(population, X, numX, u, &itemsInPopulation, POP_SIZE);
        }

        if (u > bestU) 
        {
            updateRanks(ranks, X, bestX, numCL, numX, true);
            bestU = u;

            for (int i = 0; i < numX; i++) 
                bestX[i] = X[i];
        }
        else
            updateRanks(ranks, X, bestX, numCL, numX, false);
    }

    // Write results
    ofstream resultsFile;
    stringstream fileName;
    fileName << "results" << GEN_SOLUTION << EVAL_SOLUTION << ".txt" << endl;
    resultsFile.open(fileName.str(), ios_base::app);
	for (int i=0; i<numX; i++) 
        resultsFile << bestX[i] << " ";
    
	resultsFile << ", " << timesPopulationSaved << ", " << bestU << ", " << getTime() - ts_start << endl;
}
