#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <sys/time.h>
#include "../Shared/RDOAlib.h"

/*  
    GenSolution
    0 - RDOA, 1 - RDOA-D
    EvalSolution
    0 - BinaryRule, 1 - PartialyBinaryRule 
*/
#define GEN_SOLUTION 0
#define EVAL_SOLUTION 0
#define ITERS 10

using namespace std;

/* Configuration */
int numDP = 100;      // Vietoviu skaicius (demand points, max 10000)
int numPF = 5;          // Esanciu objektu skaicius (preexisting facilities)
int numF  = 3;          // Esanciu imoniu skaicius (firms)
int numCL = 25;         // Kandidatu naujiems objektams skaicius (candidate locations)
int numX  = 3;          // Nauju objektu skaicius

double **demandPoints, **distances;
int *X, *bestX, *ranks;

void updateRanks(bool success);

int main() {
    double ts_start = getTime();

    //New seed on every run
    srand((unsigned)time(0));
    
	loadDemandPoints(numDP, &demandPoints);
	calculateDistances(numDP, &distances, demandPoints);
	
    X = new int[numX];
	bestX = new int[numX];
	double bestU = -1;
    
    //Init ranks
    ranks = new int[numCL];
    for (int i = 0; i < numCL; i++)
        ranks[i] = 1;

    randomSolution(numCL, numX, X);
    double u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, 2, EVAL_SOLUTION);
    bestU = u;
    for (int i = 0; i < numX; i++) 
        bestX[i] = X[i];
	
    for (int iters = 0; iters < ITERS; iters++) {
        printf("iteration - %d \n", iters);
        generateSolution(numX, numCL, X, bestX, ranks, distances, GEN_SOLUTION);
        u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, 2, EVAL_SOLUTION);

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

	resultsFile << ", " << bestU << ", " << getTime() - ts_start << endl;
}
