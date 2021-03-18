#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "../Shared/parallelRDOAlib.h"

/*  
    GenSolution
    0 - RDOA, 1 - RDOA-D
    EvalSolution
    0 - BinaryRule, 1 - PartialyBinaryRule 
*/
#define GEN_SOLUTION 0
#define EVAL_SOLUTION 0
#define ITERS 10000

using namespace std;

/* Configuration */
int numDP = 10000;      // Vietoviu skaicius (demand points, max 10000)
int numPF = 5;          // Esanciu objektu skaicius (preexisting facilities)
int numF  = 3;          // Esanciu imoniu skaicius (firms)
int numCL = 25;         // Kandidatu naujiems objektams skaicius (candidate locations)
int numX  = 3;          // Nauju objektu skaicius

double **demandPoints, **distances;
int *X, *bestX, *ranks;

void updateRanks(int success);

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
        ranks[i] = i + 1;

    randomSolution(numCL, numX, X);
    double u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, EVAL_SOLUTION);
    bestU = u;
    for (int i = 0; i < numX; i++) 
        bestX[i] = X[i];
	
    for (int iters = 0; iters < ITERS; iters++) {
        printf("iteration - %d \n", iters);
        generateSolution(numX, numCL, X, bestX, ranks, distances, GEN_SOLUTION);
        u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, EVAL_SOLUTION);

        if (u > bestU) 
        {
            updateRanks(1);
            bestU = u;

            /*  Note:
                Shouldn't we generate the solution BEFORE assigning X = X'?
                If we do it affter assigning, the locations that were in X but weren't in X'
                become available to pick.
            */
            for (int i = 0; i < numX; i++) 
                bestX[i] = X[i];
        }
        else
            updateRanks(0);
    }

    // Write results
    ofstream resultsFile;
    resultsFile.open("results.txt", ios_base::app);
	for (int i=0; i<numX; i++) 
        resultsFile << bestX[i] << " ";

	resultsFile << "(" << bestU << "), " << getTime() - ts_start << endl;
}

void updateRanks(int success)
{
    int zeroExists = 0;

    for (int i = 0; i < numX; i++)
    {
        if (success == 1)
        {
            ranks[X[i]]++;

            //Search if X contains bestX[i]
            int contains = 0;
            for (int j = 0; j < numX; j++)
            {
                if (X[j] == bestX[i])
                {
                    contains = 1;
                    break;
                }
            }

            if (contains == 0)
            {
                ranks[bestX[i]]--;
                if (ranks[bestX[i]] == 0)
                    zeroExists = 1;
            }
        }

        if (success == 0)
        {
            //Search if bestX contains X[i]
            int contains = 0;
            for (int j = 0; j < numX; j++)
            {
                if (X[i] == bestX[j])
                {
                    contains = 1;
                    break;
                }
            }

            if (contains == 0)
            {
                ranks[X[i]]--;
                if (ranks[X[i]] == 0)
                    zeroExists = 1;
            }
        }
    }

    if (zeroExists)
        for(int i = 0; i < numCL; i++)
            ranks[i]++;
}