#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>
#include "../Shared/RDOAlib.h"
#include "../Shared/populationlib.h"
#include "../Shared/parallellib.h"

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

double **demandPoints, *distances;
int *X, *bestX, *ranks;

// Population variables
populationItem* population;
int itemsInPopulation = 0;
int timesPopulationSaved = 0;

// Parallel variables
int id, numProcs, offset, procChunkSize;
MPI_Status *statuses;
MPI_Request *requests;
MPI_Status status;

void updateRanks(bool success);

int main(int argc , char * argv []) {
    MPI_Init(&argc , &argv);
    MPI_Comm_rank(MPI_COMM_WORLD ,&id);
    MPI_Comm_size(MPI_COMM_WORLD ,&numProcs);
    
    double ts_start = getTime();

    //New seed on every run
    srand((unsigned)time(0));

	loadDemandPoints(numDP, &demandPoints);
	calculateDistancesAsync(numDP, numProcs, id, &distances, demandPoints);

    if (id == 0)
    {
        ofstream distancesFile;
        distancesFile.open("distances.txt", ios_base::app);
        for (int i = 0; i < numDP; i ++){
            for (int j = 0; j < numDP; j++)
                distancesFile << distances[numDP * i + j] << ", ";
            distancesFile << endl;
        }
    }

    initPopulation(&population, POP_SIZE, numX);
	
    // X = new int[numX];
	// bestX = new int[numX];
	// double bestU = -1;
    
    // //Init ranks
    // ranks = new int[numCL];
    // for (int i = 0; i < numCL; i++)
    //     ranks[i] = i + 1;

    // randomSolution(numCL, numX, X);
    // double u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, EVAL_SOLUTION);
    // bestU = u;
    // for (int i = 0; i < numX; i++) 
    //     bestX[i] = X[i];
	
    // for (int iters = 0; iters < ITERS; iters++) {
    //     printf("iteration - %d \n", iters);
        
    //     generateSolution(numX, numCL, X, bestX, ranks, distances, GEN_SOLUTION);

    //     //Search for solution in population
    //     populationItem popItem = search(population, POP_SIZE, X, numX);
    //     if (popItem.solution > -1.0)
    //     {
    //         /* Generated solution was found in population */
    //         u = popItem.solution;
    //         timesPopulationSaved++;
    //     }
    //     else
    //     {
    //         /* Generated solution was not found in population */
    //         u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, EVAL_SOLUTION);
    //         insert(population, X, numX, u, &itemsInPopulation, POP_SIZE);
    //     }

    //     if (u > bestU) 
    //     {
    //         updateRanks(true);
    //         bestU = u;

    //         /*  Note:
    //             Shouldn't we generate the solution BEFORE assigning X = X'?
    //             If we do it affter assigning, the locations that were in X but weren't in X'
    //             become available to pick.
    //         */
    //         for (int i = 0; i < numX; i++) 
    //             bestX[i] = X[i];
    //     }
    //     else
    //         updateRanks(false);
    // }

    // // Write results
    // ofstream resultsFile;
    // resultsFile.open("results.txt", ios_base::app);
	// for (int i=0; i<numX; i++) 
    //     resultsFile << bestX[i] << ", ";
    
	// resultsFile << timesPopulationSaved << ", " << bestU << ", " << getTime() - ts_start << endl;

    MPI_Finalize();
}

void updateRanks(bool success)
{
    int zeroExists = 0;

    for (int i = 0; i < numX; i++)
    {
        if (success)
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

        if (!success)
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