#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>
#include <sstream>
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
#define ITERS 10000
#define POP_SIZE 50

using namespace std;

/* Configuration */
int numDP   = 10000;       // Vietoviu skaicius (demand points, max 10000)
int numPF   = 5;           // Esanciu objektu skaicius (preexisting facilities)
int numF    = 3;           // Esanciu imoniu skaicius (firms)
int numCL   = 25;          // Kandidatu naujiems objektams skaicius (candidate locations)
constexpr int numX    = 3; // Nauju objektu skaicius

double **demandPoints, *distances;
int *X, *bestX, *ranks;
double bestU;

// Population variables
populationItem* population;
int itemsInPopulation = 0;
int timesPopulationSaved = 0;

// Parallel variables
int id, numProcs, offset, procChunkSize;
double *pop_sendBuff, *pop_recvBuff;
MPI_Datatype population_dt;
MPI_Status *statuses;
MPI_Request ibcastRequest;
MPI_Status status;

//-------

void updateRanks(bool success);

void exchangeFirstSolutions();
void checkForNewNotifications();
void notifyBetterFound();

int main(int argc , char * argv []) {
    MPI_Init(&argc , &argv);
    MPI_Comm_rank(MPI_COMM_WORLD ,&id);
    MPI_Comm_size(MPI_COMM_WORLD ,&numProcs);
    
    double ts_start = getTime();

    //New seed on every run
    srand((unsigned)time(0) * (id + 1));

	loadDemandPoints(numDP, &demandPoints);
	calculateDistancesAsync(numDP, numProcs, id, &distances, demandPoints);

    initPopulationDataTypeToMPI(&population_dt, numX);
    initPopulation(&population, POP_SIZE, numX);
    pop_sendBuff = new double[numX + 1];
    pop_recvBuff = new double[numProcs * (numX + 1)];

    X = new int[numX];
	bestX = new int[numX];
	bestU = -1;
    double u;
    
    //Init ranks
    ranks = new int[numCL];
    for (int i = 0; i < numCL; i++)
        ranks[i] = 1;

    exchangeFirstSolutions();
	
    for (int iters = 0; iters < ITERS; iters++) {
        //printf("iteration - %d \n", iters);
        
        generateSolution_1D(numX, numDP, numCL, X, bestX, ranks, distances, 1, GEN_SOLUTION);

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
            u = evaluateSolution_1D(numX, numDP, numPF, numF, X, demandPoints, distances, 1, EVAL_SOLUTION);
            insert(population, X, numX, u, &itemsInPopulation, POP_SIZE);
        }

        /* Check if other processor found better solution */
        checkForNewNotifications();

        if (u > bestU) 
        {
            updateRanks(true);

            bestU = u;
            for (int i = 0; i < numX; i++) 
                bestX[i] = X[i];

            notifyBetterFound();
        }
        else
            updateRanks(false);
    }

    // Write results
    ofstream resultsFile;
    stringstream fileName;
    fileName << "resultsProc" << id << ".txt" << endl;
    resultsFile.open(fileName.str(), ios_base::app);
	for (int i=0; i<numX; i++) 
        resultsFile << bestX[i] << ", ";
    
	resultsFile << bestU << ", " << timesPopulationSaved << ", " << getTime() - ts_start << endl;
    resultsFile.close();

    MPI_Finalize();
    return 0;
}

void exchangeFirstSolutions()
{
    randomSolution(numCL, numX, X);
    double u = evaluateSolution_1D(numX, numDP, numPF, numF, X, demandPoints, distances, 1, EVAL_SOLUTION);

    // Convert populationItem to array
    // First element - solution, rest - locations
    for (int i = 0; i < numX + 1; i++)
    {
        if (i == 0)
            pop_sendBuff[i] = u;
        else
            pop_sendBuff[i] = (double)X[i - 1];
    }

    MPI_Allgather(pop_sendBuff, 1, population_dt, pop_recvBuff, 1, population_dt, MPI_COMM_WORLD);
    
    int bestIndex = -1;
    // Revert from array to populationItem
    for (int i = 0; i < numProcs; i++)
    {
        population[itemsInPopulation].solution = pop_recvBuff[i * (numX + 1)];
        if (population[itemsInPopulation].solution > bestU)
        {
            bestIndex = i;
            bestU = population[itemsInPopulation].solution;
        }

        for (int j = 0; j < numX; j++)
            population[itemsInPopulation].locations[j] = (int)pop_recvBuff[i * (numX + 1) + (j + 1)];
        
        itemsInPopulation++;
    }

    // All proccesors save the best solution
    for (int i = 0; i < numX; i++) 
        bestX[i] = population[bestIndex].locations[i];
}

void checkForNewNotifications()
{
    for (int i = 0; i < numProcs; i++)
    {
        if (i == id)
            continue;

        int flag;
        MPI_Status bcastStatus;

        MPI_Iprobe(i, 0, MPI_COMM_WORLD, &flag, &status);

        if (flag)
        {
            cout << id << ": got message from " << i << "??" << endl;
        }
    }
}

void notifyBetterFound()
{
    // First element - solution, rest - locations
    for (int i = 0; i < numX + 1; i++)
    {
        if (i == 0)
            pop_sendBuff[i] = bestU;
        else
            pop_sendBuff[i] = (double)bestX[i - 1];
    }

    MPI_Ibcast(pop_sendBuff, 1, population_dt, id, MPI_COMM_WORLD, &ibcastRequest);
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
