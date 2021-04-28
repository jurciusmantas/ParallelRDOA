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
using namespace std;

/* Command line parameters */
int genSolution = 0, 
    evalSolution = 0,
    iterations = 0;

/* Configuration */
int numDP,      // demand point locations count, max 10000
	numPF,      // preexisting facilities count
	numF,       // preexisting firm count
	numCL,      // candidate locations count
	numX,       // new location count
	popSize;    // population size

/* Algorithm variables */
double **demandPoints, *distances;
int *X, *bestX, *ranks;
double bestU;

/* Population variables */
populationItem* population;
int itemsInPopulation = 0;
int timesPopulationSaved = 0;

/* Parallel variables */
int id, numProcs, offset, procChunkSize;
double *pop_sendBuff, *pop_recvBuff;
MPI_Datatype population_dt;
MPI_Status *statuses;
MPI_Request ibcastRequest;
MPI_Status status;

//-------

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

    genSolution = atoi(argv[1]);
    evalSolution = atoi(argv[2]);
    iterations = atoi(argv[3]);

    int* params[6] = { &numDP, &numPF, &numF, &numCL, &numX, &popSize };
	readConfig(params, 6);

	loadDemandPoints(numDP, &demandPoints);
	calculateDistancesAsync(numDP, numProcs, id, &distances, demandPoints);

    initPopulationStructToMPI(&population_dt, numX);
    initPopulation(&population, popSize, numX);
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
	
    for (int iters = 0; iters < iterations; iters++) {
        printf("iteration - %d \n", iters);
        
        generateSolution(numX, numDP, numCL, X, bestX, ranks, distances, 1, genSolution);

        //Search for solution in population
        populationItem popItem = search(population, popSize, X, numX);
        if (popItem.solution > -1.0)
        {
            /* Generated solution was found in population */
            u = popItem.solution;
            timesPopulationSaved++;
        }
        else
        {
            /* Generated solution was not found in population */
            u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, 1, evalSolution);
            insert(population, X, numX, u, &itemsInPopulation, popSize);
        }

        /* Check if other processor found better solution */
        checkForNewNotifications();

        if (u > bestU) 
        {
            updateRanks(ranks, X, bestX, numCL, numX, true);

            bestU = u;
            for (int i = 0; i < numX; i++) 
                bestX[i] = X[i];

            notifyBetterFound();
        }
        else
            updateRanks(ranks, X, bestX, numCL, numX, false);
    }

    // Write results
    ofstream resultsFile;
    stringstream fileName;
    fileName << "resultsProc" << id << ".txt";
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
    double u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, 1, evalSolution);

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
