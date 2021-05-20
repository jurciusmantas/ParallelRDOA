#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>
#include <sstream>
#include "../Shared/RDOAlib.h"
#include "../Shared/parallellib.h"

/*  
    EvalSolution
    0 - BinaryRule, 1 - PartialyBinaryRule 
*/
using namespace std;

/* Command line parameters */
int evalSolution = 0;

/* Configuration */
int numDP,      // demand point locations count, max 10000
	numPF,      // preexisting facilities count
	numF,       // preexisting firm count
	numCL,      // candidate locations count
	numX;       // new location count

/* Algorithm variables */
double **demandPoints, *distances;
int *X, *bestX;
double bestU;

/* Parallel variables */
int id, numProcs, offset, procChunkSize;
double *pop_sendBuff, *pop_recvBuff;
MPI_Datatype population_dt;

bool increaseX(int *X, int index, int maxindex);
void getAndEvaluateResults(int numProcs);
void getBestSolutions();

int main(int argc , char * argv []) {
	MPI_Init(&argc , &argv);
    MPI_Comm_rank(MPI_COMM_WORLD ,&id);
    MPI_Comm_size(MPI_COMM_WORLD ,&numProcs);

	double ts_start = getTime();

    if (argc != 2)
    {
        cout << "Missing command line parameters or too many provided" << endl;
        MPI_Finalize();
        abort();
    }

    evalSolution = atoi(argv[1]);

    int* params[5] = { &numDP, &numPF, &numF, &numCL, &numX };
	readConfig(params, 5);

	loadDemandPoints(numDP, &demandPoints);
    calculateDistancesAsync(numDP, numProcs, id, &distances, demandPoints);

    initPopulationStructToMPI(&population_dt, numX);
    pop_sendBuff = new double[numX + 1];
    pop_recvBuff = new double[numProcs * (numX + 1)];

    X = new int[numX];
    bestX = new int[numX];
	bestU = -1;
    double u;

    // Init solution: [0, 1, 2, 3, ...]
    for (int i=0; i<numX; i++) {
        X[i] = i;
        bestX[i] = i;
    }
    
    // Offset by id
    for (int i = 0; i < id; i++)
        increaseX(X, numX - 1, numCL);
    
    bestU = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, 1, evalSolution);

    int i = 0;
    while (true) {
        cout << id << ": iteration - " << i << endl;
        bool breakMain = false;
        // Offset by numProcs
        for (int i = 0; i < numProcs; i++)
        {
            if (!increaseX(X, numX - 1, numCL))
            {
                breakMain = true;
                break;
            }
        }

        if (breakMain)
            break;

        u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, 1, evalSolution);
        if (u > bestU) 
        {
            bestU = u;
            for (int i=0; i<numX; i++) 
                bestX[i] = X[i];
        }

        i++;
    }

    getBestSolutions();

    // Write results
    if (id == 0)
    {
        ofstream resultsFile;
        stringstream fileName;
        fileName << "results" << evalSolution << "-n" << numProcs << ".txt";
        resultsFile.open(fileName.str(), ios_base::app);
        for (int i=0; i<numX; i++) 
            resultsFile << bestX[i] << " ";
        
        resultsFile << ", " << bestU << ", " << getTime() - ts_start << endl;
        resultsFile.close();
    }

    MPI_Finalize();
    return 0;
}

void getBestSolutions()
{
    for (int i = 0; i < numX + 1; i++)
    {
        if (i == 0)
            pop_sendBuff[i] = bestU;
        else
            pop_sendBuff[i] = (double)bestX[i - 1];
    }

    MPI_Gather(pop_sendBuff, 1, population_dt, pop_recvBuff, 1, population_dt, 0, MPI_COMM_WORLD);

    // Find best from data
    for (int i = 0; i < numProcs; i++)
    {
        if (pop_recvBuff[i * (numX + 1)] > bestU)
        {
            bestU = pop_recvBuff[i * (numX + 1)];
            for (int j = 0; j < numX; j++)
                bestX[j] = (int)pop_recvBuff[i * (numX + 1) + (j + 1)];
        }
    }
}

bool increaseX(int *X, int index, int maxindex) {
	if (X[index]+1 < maxindex-(numX-index-1)) 
		X[index]++;
	
	else 
    {		 
		if ((index == 0) && (X[index]+1 == maxindex-(numX-index-1))) 
			return false;
		
		else 
        {
			if (increaseX(X, index-1, maxindex)) 
                X[index] = X[index-1]+1;
			else 
                return false;
		}	
	}
	return true;
}