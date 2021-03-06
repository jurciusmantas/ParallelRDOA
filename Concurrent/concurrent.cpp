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
	numX;       // new location count

/* Algorithm variables */
double **demandPoints, **distances;
int *X, *bestX, *ranks;

int main(int argc, char* argv[]) {
    double ts_start = getTime();

    //New seed on every run
    srand((unsigned)time(0));

    if (argc != 4)
    {
        cout << "Missing command line parameters or too many provided" << endl;
        abort();
    }

    genSolution = atoi(argv[1]);
    evalSolution = atoi(argv[2]);
    iterations = atoi(argv[3]);

    int* params[5] = { &numDP, &numPF, &numF, &numCL, &numX };
	readConfig(params, 5);
    
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
    double u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, 2, evalSolution);
    bestU = u;
    for (int i = 0; i < numX; i++) 
        bestX[i] = X[i];
	
    for (int iters = 0; iters < iterations; iters++) {
        generateSolution(numX, numDP, numCL, X, bestX, ranks, distances, 2, genSolution);
        u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, 2, evalSolution);

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
    fileName << "results" << genSolution << evalSolution << ".txt";
    resultsFile.open(fileName.str(), ios_base::app);
	for (int i=0; i<numX; i++) 
        resultsFile << bestX[i] << " ";

	resultsFile << ", " << bestU << ", " << getTime() - ts_start << endl;
}
