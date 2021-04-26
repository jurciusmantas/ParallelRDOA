#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "../Shared/RDOAlib.h"

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
double **demandPoints, **distances;
int *X, *bestX;

int increaseX(int *X, int index, int maxindex);

int main(int argc, char* argv[]) {
	double ts = getTime();

	if (argc != 2)
    {
        cout << "Missing command line parameters or too many provided" << endl;
        abort();
    }

    evalSolution = atoi(argv[1]);
	int* params[5] = { &numDP, &numPF, &numF, &numCL, &numX };
	readConfig(params, 5);

	loadDemandPoints(numDP, &demandPoints);
	calculateDistances(numDP, &distances, demandPoints);
	
	// Sudarom pradini sprendini: [0, 1, 2, 3, ...]
	X = new int[numX];
	bestX = new int[numX];
	for (int i=0; i<numX; i++) {
		X[i] = i;
		bestX[i] = i;
	}
	double u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, 2, evalSolution);
	double bestU = u;
	int iteration = 1;

	while (true) {
		cout << "iteration = " << iteration << endl;
		iteration++;
		if (increaseX(X, numX-1, numCL)) {
			u = evaluateSolution(numX, numDP, numPF, numF, X, demandPoints, distances, 2, evalSolution);
			if (u > bestU) {
				bestU = u;
				for (int i=0; i<numX; i++) bestX[i] = X[i];
			}
		}
		else break;
	}

	// Write results
    ofstream resultsFile;
    resultsFile.open("results.txt", ios_base::app);
	resultsFile << "numX = " << numX << ", numDP = " << numDP << ", evalSol = " << evalSolution << ", numCL = " << numCL << " | " << bestU << " (";
	for (int i=0; i<numX; i++) 
        resultsFile << bestX[i] << ", ";

	resultsFile << ") | time = " << getTime() - ts << endl;
	resultsFile.close();
}

//=============================================================================

int increaseX(int *X, int index, int maxindex) {
	if (X[index]+1 < maxindex-(numX-index-1)) {
		X[index]++;
	}
	else {		 
		if ((index == 0) && (X[index]+1 == maxindex-(numX-index-1))) {
			return 0;
		}
		else {
			if (increaseX(X, index-1, maxindex)) X[index] = X[index-1]+1;
			else return 0;
		}	
	}
	return 1;
}