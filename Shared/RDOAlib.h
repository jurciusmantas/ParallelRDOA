double getTime();
double HaversineDistance(double* a, double* b);
void readConfig(int** params, int count);
void loadDemandPoints(int numDP, double*** demandPoints);
void calculateDistances(int numDP, double*** distances, double** demandPoints);
void updateRanks(int* ranks, int* X, int* bestX, int numCL, int numX, bool success);

void randomSolution(int numCL, int numX, int* X);
void generateSolution(int numX, int numDP, int numCL, int* X, int* bestX, int* ranks, void* distances, int distancesDim, int gen_solution_mode);

double evaluateSolution(int numX, int numDP, int numPF, int numF, int* X, double** demandPoints, void* distances, int distancesDim, int eval_solution_mode);