double getTime();
double HaversineDistance(double* a, double* b);
void loadDemandPoints(int numDP, double*** demandPoints);
void calculateDistances(int numDP, double*** distances, double** demandPoints);

void randomSolution(int numCL, int numX, int* X);
int rouletteWheel(double* probabilities, int count);
void generateSolution(int numX, int numCL, int* X, int* bestX, int* ranks, double** distances, int gen_solution_mode);

double evaluateSolution(int numX, int numDP, int numPF, int numF, int* X, double** demandPoints, double** distances, int eval_solution_mode);