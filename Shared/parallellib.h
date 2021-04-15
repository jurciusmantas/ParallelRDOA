void initPopulationStructToMPI(MPI_Datatype* population_dt, int numX);

void calculateDistancesAsync(int numDP, int numProcs, int id, double** distances, double** demandPoints);

double* popItemToArray(populationItem item, int numX);