#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double energie_du_systeme(double J, int M, int N, double matrix[M][N])
{
    double E = 0;
    for (int i = 0; i < M - 1; i++)
    {
        for (int j = 0; j < N-1; j++){
            E = E + matrix[i][j] * matrix[i][j + 1];
            E = E + matrix[i][j] * matrix[i + 1][j];
            E = E + matrix[i][0] * matrix[i][N - 1];
            E = E + matrix[0][j] * matrix[M-1][j];
        }
    }
    return -J * E;
}

int main(){
    int M = 1;
    int N = 1;       // nombre de spins
    int n_iter = 150000; // nombre d'iterations
    double J = 1;      // constant d'Energie generique
    double kb = 1;     // constant de Boltzman
    double results[3][20];
    // Seed the random number generator with the current time
    srand(time(NULL));
    // Create a matrix of random numbers (-1 or 1)
    double matrice[M][N];

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            // Generate a random number between 0 and 1
            int random_value = rand() % 2;

            // Map 0 to -1 and 1 to 1
            matrice[i][j] = (random_value == 0) ? -1.0 : 1.0;
        }
    } // random matrix

    // Display the matrix
    printf("Initial Matrix: \n");
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%lf ", matrice[i][j]);
        }
        printf("\n");
    }
    printf("Energy Value = %lf", energie_du_systeme(J, M, N, matrice));
    return 0;
}