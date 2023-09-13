#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double energie_du_systeme(double J, int M, int N, int O, double matrix[M][N][O])
{
    double E = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N-1; j++)
        {
            for(int k = 0; k < O-1; k++){
                E = E + matrix[i][j][k] * matrix[i][j][k+1];
                //printf("Energy value for i = %d, j = %d, k = %d : %lf\n", i ,j, k, E);
            }
            E = E + matrix[i][j][0] * matrix[i][j][O - 1];
            //printf("Energy value = %lf\n", E);
        }
        E = E + matrix[i][0][0] * matrix[i][N-1][O - 1];
    }
    for (int j = 0; j < N; j++){
        for (int k = 0; k < O-1; k++){
            for (int i = 0; i < M-1; i++){
                E = E + matrix[i][j][k] * matrix[i+1][j][k];
                //printf("Energy value for i = %d, j = %d, k = %d : %lf\n", i ,j, k, E);
            }
            E = E + matrix[0][j][k] * matrix[M-1][j][k];
            //printf("Energy value  = %lf\n", E);
        }
        E = E + matrix[0][j][0] * matrix[M - 1][j][O - 1];
    }
    for (int k = 0; k < O; k++){
        for (int i = 0; i < M-1; i++){
            for (int j = 0; j < N-1; j++){
                E = E + matrix[i][j][k] * matrix[i][j+1][k];
                //printf("Energy value for i = %d, j = %d, k = %d : %lf\n", i ,j, k, E);
            }
            E = E + matrix[i][0][k] * matrix[i][N-1][k];
            //printf("Energy value  = %lf\n", E);
        }
        E = E + matrix[0][0][k] * matrix[M - 1][N - 1][k];
    }
    return -J * E;
}

int main(){
    int M = 2;
    int N = 2;
    int O = 2;     // nombre de spins
    int n_iter = 1000; // nombre d'iterations
    double J = 1;      // constant d'Energie generique
    double kb = 1;     // constant de Boltzman
    double results[4][50];

    // Seed the random number generator with the current time
    srand(time(NULL));
    // Create a matrix of random numbers (-1 or 1)
    double matrice[M][N][O];

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < O; k++)
            {
                // Generate a random number between 0 and 1
                int random_value = rand() % 2;

                // Map 0 to -1 and 1 to 1
                matrice[i][j][k] = (random_value == 0) ? -1.0 : 1.0;
            }
        }
    } // random matrix

    // Display the matrix
    printf("Initial Matrix: \n");
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < O; k++)
            {
                printf("i = %d, j = %d, k = %d, and matrix = %lf\n", i,j,k ,matrice[i][j][k]);
            }
        }
        printf("\n");
    }
    printf("Energy Value = %lf", energie_du_systeme(J, M, N, O, matrice));
    return 0;
}