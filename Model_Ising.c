#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

double energie_du_systeme_1d(double J, int N, double matrix[1][N])
{
    double E = 0;
    for (int i = 0; i < N - 1; i++)
    {
        E = E + matrix[0][i] * matrix[0][i + 1];
    }
    E = E + matrix[0][0] * matrix[0][N - 1]; // boundary conditions
    return -J * E;
}
int methode_ising_1d(int N, int n_iter, double J, double kb)
{
    // nombre de spins
    // nombre d'iterations
    // constant d'Energie generique
    // constant de Boltzman
    double results[3][20];

    // Open a file for writing (you can change "output.txt" to your desired file name)
    FILE *file = fopen("Resultats por 1 dimension.csv", "w");

    // Check if the file was opened successfully
    if (file == NULL)
    {
        printf("Error opening the file.\n");
        return 1; // Exit with an error code
    }
    // Seed the random number generator with the current time
    srand(time(NULL));
    // Create a matrix of random numbers (-1 or 1)
    double matrice[1][N];

    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < N; j++)
        {
            // Generate a random number between 0 and 1
            int random_value = rand() % 2;

            // Map 0 to -1 and 1 to 1
            matrice[i][j] = (random_value == 0) ? -1.0 : 1.0;
        }
    } // random matrix

    /*// Display the matrix
    printf("Initial Matrix: ");
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%lf ", matrice[i][j]);
        }
        printf("\n");
    }*/
    for (int T = 0; T <= 20; T++)
    {
        double matrix[1][N];
        for (int i = 0; i < N; i++)
        {
            matrix[0][i] = matrice[0][i];
        }

        // printf("Energy Value = %lf", energie_du_systeme(J, N, matrix) / (J * N));
        double energie_moyenne = 0;
        double energie_moyenne_carre = 0;
        // Monte Carlo Method
        for (int i = 0; i < n_iter; i++)
        {
            int random_spin = rand() % N;
            double energy_before = energie_du_systeme_1d(J, N, matrix);
            matrix[0][random_spin] *= -1.0; // matrix[0][j] = matrix[0][j]*(-1,0)
            double energy_after = energie_du_systeme_1d(J, N, matrix);
            double delta_energy = energy_after - energy_before;
            if (delta_energy > 0)
            {
                double P = exp(- delta_energy / (kb * T)); // introduire la probabilite
                double random = rand() * 1.0 / RAND_MAX;
                if (random > P)
                {
                    matrix[0][random_spin] *= -1.0;
                }
            }
            double energie_aux = energie_du_systeme_1d(J, N, matrix);
            energie_moyenne += energie_aux;
            energie_moyenne_carre += pow(energie_aux, 2);
        }

        energie_moyenne = energie_moyenne / n_iter;
        energie_moyenne_carre = energie_moyenne_carre / n_iter;
        double Cv = (1 / (kb * pow(T, 2))) * (energie_moyenne_carre - pow(energie_moyenne, 2));
        results[0][T] = T;
        results[1][T] = energie_moyenne / (J * N);
        results[2][T] = Cv / (kb * N);
        // Write the numbers to the file
        fprintf(file, "%lf, %lf, %lf\n", results[0][T], results[1][T], results[2][T]);

    }
    // Close the file
    fclose(file);
    return 0;
}


int main(){
    double J = 1;
    double kb = 1;
    char dimension[2];
    printf("Dans quelle dimension souhaitez-vous simuler (1D, 2D ou 3D) ?");
    scanf("%s", &dimension);
    if(strcmp(dimension, "1D") == 0){
        int N;
        printf("En 1 dimension, vous devez simplement saisir le nombre de spins : ");
        scanf("%d", &N);
        int n_iter = 500*N;
        methode_ising_1d(N, n_iter, J, kb);
    }
    else if(strcmp(dimension, "2D") == 0){
        printf("2D");
    }
    else if(strcmp(dimension, "3D") == 0){
        printf("3D");
    }
    else{
        printf("Vous devez avoir mal tape. Reessayez avec les options suivantes : 1D, 2D et 3D.");
        return 0;
    }
    //printf("cets %s", dimension);
    //methode_ising_1d();
    return 0;
}