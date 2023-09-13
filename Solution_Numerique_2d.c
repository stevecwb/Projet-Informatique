#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double energie_du_systeme(double J, int M, int N, double matrix[M][N])
{
    double E = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N-1; j++){
            E = E + matrix[i][j] * matrix[i][j + 1];
            //printf("Energy value for i = %d, j = %d : %lf\n", i ,j, E);
        }
        E = E + matrix[i][0] * matrix[i][N - 1];
        //printf("Energy value = %lf\n", E);
    }
    for (int j = 0; j < N; j++){
        for (int i = 0; i < M-1; i++){
            E = E + matrix[i][j] * matrix[i + 1][j];
            //printf("Energy value for i = %d, j = %d : %lf\n", i ,j, E);
        }
        E = E + matrix[0][j] * matrix[M-1][j];
        //printf("Energy value  = %lf\n", E);
    }
    return -J * E;
}
double magnetization(int M, int N, double matrix[M][N]){
    double Mg = 0;
    for (int i = 0; i < M; i++)
        {
            for(int j = 0; j < N; j++){
                Mg = Mg + matrix[i][j];
            }
        }
    return Mg;
}

int main(){
    int M = 20;
    int N = 25;       // nombre de spins
    int n_iter = 1000000; // nombre d'iterations
    double J = 1;      // constant d'Energie generique
    double kb = 1;     // constant de Boltzman
    double results[4][50];

    // Open a file for writing (you can change "output.txt" to your desired file name)
    FILE *file = fopen("Energy Values Numerical 2d.csv", "w");

    // Check if the file was opened successfully
    if (file == NULL)
    {
        printf("Error opening the file.\n");
        return 1; // Exit with an error code
    }

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

    /*// Display the matrix
    printf("Initial Matrix: \n");
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%lf ", matrice[i][j]);
        }
        printf("\n");
    }*/
    printf("Energy Value = %lf", energie_du_systeme(J, M, N, matrice));

    for (int t = 0; t <= 50; t++)
    {
        double T = t/2.5;
        double matrix[M][N];
        for (int i = 0; i < M; i++)
        {
            for(int j = 0; j < N; j++){
                matrix[i][j] = matrice[i][j];
            }
        }
        /*for (int i = 0; i < 1; i++) {
            for (int j = 0; j < N; j++) {
                // Generate a random number between 0 and 1
                double random = rand()*1.0/RAND_MAX;
                    if(random > 0.75){
                        matrix[i][j] = -1.0;
                    }
                    else{
                    matrix[i][j] = +1.0;
                    }
            }
        }*/

        // printf("Energy Value = %lf", energie_du_systeme(J, N, matrix) / (J * N));
        double energie_moyenne = 0;
        double energie_moyenne_carre = 0;
        double Mg = 0;
        // Monte Carlo Method
        for (int i = 0; i < n_iter; i++)
        {
            int random_spin_lines = rand() % M;
            int random_spin_columns = rand() % N;
            double energy_before = energie_du_systeme(J, M, N, matrix);
            matrix[random_spin_lines][random_spin_columns] *= -1.0; // matrix[0][j] = matrix[0][j]*(-1,0)
            double energy_after = energie_du_systeme(J, M, N, matrix);
            double delta_energy = energy_after - energy_before;
            if (delta_energy > 0)
            {
                double P = exp(- delta_energy / (kb * T)); // introduire la probabilite
                double random = rand() * 1.0 / RAND_MAX;
                if (random > P)
                {
                    matrix[random_spin_lines][random_spin_columns] *= -1.0;
                }
            }
            double energie_aux = energie_du_systeme(J, M, N, matrix);
            energie_moyenne += energie_aux;
            energie_moyenne_carre += pow(energie_aux, 2);
            Mg += magnetization(M, N, matrix);
            /*// Display the matrix
            printf("Matrix changed of iteration %d : ", i);
            for (int i = 0; i < 1; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    printf("%lf ", matrix[i][j]);
                }
                printf("\n");
            }*/
        }
        /*for(int i = 0; i < n_iter; i++){
            for(int j = 0; j < N; j++){
                double energy_before = energie_du_systeme(J, N, matrix);
                printf("energy before %lf\n", energy_before);
                matrix[0][j] *= -1.0; // matrix[0][j] = matrix[0][j]*(-1,0)
                printf("element 1 column %d, value = %lf",j,matrix[0][j]);
                double energy_after = energie_du_systeme(J, N, matrix);
                printf("energy after %lf\n", energy_after);
                double delta_energy = energy_after - energy_before;
                if(delta_energy >0){
                    double P = exp(-delta_energy/(kb*T)); //introduire la probabilite
                    double random = rand();
                    if(random > P){
                        matrix[0][j] *= -1.0;
                    }
                }
            }
        }*/
        energie_moyenne = energie_moyenne / n_iter;
        energie_moyenne_carre = energie_moyenne_carre / n_iter;
        double Cv = (1 / (kb * pow(T, 2))) * (energie_moyenne_carre - pow(energie_moyenne, 2));
        results[0][t] = T;
        results[1][t] = energie_moyenne / (J * M*N);
        results[2][t] = Cv / (kb * M*N);
        results[3][t] = fabs(Mg/(M*N*n_iter));
        // Write the numbers to the file
        fprintf(file, "%lf, %lf, %lf, %lf\n", results[0][t], results[1][t], results[2][t], results[3][t]);

        // printf("\nAverage Energy Value = %lf\n", energie_moyenne / (J * N));
        // printf("\nHeat Capacity = %lf\n", Cv / (kb * N));
        /*// Display the matrix
        for (int i = 0; i < 1; i++) {
            for (int j = 0; j < N; j++) {
                printf("%lf ", matrix[i][j]);
            }
            printf("\n");
        }*/
    }
    // Close the file
    fclose(file);
    return 0;
}