#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double energie_du_systeme(double J, int M, int N, int O, double matrix[M][N][O])
{
    double E = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for(int k = 0; k < O-1; k++){
                E = E + matrix[i][j][k] * matrix[i][j][k+1];
                //printf("Energy value for i = %d, j = %d, k = %d : %lf\n", i ,j, k, E);
            }
            E = E + matrix[i][j][0] * matrix[i][j][O - 1];
            //printf("Energy value for i = %d, j = %d, = %lf\n", i, j, E);
        }
    }
    for (int j = 0; j < N; j++){
        for (int k = 0; k < O; k++){
            for (int i = 0; i < M-1; i++){
                E = E + matrix[i][j][k] * matrix[i+1][j][k];
                //printf("Energy value for i = %d, j = %d, k = %d : %lf\n", i ,j, k, E);
            }
            E = E + matrix[0][j][k] * matrix[M-1][j][k];
            //printf("Energy value for j = %d, k = %d, = %lf\n", j, k, E);
        }
    }
    for (int k = 0; k < O; k++){
        for (int i = 0; i < M; i++){
            for (int j = 0; j < N-1; j++){
                E = E + matrix[i][j][k] * matrix[i][j+1][k];
                //printf("Energy value for i = %d, j = %d, k = %d : %lf\n", i ,j, k, E);
            }
            E = E + matrix[i][0][k] * matrix[i][N-1][k];
            //printf("Energy value for i = %d, k = %d, = %lf\n", i, k, E);
        }
    }
    return -J * E;
}
double magnetization(int M, int N, int O, double matrix[M][N][O]){
    double Mg = 0;
    for (int i = 0; i < M; i++)
        {
            for(int j = 0; j < N; j++){
                for(int k =0; k < O; k++)
                {
                    Mg = Mg + matrix[i][j][k];
                }
            }
        }
    return Mg;
}

int main(){
    int M = 10;
    int N = 10;
    int O = 5;     // nombre de spins
    int n_iter = 750000; // nombre d'iterations
    double J = 1;      // constant d'Energie generique
    double kb = 1;     // constant de Boltzman
    double results[4][100];

    // Open a file for writing (you can change "output.txt" to your desired file name)
    FILE *file = fopen("3D Numerical.csv", "w");

    // Check if the file was opened successfully
    if (file == NULL)
    {
        printf("Error opening the file.\n");
        return 1; // Exit with an error code
    }

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
    /*printf("Initial Matrix: \n");
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
    }*/
    printf("Energy Value Global = %lf", energie_du_systeme(J, M, N, O, matrice));

    for (int t = 0; t <= 100; t++)
    {
        double T = t/10.0;
        double matrix[M][N][O];
        for (int i = 0; i < M; i++)
        {
            for(int j = 0; j < N; j++)
            {
                for(int k = 0; k < O; k++)
                {
                    matrix[i][j][k] = matrice[i][j][k];
                }
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
            int random_spin_width = rand() % O;
            double energy_before = energie_du_systeme(J, M, N, O, matrix);
            matrix[random_spin_lines][random_spin_columns][random_spin_width] *= -1.0; // matrix[0][j] = matrix[0][j]*(-1,0)
            double energy_after = energie_du_systeme(J, M, N, O, matrix);
            double delta_energy = energy_after - energy_before;
            if (delta_energy > 0)
            {
                double P = exp(- delta_energy / (kb * T)); // introduire la probabilite
                double random = rand() * 1.0 / RAND_MAX;
                if (random > P)
                {
                    matrix[random_spin_lines][random_spin_columns][random_spin_width] *= -1.0;
                }
            }
            double energie_aux = energie_du_systeme(J, M, N, O, matrix);
            energie_moyenne += energie_aux;
            energie_moyenne_carre += pow(energie_aux, 2);
            Mg += magnetization(M, N, O, matrix);
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
        results[1][t] = energie_moyenne / (J * M*N*O);
        results[2][t] = Cv / (kb * M*N*O);
        results[3][t] = fabs(Mg/(M*N*O*n_iter));
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