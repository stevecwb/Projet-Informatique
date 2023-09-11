#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double energie_du_systeme(double J, int N, double matrix[1][N]){
    double E = 0;
    for(int i = 0; i < N-1; i++ ){
        E = E + matrix[0][i]*matrix[0][i+1];
    }
    return -J*E;
}
int main(){
    int N = 500; // nombre de spins
    int n_iter = 100; // nombre d'iterations
    double J = 1.380649*pow(10,-23); // constant d'Energie generique
    double T = 4.0; // K Temperature
    double kb = 1.380649*pow(10,-23); //constant de Boltzman

    // Seed the random number generator with the current time
    srand(time(NULL));
    // Create a matrix of random numbers (-1 or 1)
    double matrix[1][N];

    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < N; j++) {
            // Generate a random number between 0 and 1
            int random_value = rand() % 2;
            
            // Map 0 to -1 and 1 to 1
            matrix[i][j] = (random_value == 0) ? -1.0 : 1.0;
        }
    }

    // Display the matrix
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }

    printf("Energy Value = %lf", energie_du_systeme(J, N, matrix)/(J*N));

    // Monte Carlo Method
    for(int i = 0; i < n_iter; i++){
        int random_spin = rand() % N;
        double energy_before = energie_du_systeme(J, N, matrix);
        matrix[0][random_spin] *= -1.0; // matrix[0][j] = matrix[0][j]*(-1,0)
        double energy_after = energie_du_systeme(J, N, matrix);
        double delta_energy = energy_after - energy_before;
        if(delta_energy > 0){
                double P = exp(-delta_energy/(kb*T)); //introduire la probabilite
                double random = rand()*1.0/RAND_MAX;
                if(random > P){
                    matrix[0][random_spin] *= -1.0;
                }
            }
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
    printf("\nEnergy Value Final = %lf\n", energie_du_systeme(J, N, matrix)/(J*N));
    // Display the matrix
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
    return 0;
}