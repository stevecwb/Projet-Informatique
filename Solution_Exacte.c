#include <stdio.h>
#include <math.h>

double analytical_energy_1d(double T){
    double J = 1;
    double kb = 1;
    double E = -J*tanh(J / (kb*T));
    return E;
}
double analytical_heat_capacity_1d(double T){
    double J = 1;
    double kb = 1;
    double Cv = ((pow(J, 2))/(kb*pow(T, 2)))*pow(1.0/cosh((J)/(kb*T)), 2);
    return Cv;
}

int main(){
    int i;
    double table[3][100];
    // Open a file for writing (you can change "output.txt" to your desired file name)
    FILE *file = fopen("1D Analytical.csv", "w");
    
    // Check if the file was opened successfully
    if (file == NULL) {
        printf("Error opening the file.\n");
        return 1; // Exit with an error code
    }
    for(i=0; i <=100; i++){
        double T = i/10.0;
        double E = analytical_energy_1d(T);
        double Cv = analytical_heat_capacity_1d(T);
        table[0][i] = T;
        table[1][i] = E;
        table[2][i] = Cv;
        printf("%lf\t%lf\t%lf\n", table[0][i],table[1][i], table[2][i]);

        // Write the numbers to the file
        fprintf(file, "%lf, %lf, %lf\n", table[0][i], table[1][i], table[2][i]);
    }

    // Close the file
    fclose(file);

    return 0;
}

