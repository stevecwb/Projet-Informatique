#include <stdio.h>
#include <math.h>

double f_energy(double x, double k){
    return 1.0/(sqrt(1 - pow(k,2)*pow(sin(x),2)));
}

double f_heat_capacity(double x, double k){
    return sqrt(1 - pow(k,2)*pow(sin(x),2));
}

double integral_calculator_energy(double k) {
    double pi = 3.14159265358979323846;
    double n = 1000;
    double h = (pi/2.0 - 0) / n;  // Width of each subinterval
    double sum = f_energy(0, k) + f_energy(pi/2.0, k);  // Initial value of the sum

    for (int i = 1; i < n; i++) {
        double x = 0 + i * h;
        if (i % 2 == 0) {
            sum += 2 * f_energy(x, k);
        } else {
            sum += 4 * f_energy(x, k);
        }
    }

    return h * sum / 3;
}

double integral_calculator_heat_capacity(double k) {
    double pi = 3.14159265358979323846;
    double n = 1000;
    double h = (pi/2.0 - 0) / n;  // Width of each subinterval
    double sum = f_heat_capacity(0, k) + f_heat_capacity(pi/2.0, k);  // Initial value of the sum

    for (int i = 1; i < n; i++) {
        double x = 0 + i * h;
        if (i % 2 == 0) {
            sum += 2 * f_heat_capacity(x, k);
        } else {
            sum += 4 * f_heat_capacity(x, k);
        }
    }

    return h * sum / 3;
}


double analytical_energy_2d(double T, int N){
    double J = 1;
    double kb = 1;
    double beta = 1.0/(kb*T);
    double k = (2*sinh(2*beta*J))/(pow(cosh(2*beta*J),2));
    double pi = 3.14159265358979323846;
    double K1 = integral_calculator_energy(k);
    double E = -2*N*J*tanh(2*beta*J) - N*J*((pow(sinh(2*beta*J),2) - 1)/(sinh(2*beta*J)*cosh(2*beta*J)))*(2*K1/pi - 1) ;
    return E;
}
double analytical_heat_capacity_2d(double T, int N){
    double J = 1;
    double kb = 1;
    double beta = 1.0/(kb*T);
    double k = (2*sinh(2*beta*J))/(pow(cosh(2*beta*J),2));
    double pi = 3.14159265358979323846;
    double K1 = integral_calculator_energy(k);
    double E1 = integral_calculator_heat_capacity(k);
    double Cv = N*kb*(4.0/pi)*pow((beta*J*(1/tanh(2*beta*J))),2)*(K1 - E1 - (1 - pow(tanh(2*beta*J),2))*(pi/2.0 + (2*pow(tanh(2*beta*J),2)-1)*K1));
    return Cv;
}
double analytical_magnatization_2d(double T){
    double J = 1;
    double kb = 1;
    double beta = 1.0/(kb*T);
    if(T < 2.269){
        return pow((1-pow((sinh(2*beta*J)),-4)),(1.0/8.0));
    }
    else{
        return 0;
    }    
}
int main(){
    int i;
    int N = 1;
    double table[4][200];
    // Open a file for writing (you can change "output.txt" to your desired file name)
    FILE *file = fopen("Energy Values 2D 2.csv", "w");
    
    // Check if the file was opened successfully
    if (file == NULL) {
        printf("Error opening the file.\n");
        return 1; // Exit with an error code
    }
    for(i=0; i <=200; i++){
        double T = i/10.0;
        double E = analytical_energy_2d(T, N);
        double Cv = analytical_heat_capacity_2d(T, N);
        double Mg = analytical_magnatization_2d(T);
        table[0][i] = T;
        table[1][i] = E;
        table[2][i] = Cv;
        table[3][i] = Mg;
        printf("%lf\t%lf\t%lf\t%lf\n", table[0][i],table[1][i], table[2][i], table[3][i]);

        // Write the numbers to the file
        fprintf(file, "%lf, %lf, %lf, %lf\n", table[0][i], table[1][i], table[2][i], table[3][i]);
    }

    // Close the file
    fclose(file);

    return 0;
}