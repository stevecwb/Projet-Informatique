/*Solution Exacte en 1 Dimension*/

#include <stdio.h>
#include <math.h>

double analytical_energy_1d(double T){
    /*Expression analytique de l'énergie. Elle reçoit une température et retourne le valeur de l'énergie.*/
    double J = 1; //constant d'énergie generique
    double kb = 1; //constant de Boltzman
    double E = -J*tanh(J / (kb*T));
    return E;
}
double analytical_heat_capacity_1d(double T){
    /*Expression analytique pour la capacité calorifique. Elle reçoit une température et retourne le valeur de Cv.*/
    double J = 1; //constant d'énergie generique
    double kb = 1; //constant de Boltzman
    double Cv = ((pow(J, 2))/(kb*pow(T, 2)))*pow(1.0/cosh((J)/(kb*T)), 2);
    return Cv;
}

int main(){
    int i;
    double table[3][100]; //tableau pour la souvegarde des données
    // Ouvrier un fichier
    FILE *file = fopen("Solution_Analytique_1D.csv", "w");
    
    // Verifie que le fichier s'ouvre
    if (file == NULL) {
        printf("Error opening the file.\n");
        return 1;
    }
    /*Calcul de l'énergie et du Cv avec sauvegarde des données dans un fichier*/
    for(i=0; i <=100; i++){
        double T = i/10.0;
        double E = analytical_energy_1d(T);
        double Cv = analytical_heat_capacity_1d(T);
        table[0][i] = T;
        table[1][i] = E;
        table[2][i] = Cv;
        printf("%lf\t%lf\t%lf\n", table[0][i],table[1][i], table[2][i]);

        // Sauvegarde des données dans un fichier
        fprintf(file, "%lf, %lf, %lf\n", table[0][i], table[1][i], table[2][i]);
    }

    // Fermeture du fichier
    fclose(file);

    return 0;
}

