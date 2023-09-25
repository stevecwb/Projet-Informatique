//Solution Numérique en 3 dimensions
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double energie_du_systeme(double J, int M, int N, int O, double matrix[M][N][O])
{       //Calcule l'énergie du systéme en 3D
    double E = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for(int k = 0; k < O-1; k++){
                E = E + matrix[i][j][k] * matrix[i][j][k+1];
            }
            E = E + matrix[i][j][0] * matrix[i][j][O - 1];  //Condition limite sur k
        }
    }
    for (int j = 0; j < N; j++){
        for (int k = 0; k < O; k++){
            for (int i = 0; i < M-1; i++){
                E = E + matrix[i][j][k] * matrix[i+1][j][k];
            }
            E = E + matrix[0][j][k] * matrix[M-1][j][k]; //Condition limite sur i
        }
    }
    for (int k = 0; k < O; k++){
        for (int i = 0; i < M; i++){
            for (int j = 0; j < N-1; j++){
                E = E + matrix[i][j][k] * matrix[i][j+1][k];
            }
            E = E + matrix[i][0][k] * matrix[i][N-1][k]; //Condition limite sur j
        }
    }
    return -J * E;
}
double aimantation(int M, int N, int O, double matrix[M][N][O]){    //Calcule l'aimantation du systéme
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
    // nombre de spins
    int M = 10;
    int N = 10;
    int O = 5;
    int n_iter = 750000; // nombre d'itérations
    double J = 1;      // constante d'énergie generique
    double kb = 1;     // constante de Boltzman
    double resultats[4][100];

    FILE *file = fopen("Solution_Numerique_3D.csv", "w");

    if (file == NULL)   //Vérifie que le fichier s'ouvre
    {
        printf("Error opening the file.\n");
        return 1;
    }

    srand(time(NULL));
    double matrice[M][N][O];    //Création aléatoire de la matrice, comme en dimension 1 et 2

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < O; k++)
            {
                int random_value = rand() % 2;

                matrice[i][j][k] = (random_value == 0) ? -1.0 : 1.0;
            }
        }
    }


    printf("Valeur de l'2nergie = %lf", energie_du_systeme(J, M, N, O, matrice));

    for (int t = 0; t <= 100; t++)      //Calculs en fonction de T
    {
        double T = t/10.0;
        double matrix[M][N][O];
        for (int i = 0; i < M; i++)
        {
            for(int j = 0; j < N; j++)
            {
                for(int k = 0; k < O; k++)
                {
                    matrix[i][j][k] = matrice[i][j][k];     //Copie de la matrice à modifier
                }
            }
        }

        double energie_moyenne = 0;
        double energie_moyenne_carre = 0;
        double Mg = 0;
        for (int i = 0; i < n_iter; i++)        //Méthode de Monte Carlo
        {
            int random_spin_lignes = rand() % M;
            int random_spin_colonnes = rand() % N;
            int random_spin_profondeur = rand() % O;
            double energie_avant = energie_du_systeme(J, M, N, O, matrix);
            matrix[random_spin_lignes][random_spin_colonnes][random_spin_profondeur] *= -1.0;
            double energie_apres = energie_du_systeme(J, M, N, O, matrix);
            double delta_energie = energie_apres - energie_avant;
            if (delta_energie > 0)
            {
                double P = exp(- delta_energie / (kb * T)); // Introduire la probabilité
                double random = rand() * 1.0 / RAND_MAX;
                if (random > P)
                {
                    matrix[random_spin_lignes][random_spin_colonnes][random_spin_profondeur] *= -1.0;
                }
            }
            double energie_aux = energie_du_systeme(J, M, N, O, matrix);
            energie_moyenne += energie_aux;
            energie_moyenne_carre += pow(energie_aux, 2);
            Mg += aimantation(M, N, O, matrix);
        }

        energie_moyenne = energie_moyenne / n_iter;
        energie_moyenne_carre = energie_moyenne_carre / n_iter;
        double Cv = (1 / (kb * pow(T, 2))) * (energie_moyenne_carre - pow(energie_moyenne, 2));
        resultats[0][t] = T;
        resultats[1][t] = energie_moyenne / (J * M*N*O);
        resultats[2][t] = Cv / (kb * M*N*O);
        resultats[3][t] = fabs(Mg/(M*N*O*n_iter));
        fprintf(file, "%lf, %lf, %lf, %lf\n", resultats[0][t], resultats[1][t], resultats[2][t], resultats[3][t]);


    }
    fclose(file);
    return 0;
}