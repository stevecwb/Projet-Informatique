//Solution Numérique 1 dimension
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double energie_du_systeme(double J, int N, double matrix[1][N])
{         //Fonction donnant l'énergie du systéme, en fonction de J, de la taille du systéme et de l'état de celui-ci
    double E = 0;
    for (int i = 0; i < N - 1; i++)
    {
        E = E + matrix[0][i] * matrix[0][i + 1];
    }
    E = E + matrix[0][0] * matrix[0][N - 1]; // Condition limite : interaction entre le premier et le dernier spin
    return -J * E;
}
int main()
{
    int N = 500;       // nombre de spins
    int n_iter = 250000; // nombre d'itérations
    double J = 1;      // constante d'Energie generique
    double kb = 1;     // constante de Boltzmann
    double resultats[3][100];


    FILE *file = fopen("Solution_Numerique_1D.csv", "w");


    if (file == NULL) //Vérifie que le fichier s'ouvre bien
    {
        printf("Error opening the file.\n");
        return 1;
    }
    srand(time(NULL));
    double matrice[1][N];       //Création aléatoire de la matrice

    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < N; j++)
        {
            int random_value = rand() % 2;  //Génération aléatoire d'un nombre entre 0 et 1

            matrice[i][j] = (random_value == 0) ? -1.0 : 1.0;   //Change les 0 avec -1
        }
    }


    for (int t = 0; t <= 100; t++)  //boucle sur la température
    {
        double T = t/10.0;
        double matrix[1][N];
        for (int i = 0; i < N; i++)
        {
            matrix[0][i] = matrice[0][i];   //Remplissage de la copie à modifier
        }

        double energie_moyenne = 0;
        double energie_moyenne_carre = 0;
        for (int i = 0; i < n_iter; i++)    //Méthode de Monte Carlo
        {
            int random_spin = rand() % N;   //Choix aléatoire du spin à modifier
            double energie_avant = energie_du_systeme(J, N, matrix);
            matrix[0][random_spin] *= -1.0; //Inversion du spin
            double energie_apres = energie_du_systeme(J, N, matrix);
            double delta_energie = energie_apres - energie_avant;       //Calcul de la différence d'énergie
            if (delta_energie > 0)
            {
                double P = exp(- delta_energie / (kb * T)); // Introduire la probabilite selon Boltzmann
                double random = rand() * 1.0 / RAND_MAX;
                if (random > P)
                {
                    matrix[0][random_spin] *= -1.0;
                }
            }
            double energie_aux = energie_du_systeme(J, N, matrix);
            energie_moyenne += energie_aux;
            energie_moyenne_carre += pow(energie_aux, 2);

        }

        energie_moyenne = energie_moyenne / n_iter;
        energie_moyenne_carre = energie_moyenne_carre / n_iter;
        double Cv = (1 / (kb * pow(T, 2))) * (energie_moyenne_carre - pow(energie_moyenne, 2));     //Calcul numérique de la capacité calorifique
        resultats[0][t] = T;
        resultats[1][t] = energie_moyenne / (J * N);
        resultats[2][t] = Cv / (kb * N);
        fprintf(file, "%lf, %lf, %lf\n", resultats[0][t], resultats[1][t], resultats[2][t]);


    }

    fclose(file);
    return 0;
}