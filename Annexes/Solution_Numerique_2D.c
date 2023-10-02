//Solution Numerique 2 dimensions
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


double energie_du_systeme(double J, int M, int N, double matrix[M][N])
{
    //Fonction permettant le calcul de l'energie du systeme en 2D, avec en entree J dependant de T,
    //la taille de la matrice NxM, et la matrice representant le systeme
    double E = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N-1; j++){
            E = E + matrix[i][j] * matrix[i][j + 1];    //Ajout de la 1re interaction du spin

        }
        E = E + matrix[i][0] * matrix[i][N - 1];    //Ajout de la condition limite selon les colonnes

    }
    for (int j = 0; j < N; j++){
        for (int i = 0; i < M-1; i++){
            E = E + matrix[i][j] * matrix[i + 1][j];    //Ajout de la 2nde interaction du spin
        }
        E = E + matrix[0][j] * matrix[M-1][j];      //Ajout de la condition limite selon les lignes

    }
    return -J * E;
}
double aimantation(int M, int N, double matrix[M][N]){
    //Fonction permettant le calcul de l'aimantation en fonction des dimensions du systeme et du systeme
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
    int M = 20;       //nombre de colonnes
    int N = 25;       // nombre de lignes
    int n_iter = 500000; // nombre d'iterations
    double J = 1;      // constante d'Energie generique
    double kb = 1;     // constante de Boltzmann
    double resultats[4][100]; //Tableau des resultats

    FILE *file = fopen("Solution_Numerique_2D.csv", "w");


    if (file == NULL)   //Verifie l'ouverture du fichier
    {
        printf("Error opening the file.\n");
        return 1;
    }

    srand(time(NULL));

    double matrice[M][N]; //Creation de la matrice aleatoire representant le systeme

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {

            int nombre_aleatoire = rand() % 2; //Genere un nombre aleatoire entre 0 et 1


            matrice[i][j] = (nombre_aleatoire == 0) ? -1.0 : 1.0; //Change les 0 en -1
        }
    }


    printf("Valeur de l'Energie = %lf", energie_du_systeme(J, M, N, matrice)); //Energie initiale du systeme

    for (int t = 0; t <= 100; t++)  //Boucle sur la temperature : 100 points pour une temperature entre 0 et 10 K
    {
        double T = t/10.0;
        double matrix[M][N]; //Copie � modifier du systeme initial
        for (int i = 0; i < M; i++)
        {
            for(int j = 0; j < N; j++){
                matrix[i][j] = matrice[i][j];
            }
        }

        double energie_moyenne = 0;         //Valeurs initiales nulles
        double energie_moyenne_carre = 0;
        double Mg = 0;
        for (int i = 0; i < n_iter; i++)    //Methode de Monte Carlo
        {
            int aleatoire_ligne = rand() % M;         //Changement aleatoire d'un spin
            int aleatoire_colonne = rand() % N;
            double energie_avant = energie_du_systeme(J, M, N, matrix);
            matrix[aleatoire_ligne][aleatoire_colonne] *= -1.0;         //Inversion du spin choisi
            double energie_apres = energie_du_systeme(J, M, N, matrix);
            double delta_energie = energie_apres - energie_avant;     //Calcul de la difference d'energie avant et apres inversion du spin
            if (delta_energie > 0)
            {
                double P = exp(- delta_energie / (kb * T)); // Calcul de la probabilite, dependante de T
                double random = rand() * 1.0 / RAND_MAX;
                if (random > P)
                {
                    matrix[aleatoire_ligne][aleatoire_colonne] *= -1.0;
                }
            }

            double energie_aux = energie_du_systeme(J, M, N, matrix);
            energie_moyenne += energie_aux;         //Ajout de l'energie afin de calculer la moyenne
            energie_moyenne_carre += pow(energie_aux, 2);
            Mg += aimantation(M, N, matrix);
        }
        energie_moyenne = energie_moyenne / n_iter;
        energie_moyenne_carre = energie_moyenne_carre / n_iter;
        double Cv = (1 / (kb * pow(T, 2))) * (energie_moyenne_carre - pow(energie_moyenne, 2));     //Calcul de Cv. La methode par derivee se fait directement sur Python
        resultats[0][t] = T;
        resultats[1][t] = energie_moyenne / (J * M*N);
        resultats[2][t] = Cv / (kb * M*N);
        resultats[3][t] = fabs(Mg/(M*N*n_iter));
        fprintf(file, "%lf, %lf, %lf, %lf\n", resultats[0][t], resultats[1][t], resultats[2][t], resultats[3][t]);      //Implementation des resultats dans le tableau

    }
    fclose(file);
    return 0;
}