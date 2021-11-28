#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "funs_linalg.h"

int main() {
    int i, j, k, n, option;
    double **A, *b, tol = 0.0001, totalTime = 0.0;
    FILE *file;
    clock_t inici, final;
    printf("Metode per a resoldre els sistemes:\n");
    printf("[1]-Gauss sense pivotatge:\n");
    printf("[2]-Gauss amb pivotatge:\n");
    scanf("%d", &option);
    if (option == 1) {
        file = fopen("gaussTemps.out", "w");
    } else {
        file = fopen("gaussPivTemps.out", "w");
    }
    for (n = 2; n < 200; n++) {
        A = (double **) malloc(sizeof(double) * n);
        b = (double *) malloc(sizeof(double) * n);
        for (i = 0; i < n; i++) {
            A[i] = (double *) malloc(n * sizeof(double));
            if (A[i] == NULL) {
                printf("No hi ha suficient memòria\n");
                exit(EXIT_FAILURE);
            }
        }
        srand(time(NULL));
        inici = clock();
        for (k = 0; k < 1000; k++) {
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    A[i][j] = rand() / (double) RAND_MAX;
                    b[i] = rand() / (double) RAND_MAX;
                }
            }

            if (option == 1) {
                gauss(n, A, b, tol);
            } else {
                gausspiv(n, A, b, tol);
            }

        }
        final = clock();
        free(A);
        free(b);
        totalTime += (double) (final - inici) / CLOCKS_PER_SEC;
        printf("1000 sistemes %dx%d resolts - Temps total: %24.16e segons\n", n, n, totalTime);
        if (option == 1) {
            file = freopen("gaussTemps.out", "a", file);
        } else {
            file = freopen("gaussPivTemps.out", "a", file);
        }
        fprintf(file, "%d\t %24.16e\n", n, totalTime);
    }
    fclose(file);
    return 0;
}
