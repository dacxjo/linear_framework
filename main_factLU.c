#include <stdio.h>
#include <stdlib.h>
#include "funs_linalg.h"

int main() {
    double **a, **acp, **L, **U, *vector, *y, *x, tol = 0.0;
    int i, j, k, n = 4;
    a = (double **) malloc(sizeof(double) * n);
    acp = (double **) malloc(sizeof(double) * n);
    L = (double **) malloc(sizeof(double) * n);
    U = (double **) malloc(sizeof(double) * n);
    vector = (double *) malloc(sizeof(double) * n);
    y = (double *) malloc(sizeof(double) * n);
    x = (double *) malloc(sizeof(double) * n);
    if (a == NULL || acp == NULL || L == NULL || U == NULL || vector == NULL || y == NULL || x == NULL) {
        printf("No hi ha suficient memoria\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < n; i++) {
        a[i] = (double *) malloc(n * sizeof(double));
        acp[i] = (double *) malloc(n * sizeof(double));
        L[i] = (double *) malloc(n * sizeof(double));
        U[i] = (double *) malloc(n * sizeof(double));
        if (a[i] == NULL || acp[i] == NULL || L[i] == NULL || U[i] == NULL) {
            printf("No hi ha suficient memoria\n");
            exit(EXIT_FAILURE);
        }
    }

    a[0][0] = 1;
    a[0][1] = 2;
    a[0][2] = 3;
    a[0][3] = 4;
    a[1][0] = -2;
    a[1][1] = 1;
    a[1][2] = 2;
    a[1][3] = 3;
    a[2][0] = -3;
    a[2][1] = -2;
    a[2][2] = 1;
    a[2][3] = 2;
    a[3][0] = -4;
    a[3][1] = -3;
    a[3][2] = -2;
    a[3][3] = 1;

    acp[0][0] = 1;
    acp[0][1] = 2;
    acp[0][2] = 3;
    acp[0][3] = 4;
    acp[1][0] = -2;
    acp[1][1] = 1;
    acp[1][2] = 2;
    acp[1][3] = 3;
    acp[2][0] = -3;
    acp[2][1] = -2;
    acp[2][2] = 1;
    acp[2][3] = 2;
    acp[3][0] = -4;
    acp[3][1] = -3;
    acp[3][2] = -2;
    acp[3][3] = 1;

  /*  a[0][0] = 4;
    a[0][1] = -2;
    a[0][2] = 1;
    a[1][0] = 20;
    a[1][1] = -7;
    a[1][2] = 12;
    a[2][0] = -8;
    a[2][1] = 13;
    a[2][2] = 17;

    acp[0][0] = 4;
    acp[0][1] = -2;
    acp[0][2] = 1;
    acp[1][0] = 20;
    acp[1][1] = -7;
    acp[1][2] = 12;
    acp[2][0] = -8;
    acp[2][1] = 13;
    acp[2][2] = 17;*/

    gaussLU(n, acp);
    genMatId(n, L);
    luDecompose(n, acp, L, U);
    printf("Matriu L: \n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%e\t", L[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("Matriu U: \n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%e\t", U[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for (k = 0; k < n; k++) {
        genVectNul(n, vector);
        genVectNul(n,y);
        genVectNul(n,x);
        vector[k] = 1.0;
        // Resolver trisup con L y B
        printf("Resolviendo:\n");
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%e\t", L[i][j]);
            }
            printf("= y%d = %e", i + 1, vector[i]);
            printf("\n");
        }
        resoltriinf(n, L, vector, y, tol);
        printf("\n");
        printf("Obtenido el vector:\n");
        for (i = 0; i < n; i++) {
            printf("y%d = %lf\n", (i + 1), y[i]);

        }
        printf("\n");
        printf("Resolviendo:\n");
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%e\t",U[i][j]);
            }
            printf("= x%d = %e", i + 1, y[i]);
            printf("\n");
        }
        resoltrisup(n,U,y,x,tol);
        printf("\n");
        printf("Obtenida la solucion:\n");
        for (i = 0; i < n; i++) {
            printf("x%d = %lf\n", i + 1, x[i]);
        }
        printf("\n");
    }
    free(a);
    free(acp);
    free(U);
    free(L);
    free(vector);
    free(x);
    free(y);
    return (0);
}

