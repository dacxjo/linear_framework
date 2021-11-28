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
    if (a == NULL || acp == NULL || L == NULL || U == NULL || vector == NULL || y == NULL) {
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
    a[0][1] = 0.5;
    a[0][2] = 0.333333;
    a[0][3] = 0.25;
    a[1][0] = 0.5;
    a[1][1] = 0.333333;
    a[1][2] = 0.25;
    a[1][3] = 0.2;
    a[2][0] = 0.333333;
    a[2][1] = 0.25;
    a[2][2] = 0.2;
    a[2][3] = 0.166667;
    a[3][0] = 0.25;
    a[3][1] = 0.2;
    a[3][2] = 0.166667;
    a[3][3] = 0.142857;

    acp[0][0] = 1;
    acp[0][1] = 0.5;
    acp[0][2] = 0.333333;
    acp[0][3] = 0.25;
    acp[1][0] = 0.5;
    acp[1][1] = 0.333333;
    acp[1][2] = 0.25;
    acp[1][3] = 0.2;
    acp[2][0] = 0.333333;
    acp[2][1] = 0.25;
    acp[2][2] = 0.2;
    acp[2][3] = 0.166667;
    acp[3][0] = 0.25;
    acp[3][1] = 0.2;
    acp[3][2] = 0.166667;
    acp[3][3] = 0.142857;

    gaussLU(n, acp);
    genMatId(n, L);
    genMatNul(n,U);
    luDecompose(n, acp, L, U);
    printf("---------Matriu L---------\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%24.16e\t", L[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("---------Matriu U---------\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%24.16e\t", U[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for (k = 0; k < n; k++) {
        genVectNul(n, vector);
        genVectNul(n,y);
        genVectNul(n,x);
        vector[k] = 1.0;
        printf("Resolent per a Y....\n");
        printf("\n");
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%24.16e\t", L[i][j]);
            }
            printf("= y%d = %24.16e", i + 1, vector[i]);
            printf("\n");
        }
        resoltriinf(n, L, vector, y, tol);
        printf("\n");
        printf("Obtingut el vector:\n");
        printf("\n");
        for (i = 0; i < n; i++) {
            printf("y%d = %24.16e\n", (i + 1), y[i]);
        }
        printf("\n");
        printf("Resolent per a X....\n");
        printf("\n");
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%24.16e\t",U[i][j]);
            }
            printf("= x%d = %24.16e", i + 1, y[i]);
            printf("\n");
        }
        resoltrisup(n,U,y,x,tol);
        printf("\n");
        printf("Obtinguda la solucio:\n");
        for (i = 0; i < n; i++) {
            printf("x%d = %24.16e\n", i + 1, x[i]);
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

