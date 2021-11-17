#include <stdlib.h>
#include <stdio.h>
#include "funs_linalg.h"

int main() {
    int i, j, n, is_solved,is_exact = 0;
    double **A, **aTemp, *b, *bTemp,*r, tol = 0.0;
    printf("Ingressi la dimensió de la matriu:\n");
    scanf("%d", &n);
    A = (double **) malloc(sizeof(double) * n);
    aTemp = (double **) malloc(sizeof(double) * n);
    b = (double *) malloc(sizeof(double) * n);
    bTemp = (double *) malloc(sizeof(double) * n);
    r = (double *) malloc(sizeof(double) * n);
    if (A == NULL || aTemp == NULL || b == NULL || r == NULL || bTemp == NULL) {
        printf("No hi ha suficient memòria\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < n; i++) {
        A[i] = (double *) malloc(n * sizeof(double));
        aTemp[i] = (double *) malloc(n * sizeof(double));
        if (A[i] == NULL) {
            printf("No hi ha suficient memòria\n");
            exit(EXIT_FAILURE);
        }
        if (aTemp[i] == NULL) {
            printf("No hi ha suficient memòria\n");
            exit(EXIT_FAILURE);
        }
    }
    printf("Ingressi els elements de la matriu:\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            scanf("%lf", &A[i][j]);
            aTemp[i][j] = A[i][j];
        }
    }
    printf("Ingressi els elements del vector B:\n");
    for (i = 0; i < n; ++i) {
        scanf("%lf", &b[i]);
        bTemp[i] = b[i];
    }
    printf("---------------Matriu----------------\n");
    printf("\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            printf("%lf\t", A[i][j]);
        }
        printf("|\t %lf", b[i]);
        printf("\n");
    }
    printf("\n");
    is_solved = gauss(n, A, b, tol);
    if (is_solved == 0) {
        printf("El sistema s'ha resolt: Status %d\n", is_solved);
        printf("-----------Matriu de gauss------------\n");
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                printf("%lf\t", A[i][j]);
            }
            printf("\n");
        }
        printf("\n");
        printf("-----------Solucions----------\n");
        for (i = 0; i < n; ++i) {
            printf("X%d: %lf\n", i + 1, b[i]);
        }
        printf("----------Comprovació---------\n");
        double *Ax = prodMatVect(aTemp, bTemp, n);
        for (i = 0; i < n; ++i) {
            if(Ax[i] != bTemp[i]){
                is_exact = 1;
                printf("%.20lf(Ax) != %.20lf(b)\n", Ax[i], bTemp[i]);
            }else{
                printf("%.20lf(Ax) = %.20lf(b)\n", Ax[i], bTemp[i]);
            }
        }
        if(is_exact == 0){
            printf("La solució és exacta\n");
        }else{
            printf("La solució és aproximada\n");
            for (i = 0; i < n; ++i) {
                r[i] += (Ax[i] - bTemp[i]);
            }
            printf("\n");
            printf("Vector Residu:\n");
            for (i = 0; i < n; ++i) {
                printf("%.20lf\n", r[i]);
            }
            printf("\n");
        }
        printf("\n");
    } else {
        printf("\033[0;31m");
        printf("El sistema no es pot resoldre\n");
        printf("\033[0;37m");
    }
    free(A);
    free(aTemp);
    free(b);
    free(bTemp);
    free(r);
    return 0;
}