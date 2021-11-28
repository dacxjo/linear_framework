#include <stdlib.h>
#include <stdio.h>
#include "funs_linalg.h"

int main() {
    int i, j, n, is_solved,is_exact,option = 0;
    double **A,*Ax, **aTemp, *b, *bTemp,*r, tol = 0.0;
    while (option != 3) {
        printf("---------Metode de Gauss amb pivotatge---------\n");
        printf("[1]. Usar matriu de l'exercici(A2)\n");
        printf("[2]. Introduir dades manualment\n");
        printf("[3]. Sortir\n");
        printf("-----------------------------------------------\n");
        scanf("%d", &option);
        switch (option) {
            case 1: {
                n = 4;
                printf("Ingressi la tolerancia acceptada:\n");
                scanf("%le", &tol);
                A = (double **) malloc(sizeof(double) * n);
                aTemp = (double **) malloc(sizeof(double) * n);
                b = (double *) malloc(sizeof(double) * n);
                bTemp = (double *) malloc(sizeof(double) * n);
                r = (double *) malloc(sizeof(double) * n);
                Ax = (double *) malloc(sizeof(double) * n);
                if (A == NULL || aTemp == NULL || b == NULL || r == NULL || bTemp == NULL) {
                    printf("No hi ha suficient memoria\n");
                    exit(EXIT_FAILURE);
                }
                for (i = 0; i < n; i++) {
                    A[i] = (double *) malloc(n * sizeof(double));
                    aTemp[i] = (double *) malloc(n * sizeof(double));
                    if (A[i] == NULL) {
                        printf("No hi ha suficient memoria\n");
                        exit(EXIT_FAILURE);
                    }
                    if (aTemp[i] == NULL) {
                        printf("No hi ha suficient memoria\n");
                        exit(EXIT_FAILURE);
                    }
                }

                A[0][0] = 1;
                A[0][1] = 0.5;
                A[0][2] = 0.333333;
                A[0][3] = 0.25;
                A[1][0] = 0.5;
                A[1][1] = 0.333333;
                A[1][2] = 0.25;
                A[1][3] = 0.2;
                A[2][0] =  0.333333;
                A[2][1] = 0.25;
                A[2][2] = 0.2;
                A[2][3] = 0.166667;
                A[3][0] = 0.25;
                A[3][1] = 0.2;
                A[2][3] = 0.166667;
                A[3][3] = 0.142857;

                b[0] = 2.083333;
                b[1] = 1.283333;
                b[2] = 0.95;
                b[3] = 0.759524;

                printf("---------------Matriu ingresada----------------\n");
                printf("\n");
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        printf("%24.16e\t", A[i][j]);
                        aTemp[i][j] = A[i][j];
                    }
                    printf("|\t %24.16e", b[i]);
                    bTemp[i] = b[i];
                    printf("\n");
                }
                printf("\n");
                is_solved = gausspiv(n, A, b, tol);
                if (is_solved == 0) {
                    printf("El sistema s'ha resolt: Status %d\n", is_solved);
                    printf("\n");
                    printf("-----------Matriu de gauss------------\n");
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            printf("%24.16e\t", A[i][j]);
                        }
                        printf("\n");
                    }
                    printf("\n");
                    printf("-----------Solucions----------\n");
                    for (i = 0; i < n; i++) {
                        printf("X%d: %24.16e\n", i + 1, b[i]);
                    }
                    printf("\n");
                    printf("----------Comprovacio---------\n");
                    prodMatVect(aTemp, b,Ax, n);
                    is_exact = 0;
                    for (i = 0; i < n; i++) {
                        if(Ax[i] != bTemp[i]){
                            is_exact = 1;
                            printf("%24.16e(Ax) != %24.16e(b)\n", Ax[i], bTemp[i]);
                        }else{
                            printf("%24.16e(Ax) = %24.16e(b)\n", Ax[i], bTemp[i]);
                        }
                    }
                    if(is_exact == 0){
                        printf("\n");
                        printf("La solucio es exacta\n");
                        printf("\n");
                    }else{
                        printf("\n");
                        printf("La solucio es aproximada\n");
                        for (i = 0; i < n; i++) {
                            r[i] += (Ax[i] - bTemp[i]);
                        }
                        printf("\n");
                        printf("Vector Residu:\n");
                        for (i = 0; i < n; i++) {
                            printf("%24.16e\n", r[i]);
                        }
                        printf("\n");
                        printf("Residu: %24.16e\n", calcNormEucl(n, r));
                        printf("\n");
                    }
                    printf("\n");
                } else {
                    printf("\033[0;31m");
                    printf("El sistema no es pot resoldre\n");
                    printf("\033[0;37m");
                }
                free(A);
                free(Ax);
                free(aTemp);
                free(b);
                free(bTemp);
                free(r);
                break;
            }
            case 2:
                printf("Ingressi la dimensio de la matriu:\n");
                scanf("%d", &n);
                A = (double **) malloc(sizeof(double) * n);
                aTemp = (double **) malloc(sizeof(double) * n);
                b = (double *) malloc(sizeof(double) * n);
                bTemp = (double *) malloc(sizeof(double) * n);
                r = (double *) malloc(sizeof(double) * n);
                Ax = (double *) malloc(sizeof(double) * n);
                if (A == NULL || aTemp == NULL || b == NULL || r == NULL || bTemp == NULL) {
                    printf("No hi ha suficient memoria\n");
                    exit(EXIT_FAILURE);
                }
                for (i = 0; i < n; i++) {
                    A[i] = (double *) malloc(n * sizeof(double));
                    aTemp[i] = (double *) malloc(n * sizeof(double));
                    if (A[i] == NULL) {
                        printf("No hi ha suficient memoria\n");
                        exit(EXIT_FAILURE);
                    }
                    if (aTemp[i] == NULL) {
                        printf("No hi ha suficient memoria\n");
                        exit(EXIT_FAILURE);
                    }
                }
                printf("Ingressi els elements de la matriu:\n");
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        scanf("%le", &A[i][j]);
                        aTemp[i][j] = A[i][j];
                    }
                }
                printf("Ingressi els elements del vector B:\n");
                for (i = 0; i < n; i++) {
                    scanf("%le", &b[i]);
                    bTemp[i] = b[i];
                }
                printf("Ingressi la tolerancia acceptada:\n");
                scanf("%le", &tol);
                printf("---------------Matriu ingresada----------------\n");
                printf("\n");
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        printf("%24.16e\t", A[i][j]);
                    }
                    printf("|\t %24.16e", b[i]);
                    printf("\n");
                }
                printf("\n");
                is_solved = gausspiv(n, A, b, tol);
                if (is_solved == 0) {
                    printf("El sistema s'ha resolt: Status %d\n", is_solved);
                    printf("\n");
                    printf("-----------Matriu de gauss------------\n");
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            printf("%24.16e\t", A[i][j]);
                        }
                        printf("\n");
                    }
                    printf("\n");
                    printf("-----------Solucions----------\n");
                    for (i = 0; i < n; i++) {
                        printf("X%d: %24.16e\n", i + 1, b[i]);
                    }
                    printf("\n");
                    printf("----------Comprovacio---------\n");
                    prodMatVect(aTemp, b,Ax, n);
                    is_exact = 0;
                    for (i = 0; i < n; i++) {
                        if(Ax[i] != bTemp[i]){
                            is_exact = 1;
                            printf("%24.16e(Ax) != %24.16e(b)\n", Ax[i], bTemp[i]);
                        }else{
                            printf("%24.16e(Ax) = %24.16e(b)\n", Ax[i], bTemp[i]);
                        }
                    }
                    if(is_exact == 0){
                        printf("\n");
                        printf("La solucio es exacta\n");
                    }else{
                        printf("\n");
                        printf("La solucio es aproximada\n");
                        for (i = 0; i < n; i++) {
                            r[i] += (Ax[i] - bTemp[i]);
                        }
                        printf("\n");
                        printf("Vector Residu:\n");
                        for (i = 0; i < n; i++) {
                            printf("%24.16e\n", r[i]);
                        }
                        printf("\n");
                        printf("Residu: %24.16e\n", calcNormEucl(n, r));
                        printf("\n");
                    }
                    printf("\n");
                } else {
                    printf("\033[0;31m");
                    printf("El sistema no es pot resoldre\n");
                    printf("\033[0;37m");
                }
                free(A);
                free(Ax);
                free(aTemp);
                free(b);
                free(bTemp);
                free(r);
                break;
            default:
                exit(EXIT_SUCCESS);
        }
    }
    return 0;
}
