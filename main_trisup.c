#include <stdio.h>
#include <stdlib.h>
#include "funs_linalg.h"


int main() {
    int is_solved,is_exact = 0, option = 0;
    double **A, *b, *x, *r, tol = 0.0;
    int n;
    while (option != 3) {
        printf("---------Triangular superior---------\n");
        printf("[1]. Usar matriu de l'exercici\n");
        printf("[2]. Introduir dades manualment\n");
        printf("[3]. Sortir\n");
        printf("-------------------------------------\n");
        scanf("%d", &option);
        switch (option) {
            case 1: {
                n = 4;
                double matriuEx[4][4] = {
                        {1.0234, 2.0981,  9.9871,  1.1},
                        {0,      -6.9876, 2.2222,  0.3333},
                        {0,      0,       -1.9870, 20.121},
                        {0,      0,       0,       1.1234}
                };
                double vecEx[4] = {0, 1, 0, 1};
                A = (double **) malloc(sizeof(double) * n);
                x = (double *) malloc(sizeof(double) * n);
                b = (double *) malloc(sizeof(double) * n);
                r = (double *) malloc(sizeof(double) * n);
                if (A == NULL || x == NULL || b == NULL) {
                    printf("No hi ha suficient memòria\n");
                    exit(EXIT_FAILURE);
                }
                for (int i = 0; i < n; i++) {
                    A[i] = (double *) malloc(n * sizeof(double));
                    if (A[i] == NULL) {
                        printf("No hi ha suficient memòria\n");
                        exit(EXIT_FAILURE);
                    }
                }
                printf("---------------Matriu----------------\n");
                printf("\n");
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        printf("%lf\t", matriuEx[i][j]);
                        A[i][j] = matriuEx[i][j];
                    }
                    printf("|\t %lf", vecEx[i]);
                    b[i] = vecEx[i];
                    printf("\n");
                }
                printf("\n");
                is_solved = resoltrisup(n, A, b, x, tol);
                if (is_solved == 0) {
                    printf("El sistema s'ha resolt: Status %d\n", is_solved);
                    printf("\n");
                    printf("-----------Solucions----------\n");
                    for (int i = 0; i < n; ++i) {
                        printf("X%d: %lf\n", i + 1, x[i]);
                    }
                    printf("\n");
                }
                printf("----------Comprovació---------\n");
                double *Ax = prodMatVect(A, x, n);
                for (int i = 0; i < n; ++i) {
                    if(Ax[i] != b[i]){
                        is_exact = 1;
                        printf("%.20lf(Ax) != %.20lf(b)\n", Ax[i], b[i]);
                    }else{
                        printf("%.20lf(Ax) = %.20lf(b)\n", Ax[i], b[i]);
                    }
                }
                if(is_exact == 0){
                    printf("La solució és exacta\n");
                }else{
                    printf("La solució és aproximada\n");
                    for (int i = 0; i < n; ++i) {
                        r[i] += (Ax[i] - b[i]);
                    }
                    printf("\n");
                    printf("Vector Residu:\n");
                    for (int i = 0; i < n; ++i) {
                        printf("%.20lf\n", r[i]);
                    }
                    printf("\n");
                }
                free(A);
                free(b);
                free(x);
                free(r);
                free(Ax);
                break;
            }
            case 2:
                printf("Ingressi la dimensió de la matriu:\n");
                scanf("%d", &n);
                A = (double **) malloc(sizeof(double) * n);
                b = (double *) malloc(sizeof(double) * n);
                x = (double *) malloc(sizeof(double) * n);
                r = (double *) malloc(sizeof(double) * n);
                if (A == NULL || b == NULL || x == NULL) {
                    printf("No hi ha suficient memòria\n");
                    exit(EXIT_FAILURE);
                }
                for (int i = 0; i < n; i++) {
                    A[i] = (double *) malloc(n * sizeof(double));
                    if (A[i] == NULL) {
                        printf("No hi ha suficient memòria\n");
                        exit(EXIT_FAILURE);
                    }
                }
                printf("Ingressi els elements de la matriu:\n");
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        scanf("%lf", &A[i][j]);
                    }
                }
                if (checktrisup(A, n)) {
                    printf("La matriu és triangular superior\n");
                } else {
                    printf("\033[0;31m");
                    printf("La matriu no és triangular superior\n");
                    exit(2);
                }
                printf("Ingressi els elements del vector B:\n");
                for (int i = 0; i < n; ++i) {
                    scanf("%lf", &b[i]);
                }
                printf("Ingressi la tolerància acceptada:\n");
                scanf("%lf", &tol);
                printf("Matriu ingressada:\n");
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        printf("%lf\t", A[i][j]);
                    }
                    printf("|\t %lf", b[i]);
                    printf("\n");
                }
                is_solved = resoltrisup(n, A, b, x, tol);
                if (is_solved == 0) {
                    printf("El sistema s'ha resolt: Status %d\n", is_solved);
                    printf("\n");
                    printf("-----------Solucions----------\n");
                    for (int i = 0; i < n; ++i) {
                        printf("X%d: %lf\n", i + 1, x[i]);
                    }
                    printf("----------Comprovació---------\n");
                    double *Ax = prodMatVect(A, x, n);
                    for (int i = 0; i < n; ++i) {
                        if(Ax[i] != b[i]){
                            is_exact = 1;
                            printf("%.20lf(Ax) != %.20lf(b)\n", Ax[i], b[i]);
                        }else{
                            printf("%.20lf(Ax) = %.20lf(b)\n", Ax[i], b[i]);
                        }
                    }
                    if(is_exact == 0){
                        printf("La solució és exacta\n");
                    }else{
                        printf("La solució és aproximada\n");
                        for (int i = 0; i < n; ++i) {
                            r[i] += (Ax[i] - b[i]);
                        }
                        printf("\n");
                        printf("Vector Residu:\n");
                        for (int i = 0; i < n; ++i) {
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
                free(b);
                free(x);
                free(r);
                break;

        }
    }
}
