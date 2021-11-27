#include <stdio.h>
#include <stdlib.h>
#include "funs_linalg.h"

int main() {
    int i, j, is_solved, is_exact, option = 0;
    double **A, *Ax, *b, *x, *r, tol = 0.0;
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
                printf("Ingressi la tolerància acceptada:\n");
                scanf("%le", &tol);
                A = (double **) malloc(sizeof(double) * n);
                Ax = (double *) malloc(sizeof(double) * n);
                x = (double *) malloc(sizeof(double) * n);
                b = (double *) malloc(sizeof(double) * n);
                r = (double *) malloc(sizeof(double) * n);
                if (A == NULL || x == NULL || b == NULL || r == NULL || Ax == NULL) {
                    printf("No hi ha suficient memòria\n");
                    exit(EXIT_FAILURE);
                }
                for (i = 0; i < n; i++) {
                    A[i] = (double *) malloc(n * sizeof(double));
                    if (A[i] == NULL) {
                        printf("No hi ha suficient memòria\n");
                        exit(EXIT_FAILURE);
                    }
                }

                A[0][0] = 1.0234;
                A[0][1] = 2.0981;
                A[0][2] = 9.9871;
                A[0][3] = 1.1;
                A[1][0] = 0;
                A[1][1] = -6.9876;
                A[1][2] = 2.2222;
                A[1][3] = 0.3333;
                A[2][0] = 0;
                A[2][1] = 0;
                A[2][2] = -1.9870;
                A[2][3] = 20.121;
                A[3][0] = 0;
                A[3][1] = 0;
                A[3][2] = 0;
                A[3][3] = 1.1234;

                b[0] = 0;
                b[1] = 1;
                b[2] = 0;
                b[3] = 1;

                printf("---------------Matriu----------------\n");
                printf("\n");
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        printf("%.24le\t", A[i][j]);
                    }
                    printf("|\t %.24le", b[i]);
                    printf("\n");
                }
                printf("\n");
                genVectNul(n, x);
                is_solved = resoltrisup(n, A, b, x, tol);
                if (is_solved == 0) {
                    genVectNul(n, r);
                    genVectNul(n, Ax);
                    printf("El sistema s'ha resolt: Status %d\n", is_solved);
                    printf("\n");
                    printf("-----------Solucions----------\n");
                    for (i = 0; i < n; i++) {
                        printf("X%d: %.24le\n", i + 1, x[i]);
                    }
                    printf("\n");
                    printf("----------Comprovació---------\n");
                    prodMatVect(A, x, Ax, n);
                    for (i = 0; i < n; i++) {
                        if (Ax[i] != b[i]) {
                            is_exact = 1;
                            printf("%.24le(Ax) != %.24le(b)\n", Ax[i], b[i]);
                        } else {
                            printf("%.24le(Ax) = %.24le(b)\n", Ax[i], b[i]);
                        }
                    }
                    if (is_exact == 0) {
                        printf("La solució és exacta\n");
                    } else {
                        printf("La solució és aproximada\n");
                        for (i = 0; i < n; i++) {
                            r[i] += (Ax[i] - b[i]);
                        }
                        printf("\n");
                        printf("Vector Residu:\n");
                        for (i = 0; i < n; i++) {
                            printf("%.24le\n", r[i]);
                        }
                        printf("\n");
                        printf("Residu: %.24le\n", calcNormEucl(n, r));
                        printf("\n");
                    }
                } else {
                    printf("\033[0;31m");
                    printf("El sistema no es pot resoldre\n");
                    printf("\033[0;37m");
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
                Ax = (double *) malloc(sizeof(double) * n);
                if (A == NULL || b == NULL || x == NULL || Ax == NULL) {
                    printf("No hi ha suficient memòria\n");
                    exit(EXIT_FAILURE);
                }
                for (i = 0; i < n; i++) {
                    A[i] = (double *) malloc(n * sizeof(double));
                    if (A[i] == NULL) {
                        printf("No hi ha suficient memòria\n");
                        exit(EXIT_FAILURE);
                    }
                }
                printf("Ingressi els elements de la matriu:\n");
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        scanf("%le", &A[i][j]);
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
                for (i = 0; i < n; i++) {
                    scanf("%le", &b[i]);
                }
                printf("Ingressi la tolerància acceptada:\n");
                scanf("%le", &tol);
                printf("\n");
                printf("Matriu ingressada:\n");
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        printf("%.24le\t", A[i][j]);
                    }
                    printf("|\t %.24le", b[i]);
                    printf("\n");
                }
                genVectNul(n, x);
                is_solved = resoltrisup(n, A, b, x, tol);
                if (is_solved == 0) {
                    genVectNul(n, r);
                    genVectNul(n, Ax);
                    printf("\n");
                    printf("El sistema s'ha resolt: Status %d\n", is_solved);
                    printf("\n");
                    printf("-----------Solucions----------\n");
                    for (i = 0; i < n; i++) {
                        printf("X%d: %.24le\n", i + 1, x[i]);
                    }
                    printf("\n");
                    printf("----------Comprovació---------\n");
                    prodMatVect(A, x, Ax, n);
                    is_exact = 0;
                    for (i = 0; i < n; i++) {
                        if (Ax[i] != b[i]) {
                            is_exact = 1;
                            printf("%.24le(Ax) != %.24le (b)\n", Ax[i], b[i]);
                        } else {
                            printf("%.24le(Ax) = %.24le (b)\n", Ax[i], b[i]);
                        }
                    }
                    printf("\n");
                    if (is_exact == 0) {
                        printf("La solució és exacta\n");
                        printf("\n");
                    } else {
                        printf("La solució és aproximada\n");
                        for (i = 0; i < n; i++) {
                            r[i] += (Ax[i] - b[i]);
                        }
                        printf("\n");
                        printf("Vector Residu:\n");
                        for (i = 0; i < n; i++) {
                            printf("%.24le\n", r[i]);
                        }
                        printf("\n");
                    }
                    printf("\n");
                    printf("Residu: %.24le\n", calcNormEucl(n, r));
                    printf("\n");
                } else {
                    printf("\033[0;31m");
                    printf("El sistema no es pot resoldre\n");
                    printf("\033[0;37m");
                }
                free(A);
                free(Ax);
                free(b);
                free(x);
                free(r);
                break;
            default:
                exit(EXIT_SUCCESS);

        }
    }
    return 0;
}
