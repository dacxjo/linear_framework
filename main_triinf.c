#include <stdio.h>
#include <stdlib.h>
#include "funs_linalg.h"

int main() {
    int i, j, is_solved, is_exact = 0, option = 0;
    double **A, **matriuEx, *vecEx, *Ax, *b, *x, *r, tol = 0.0;
    int n;
    while (option != 3) {
        printf("---------Triangular inferior---------\n");
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
                matriuEx = (double **) malloc(sizeof(double) * n);
                vecEx = (double *) malloc(sizeof(double) * n);
                A = (double **) malloc(sizeof(double) * n);
                x = (double *) malloc(sizeof(double) * n);
                b = (double *) malloc(sizeof(double) * n);
                r = (double *) malloc(sizeof(double) * n);

                if (A == NULL || x == NULL || b == NULL || r == NULL || matriuEx == NULL || vecEx == NULL) {
                    printf("No hi ha suficient memòria\n");
                    exit(EXIT_FAILURE);
                }
                for (i = 0; i < n; i++) {
                    A[i] = (double *) malloc(n * sizeof(double));
                    matriuEx[i] = (double *) malloc(n * sizeof(double));
                    if (A[i] == NULL || matriuEx[i] == NULL) {
                        printf("No hi ha suficient memòria\n");
                        exit(EXIT_FAILURE);
                    }
                }

                matriuEx[0][0] = 1.0234;
                matriuEx[0][1] = 2.0981;
                matriuEx[0][2] = 9.9871;
                matriuEx[0][3] = 1.1;
                matriuEx[1][0] = 0;
                matriuEx[1][1] = -6.9876;
                matriuEx[1][2] = 2.2222;
                matriuEx[1][3] = 0.3333;
                matriuEx[2][0] = 0;
                matriuEx[2][1] = 0;
                matriuEx[2][2] = -1.9870;
                matriuEx[2][3] = 20.121;
                matriuEx[3][0] = 0;
                matriuEx[3][1] = 0;
                matriuEx[3][2] = 0;
                matriuEx[3][3] = 1.1234;

                vecEx[0] = 0;
                vecEx[1] = 1;
                vecEx[2] = 0;
                vecEx[3] = 1;
                printf("---------------Matriu----------------\n");
                printf("\n");
                for (i = 0; i < n; ++i) {
                    for (j = 0; j < n; ++j) {
                        printf("%le\t", matriuEx[i][j]);
                        A[i][j] = matriuEx[i][j];
                    }
                    printf("|\t %le", vecEx[i]);
                    b[i] = vecEx[i];
                    printf("\n");
                }
                A = transposar(A, n, n);
                printf("----------Matriu transposada----------\n");
                printf("\n");
                for (i = 0; i < n; ++i) {
                    for (j = 0; j < n; ++j) {
                        printf("%le\t", A[i][j]);

                    }
                    printf("|\t %le", vecEx[i]);
                    printf("\n");
                }
                printf("\n");
                is_solved = resoltriinf(n, A, b, x, tol);
                if (is_solved == 0) {
                    printf("El sistema s'ha resolt: Status %d\n", is_solved);
                    printf("\n");
                    printf("-----------Solucions----------\n");
                    for (i = 0; i < n; ++i) {
                        printf("X%d: %le\n", i + 1, x[i]);
                    }
                    printf("\n");
                    printf("----------Comprovació---------\n");
                    Ax = prodMatVect(A, x, n);
                    for (i = 0; i < n; ++i) {
                        if (Ax[i] != b[i]) {
                            is_exact = 1;
                            printf("%.20lf(Ax) != %.20lf(b)\n", Ax[i], b[i]);
                        } else {
                            printf("%.20lf(Ax) = %.20lf(b)\n", Ax[i], b[i]);
                        }
                    }
                    if (is_exact == 0) {
                        printf("La solució és exacta\n");
                    } else {
                        printf("La solució és aproximada\n");
                        for (i = 0; i < n; ++i) {
                            r[i] += (Ax[i] - b[i]);
                        }
                        printf("\n");
                        printf("Vector Residu:\n");
                        for (i = 0; i < n; ++i) {
                            printf("%.20lf\n", r[i]);
                        }
                        printf("\n");
                    }
                    free(Ax);
                } else {
                    printf("\033[0;31m");
                    printf("El sistema no es pot resoldre\n");
                    printf("\033[0;37m");
                }
                free(A);
                free(b);
                free(x);
                free(r);
                free(matriuEx);
                free(vecEx);
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
                for (i = 0; i < n; i++) {
                    A[i] = (double *) malloc(n * sizeof(double));
                    if (A[i] == NULL) {
                        printf("No hi ha suficient memòria\n");
                        exit(EXIT_FAILURE);
                    }
                }
                printf("Ingressi els elements de la matriu:\n");
                for (i = 0; i < n; ++i) {
                    for (j = 0; j < n; ++j) {
                        scanf("%le", &A[i][j]);
                    }
                }
                if (checktriinf(A, n)) {
                    printf("La matriu és triangular inferior\n");
                } else {
                    printf("\033[0;31m");
                    printf("La matriu no és triangular inferior\n");
                    exit(EXIT_FAILURE);
                }
                printf("Ingressi els elements del vector B:\n");
                for (i = 0; i < n; ++i) {
                    scanf("%le", &b[i]);
                }
                printf("Ingressi la tolerància acceptada:\n");
                scanf("%le", &tol);
                printf("Matriu ingressada:\n");
                for (i = 0; i < n; ++i) {
                    for (j = 0; j < n; ++j) {
                        printf("%le\t", A[i][j]);
                    }
                    printf("|\t %le", b[i]);
                    printf("\n");
                }
                is_solved = resoltriinf(n, A, b, x, tol);
                if (is_solved == 0) {
                    printf("El sistema s'ha resolt: Status %d\n", is_solved);
                    printf("\n");
                    printf("-----------Solucions----------\n");
                    for (i = 0; i < n; ++i) {
                        printf("X%d: %le\n", i + 1, x[i]);
                    }
                    printf("-------------------------------\n");
                    Ax = prodMatVect(A, x, n);
                    for (i = 0; i < n; ++i) {
                        if (Ax[i] != b[i]) {
                            is_exact = 1;
                            printf("%.20lf(Ax) != %.20lf(b)\n", Ax[i], b[i]);
                        } else {
                            printf("%.20lf(Ax) = %.20lf(b)\n", Ax[i], b[i]);
                        }
                    }
                    if (is_exact == 0) {
                        printf("La solució és exacta\n");
                    } else {
                        printf("La solució és aproximada\n");
                        for (i = 0; i < n; ++i) {
                            r[i] += (Ax[i] - b[i]);
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
