#include <stdlib.h>
#include <stdio.h>
#include "funs_linalg.h"

int main() {
    int i, j, n, option, is_solved, is_exact;
    double **A, *Ax, **aTemp, *b, *bTemp, *r, tol = 0.0001;
    printf("Metod per resoldre els sistemes:\n");
    printf("[1]-Gauss sense pivotatge:\n");
    printf("[2]-Gauss amb pivotatge:\n");
    scanf("%d", &option);
    A = (double **) malloc(sizeof(double) * n);
    aTemp = (double **) malloc(sizeof(double) * n);
    b = (double *) malloc(sizeof(double) * n);
    bTemp = (double *) malloc(sizeof(double) * n);
    r = (double *) malloc(sizeof(double) * n);
    Ax = (double *) malloc(sizeof(double) * n);
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

    for (n = 2; n <= 24; n++) {
        genMatHessenberg(n, A);
        genMatHessenberg(n, aTemp);
        genVectId(n, b);
        genVectId(n, bTemp);
        prodMatVect(A, bTemp, b, n);
        printf("Matriu de Hessenberg per n = %d\n", n);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%lf \t", A[i][j]);

            }
            printf("\n");
        }
        printf("Vector B\n");
        for (i = 0; i < n; i++) {
            bTemp[i] = b[i];
            printf("%lf \t", b[i]);
            printf("\n");
        }
        if (option == 1) {
            is_solved = gauss(n, A, b, tol);
        }else{
            is_solved = gausspiv(n, A, b, tol);
        }
        if (is_solved == 0) {
            printf("El sistema s'ha resolt: Status %d\n", is_solved);
            printf("-----------Matriu de gauss------------\n");
            for (i = 0; i < n; ++i) {
                for (j = 0; j < n; ++j) {
                    printf("%le\t", A[i][j]);
                }
                printf("\n");
            }
            printf("\n");
            printf("-----------Solucions----------\n");
            for (i = 0; i < n; ++i) {
                printf("X%d: %le\n", i + 1, b[i]);
            }
            printf("----------Comprovació---------\n");
            prodMatVect(aTemp, b, Ax, n);
            is_exact = 0;
            for (i = 0; i < n; ++i) {
                if (Ax[i] != bTemp[i]) {
                    is_exact = 1;
                    printf("%.24le(Ax) != %.24le(b)\n", Ax[i], bTemp[i]);
                } else {
                    printf("%.24le(Ax) = %.24le(b)\n", Ax[i], bTemp[i]);
                }
            }
            if (is_exact == 0) {
                printf("La solució és exacta\n");
            } else {
                printf("La solució és aproximada\n");
                for (i = 0; i < n; ++i) {
                    r[i] += (Ax[i] - bTemp[i]);
                }
                printf("\n");
                printf("Vector Residu:\n");
                for (i = 0; i < n; ++i) {
                    printf("%.24le\n", r[i]);
                }
                printf("\n");
            }
            printf("\n");
        } else {
            printf("\033[0;31m");
            printf("El sistema no es pot resoldre\n");
            printf("\033[0;37m");
        }

    }

    free(A);
    free(Ax);
    free(aTemp);
    free(b);
    free(bTemp);
    free(r);
    return 0;
}
