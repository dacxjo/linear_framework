#include <stdlib.h>
#include <stdio.h>
#include "funs_linalg.h"

int main(){
    int n,is_solved;
    double **A, *b, *x, tol = 0.0;
    printf("Ingressi la dimensió de la matriu:\n");
    scanf("%d", &n);
    A = (double **) malloc(sizeof(double) * n);
    x = (double *) malloc(sizeof(double) * n);
    b = (double *) malloc(sizeof(double) * n);
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
    printf("Ingressi els elements de la matriu:\n");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            scanf("%lf", &A[i][j]);
        }
    }
    printf("Ingressi els elements del vector B:\n");
    for (int i = 0; i < n; ++i) {
        scanf("%lf", &b[i]);
    }
    printf("Matriu ingressada:\n");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%lf\t", A[i][j]);
        }
        printf("|\t %lf",b[i]);
        printf("\n");
    }
    is_solved = gauss(n,A,b,tol);
    printf("Matriu gauss:\n");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%lf\t", A[i][j]);
        }
        printf("|\t %lf",b[i]);
        printf("\n");
    }
    if(is_solved == 0){
        resoltrisup(n,A,b,x,tol);
        printf("\n");
        printf("-----------Solucions----------\n");
        for (int i = 0; i < n; ++i) {
            printf("X%d: %lf\n",i+1, x[i]);
        }
        printf("-------------------------------\n");
        double *Ax = prodMatVect(A,x,n);
        printf("Comprovació:\n");
        for (int i = 0; i < n; ++i) {
            printf("%lf(Ax) = %lf(b)\n",Ax[i],b[i]);
        }
    }


    return 0;
}