#include <stdlib.h>
#include <stdio.h>
#include "funs_linalg.h"

int main() {
    int i, j, n, option, is_solved;
    FILE *file;
    double **A, *Ax, **aTemp, *b, *bTemp, tol = 0.0001;
    printf("Metode per a resoldre els sistemes:\n");
    printf("[1]-Gauss sense pivotatge:\n");
    printf("[2]-Gauss amb pivotatge:\n");
    scanf("%d", &option);
    if(option == 1) {
        file = fopen("frankG.out","w");
    }else{
        file = fopen("frankGP.out","w");
    }
    for (n = 2; n <= 24; n++) {
        A = (double **) malloc(sizeof(double) * n);
        aTemp = (double **) malloc(sizeof(double) * n);
        b = (double *) malloc(sizeof(double) * n);
        bTemp = (double *) malloc(sizeof(double) * n);
        Ax = (double *) malloc(sizeof(double) * n);

        if (A == NULL || aTemp == NULL || b == NULL || bTemp == NULL || Ax == NULL) {
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
        genMatFrankHessenberg(n, A);
        genMatFrankHessenberg(n, aTemp);
        genVectId(n, b);
        genVectId(n, bTemp);
        prodMatVect(A, bTemp, b, n);
        printf("-------------------------------\n");
        printf("Matriu de Hessenberg per n = %d\n", n);
        printf("-------------------------------\n");
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%24.16e \t", A[i][j]);

            }
            printf("\n");
        }
        printf("--------\n");
        printf("Vector B\n");
        printf("--------\n");
        for (i = 0; i < n; i++) {
            bTemp[i] = b[i];
            printf("%24.16e \t", b[i]);
        }
        printf("\n");
        if (option == 1) {
            is_solved = gauss(n, A, b, tol);
            file = freopen("frankG.out","a",file);
        }else{
            is_solved = gausspiv(n, A, b, tol);
            file = freopen("frankGP.out","a",file);
        }
        if (is_solved == 0) {
            printf("-----------Solucions----------\n");
            for (i = 0; i < n; ++i) {
                printf("X%d: %24.16e\n", i + 1, b[i]);
            }
            prodMatVect(aTemp, b, Ax, n);
            restVectVect(Ax,bTemp,Ax,n);
            fprintf(file,"%d\t %24.16e\n",n,calcNormEucl(n,Ax)/calcNormEucl(n,bTemp));
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
        fclose(file);
    }

    return 0;
}
