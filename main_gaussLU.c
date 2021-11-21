#include<stdio.h>
#include <stdlib.h>
#include "funs_linalg.h"

int main() {
    double **a,**acp,result;
    int i, j, n = 4;
    a = (double **) malloc(sizeof(double) * n);
    acp = (double **) malloc(sizeof(double) * n);
    if(a == NULL || acp == NULL){
        printf("No hi ha suficient memoria\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < n; i++) {
        a[i] = (double *) malloc(n * sizeof(double));
        acp[i] = (double *) malloc(n * sizeof(double));
        if (a[i] == NULL || acp[i] == NULL) {
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


/*
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
*/
    printf("\n");
    gaussLU(n,acp);
    printf("LU Decomposition: \n");
    printf("\n");
    printf("Matriu ACP: \n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%e\t",acp[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("Matriu A: \n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%e\t",a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    result = checkLU(n,a,acp);
    printf("Maxim de B:  %e\n",result);
    if(result == 0.0){
        printf("La factorizacio es correcta!!!");
    }else{
        printf("La factorizacio es incorrecta...");
    }
    free(a);
    free(acp);
    return (0);
}

