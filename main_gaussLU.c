#include<stdio.h>
#include <stdlib.h>
#include "funs_linalg.h"

int main() {
    double **a,**acp,result;
    int i, j, n,option = 0;
    while (option != 3) {
        printf("---------Factorizacio LU---------\n");
        printf("[1]. Usar matriu de l'exercici(A2)\n");
        printf("[2]. Introduir dades manualment\n");
        printf("[3]. Sortir\n");
        printf("--------------------------------- \n");
        scanf("%d", &option);
        switch (option) {
            case 1: {
                n = 4;
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

                printf("\n");
                gaussLU(n,acp);
                printf("LU Decomposition: \n");
                printf("\n");
                printf("Matriu ACP: \n");
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        printf("%24.16e\t",acp[i][j]);
                    }
                    printf("\n");
                }
                printf("\n");
                printf("Matriu A: \n");
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        printf("%24.16e\t",a[i][j]);
                    }
                    printf("\n");
                }
                printf("\n");
                result = checkLU(n,a,acp);
                printf("Maxim de B:  %24.16e\n",result);
                if(result == 0.0){
                    printf("La factorizacio es correcta!!!");
                    printf("\n");
                }else{
                    printf("La factorizacio es incorrecta...");
                    printf("\n");
                }
                printf("\n");
                free(a);
                free(acp);
                break;
            }
            case 2:
                printf("Ingressi la dimensio de la matriu:\n");
                scanf("%d", &n);
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
                printf("Ingressi els elements de la matriu:\n");
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        scanf("%le", &a[i][j]);
                        acp[i][j] = a[i][j];
                    }
                }
                printf("\n");
                gaussLU(n,acp);
                printf("LU Decomposition: \n");
                printf("\n");
                printf("Matriu ACP: \n");
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        printf("%24.16e\t",acp[i][j]);
                    }
                    printf("\n");
                }
                printf("\n");
                printf("Matriu A: \n");
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        printf("%24.16e\t",a[i][j]);
                    }
                    printf("\n");
                }
                printf("\n");
                result = checkLU(n,a,acp);
                printf("Maxim de B:  %24.16e\n",result);
                if(result == 0.0){
                    printf("La factorizacio es correcta!!!");
                    printf("\n");
                }else{
                    printf("La factorizacio es incorrecta...");
                    printf("\n");
                }
                printf("\n");
                free(a);
                free(acp);
                break;
            default:
                exit(EXIT_SUCCESS);

        }
    }

    return (0);
}

