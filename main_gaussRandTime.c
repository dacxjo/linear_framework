#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "funs_linalg.h"

int main() {
    int i, j, k, n, is_resolved,option;
    double **A, *b, tol = 0.0001, temps;
    FILE *file;
    clock_t start, end;
    printf("Metod per resoldre els sistemes:\n");
    printf("[1]-Gauss sense pivotatge:\n");
    printf("[2]-Gauss amb pivotatge:\n");
    scanf("%d", &option);

    if(option == 1) {
        file = fopen("../gaussTemps.out","w");
    }else{
        file = fopen("../gaussPivTemps.out","w");
    }


    for (n = 2; n < 200; n++) {
        start = clock();
        A = (double **) malloc(sizeof(double) * n);
        b = (double *) malloc(sizeof(double) * n);
        for (i = 0; i < n; i++) {
            A[i] = (double *) malloc(n * sizeof(double));
            if (A[i] == NULL) {
                printf("No hi ha suficient memòria\n");
                exit(EXIT_FAILURE);
            }
        }
        srand(time(NULL));
        for (k = 0; k < 1000; k++) {
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    A[i][j] = rand() / (double) RAND_MAX;
                    b[i] = rand() / (double) RAND_MAX;
                }
            }
            if(option == 1){
                is_resolved = gauss(n,A,b,tol);
            }else{
                is_resolved = gausspiv(n,A,b,tol);
            }
        }
        free(A);
        free(b);
        end = clock();
        temps = (double) (end - start) / CLOCKS_PER_SEC;
        printf("1000 sistemas %dx%d resueltos en %.24lf segundos\n",n,n, temps);
        if(option == 1){
            file = freopen("../gaussTemps.out","a",file);
        }else{
            file = freopen("../gaussPivTemps.out","a",file);
        }

        fprintf(file,"%d\t %.24lf\n",n,temps);
    }
    fclose(file);
    return 0;
}