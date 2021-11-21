#include<stdio.h>
#include<math.h>
#include <stdlib.h>

int main() {
    double **a, **l,term;
    int i, j, n = 4,k;
    a = (double **) malloc(sizeof(double) * n);
    l = (double **) malloc(sizeof(double) * n);
    for (i = 0; i < n; i++) {
        a[i] = (double *) malloc(n * sizeof(double));
        l[i] = (double *) malloc(n * sizeof(double));
        if (a[i] == NULL || l[i] == NULL) {
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

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                l[i][j] = 1;
            } else {
                l[i][j] = 0;
            }
        }
    }
    printf("\nL BEFORE FUNCTION:\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%lf\t", l[i][j]);
        }
        printf("\n");
    }
    for (k = 0; k <= n - 1; k++) {
        for (i = k + 1; i < n; i++) {
            if (fabs(a[k][k]) == 0) {
                exit(EXIT_FAILURE);
            }
            term = a[i][k] / a[k][k];
            for (j = k; j <= n; j++) {
                a[i][j] = a[i][j] - (term * a[k][j]);
                if(i > j){
                    l[i][j] = term;
                }
            }
        }
    }

    printf("\nL AFTER FUNCTION:\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%lf\t", l[i][j]);
        }
        printf("\n");
    }
    printf("\nU:\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%lf\t", a[i][j]);
        }
        printf("\n");
    }
    free(a);
    free(l);
    return (0);
}

