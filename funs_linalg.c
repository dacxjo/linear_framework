#include "funs_linalg.h"
#include "stdlib.h"
#include "math.h"
#include "stdio.h"


double prod_esc(int n, double *x, double *y) {
    double result = 0.0;
    int i;
    for (i = 0; i < n; ++i) {
        result += x[i] * y[i];
    }
    return result;
}

void prodMatVect(double **M, double *x, double *destination, int n) {
    int i, j;
    genVectNul(n,destination);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            destination[i] += M[i][j] * x[j];
        }
    }
}


void prodMatMat(double **A, double **B, double **destination, int n, int m) {
    int i, j, k;
    double sum;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            for (k = 0; k < m; k++) {
                sum = sum + A[i][k] * B[k][j];
            }
            destination[i][j] = sum;
            sum = 0;
        }
    }
}

void sumMatMat(double **A, double **B, double **destination, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            destination[i][j] = A[i][j] + B[i][j];
        }
    }
}


void restMatMat(double **A, double **B, double **destination, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            destination[i][j] = A[i][j] - B[i][j];
        }
    }
}
void restVectVect(double *V, double *W,double *destination,int n){
    int i;
    for (i = 0; i < n; i++) {
            destination[i] = V[i] - W[i];
    }
}

double **transposar(double **a, int m, int n) {
    double **temp = (double **) malloc(sizeof(double) * m);
    int i, j;
    for (i = 0; i < n; i++) {
        temp[i] = (double *) malloc(n * sizeof(double));
        if (temp[i] == NULL) {
            exit(EXIT_FAILURE);
        }
    }
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            temp[i][j] = a[j][i];
        }
    }
    return temp;
}

void genMatId(int n, double **M) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                M[i][j] = 1;
            } else {
                M[i][j] = 0;
            }
        }
    }
}

void genMatNul(int n, double **M) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            M[i][j] = 0.0;
        }
    }
}

void genMatFrankHessenberg(int n, double **M) {
    int i,j;
    genMatNul(n, M);
    for ( i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            if (j <= (i - 2)) {
                continue;
            } else if (j == (i - 1)) {
                M[i - 1][j - 1] = n + 1 - i;
            } else if (j >= i) {
                M[i - 1][j - 1] = n + 1 - j;
            }
        }
    }
}

void genVectNul(int n, double *V) {
    int i;
    for (i = 0; i < n; i++) {
        V[i] = 0.0;
    }
}

void genVectId(int n, double *V) {
    int i;
    for (i = 0; i < n; i++) {
        V[i] = 1.0;
    }
}

double calcNormEucl(int n, double *V) {
    return sqrt(prod_esc(n,V,V));
}


int resoltrisup(int n, double **A, double *b, double *x, double tol) {
    int i, j;
    if(A[n - 1][n - 1] == 0.0){
        return 1;
    }
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (i = n - 2; i >= 0; --i) {
        double sum = 0.0;
        for (j = i + 1; j <= n - 1; ++j) {
            sum += (A[i][j] * x[j]);
        }
        if (fabs(A[i][i]) < tol || fabs(A[i][i]) == 0) {
            return 1;
        }
        x[i] += (b[i] - sum) / A[i][i];
    }
    return 0;
}

int resoltriinf(int n, double **A, double *b, double *x, double tol) {
    int i, j;
    if(A[0][0] == 0.0){
        return 1;
    }
    x[0] = b[0] / A[0][0];
    for (i = 1; i < n; i++) {
        double sum = 0.0;
        for (j = 0; j <= (i - 1); j++) {
            sum += (A[i][j] * x[j]);
        }
        if (fabs(A[i][i]) < tol || fabs(A[i][i]) == 0.0) {
            return 1;
        }
        x[i] += (b[i] - sum) / A[i][i];
    }
    return 0;
}

int checktrisup(double **A, int n) {
    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (j < i && A[i][j] != 0) {
                return 0;
            }
        }
    }
    return 1;
}

int checktriinf(double **A, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            if (A[i][j] != 0) {
                return 0;
            }
        }
    }
    return 1;
}

int gauss(int n, double **A, double *b, double tol) {
    int i, j, k;
    double *x, ratio;
    x = (double *) malloc(sizeof(double) * n);
    genVectNul(n, x);
    for (k = 0; k < n - 1; k++) {
        for (i = k + 1; i < n; i++) {
            if (fabs(A[k][k]) == 0.0 || fabs(A[k][k]) < tol) {
                return 1;
            }
            ratio = A[i][k] / A[k][k];
            b[i] = b[i] - ratio * b[k];
            for (j = k; j < n; j++) {
                A[i][j] = A[i][j] - (ratio * A[k][j]);
            }
        }
    }
    if (resoltrisup(n, A, b, x, tol) == 1) {
        return 1;
    }
    for (i = 0; i < n; i++) {
        b[i] = x[i];
    }
    free(x);
    return 0;
}


int gausspiv(int n, double **A, double *b, double tol) {
    int i, j, k;
    double *x, ratio, temp, tempB;
    x = (double *) malloc(sizeof(double) * n);
    genVectNul(n, x);
    for (k = 0; k < n - 1; k++) {
        for (i = k + 1; i < n; i++) {
            if (fabs(A[k][k]) < fabs(A[i][k])) {
                tempB = b[k];
                b[k] = b[i];
                b[i] = tempB;
                for (j = 0; j < n; j++) {
                    temp = A[k][j];
                    A[k][j] = A[i][j];
                    A[i][j] = temp;
                }
            }
        }
        for (i = k + 1; i < n; i++) {
            if (fabs(A[k][k]) == 0.0) {
                return 1;
            }
            ratio = A[i][k] / A[k][k];
            b[i] = b[i] - ratio * b[k];
            for (j = k; j < n; j++) {
                A[i][j] = A[i][j] - (ratio * A[k][j]);
            }
        }
    }
    if (resoltrisup(n, A, b, x, tol) == 1) {
        return 1;
    }
    for (i = 0; i < n; i++) {
        b[i] = x[i];
    }
    free(x);
    return 0;
}


void gaussLU(int n, double **A) {
    int i, j, k;
    double ratio, **temp;
    temp = (double **) malloc(sizeof(double) * n);
    if (temp == NULL) {
        printf("No hi ha suficient memoria\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < n; i++) {
        temp[i] = (double *) malloc(n * sizeof(double));
        if (temp[i] == NULL) {
            printf("No hi ha suficient memoria\n");
            exit(EXIT_FAILURE);
        }
    }
    genMatNul(n,temp);
    for (k = 0; k < n - 1; k++) {
        for (i = k + 1; i < n; i++) {
            if (fabs(A[k][k]) == 0.0) {
                exit(EXIT_FAILURE);
            }
            ratio = A[i][k] / A[k][k];
            for (j = k; j < n; j++) {
                A[i][j] = A[i][j] - (ratio * A[k][j]);
                if (i > j) {
                    temp[i][j] = ratio;
                } else {
                    temp[i][j] = 0.0;
                }
            }
        }
    }
    sumMatMat(A, temp, A, n);
    free(temp);
}

void luDecompose(int n, double **acp, double **L, double **U) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i > j) {
                L[i][j] = acp[i][j];
            } else {
                U[i][j] = acp[i][j];
            }
        }
    }
}

double checkLU(int n, double **a, double **acp) {
    int i, j;
    double **L, **U, **B, maxB;
    L = (double **) malloc(sizeof(double) * n);
    U = (double **) malloc(sizeof(double) * n);
    B = (double **) malloc(sizeof(double) * n);
    if (L == NULL || U == NULL || B == NULL) {
        printf("No hi ha suficient memoria\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < n; i++) {
        L[i] = (double *) malloc(n * sizeof(double));
        U[i] = (double *) malloc(n * sizeof(double));
        B[i] = (double *) malloc(n * sizeof(double));
        if (L[i] == NULL || U[i] == NULL || B[i] == NULL) {
            printf("No hi ha suficient memoria\n");
            exit(EXIT_FAILURE);
        }
    }
    genMatId(n, L);
    genMatNul(n,U);
    genMatNul(n,B);
    luDecompose(n, acp, L, U);
    printf("Matriu L: \n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%e\t", L[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("Matriu U: \n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%e\t", U[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    /** Producte de L per U, acp es el resultat */
    prodMatMat(L, U, acp, n, n);
    /** Resta  de A menys LU,B es el resultat */
    restMatMat(a, acp, B, n);

    printf("Matriu LU: \n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%e\t", acp[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("Matriu B: \n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%e\t", B[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    maxB = fabs(B[0][0]);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (fabs(B[i][j]) > maxB) {
                maxB = fabs(B[i][j]);
            }
        }
    }
    free(L);
    free(U);
    free(B);
    return maxB;
}
