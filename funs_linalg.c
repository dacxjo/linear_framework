#include "funs_linalg.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"


int resoltrisup(int n, double **A, double *b, double *x, double tol) {
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j <= n - 1; ++j) {
            sum += (A[i][j] * x[j]);
        }
        x[i] += (b[i] - sum) / A[i][i];
    }
    return 0;
}

int resoltriinf(int n, double **A, double *b, double *x, double tol) {
    x[0] = b[0] / A[0][0];
    for (int i = 1; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j <= (i - 1); ++j) {
            sum += (A[i][j] * x[j]);
        }
        x[i] += (b[i] - sum) / A[i][i];
    }
    return 0;
}

_Bool checktrisup(double **A, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j < i && A[i][j] != 0) {
                return 0;
            }
        }
    }
    return 1;
}

_Bool checktriinf(double **A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (A[i][j] != 0) {
                return 0;
            }
        }
    }
    return 1;
}

double prod_esc(int n, double *x, double *y) {
    double result = 0.0;
    for (int i = 0; i < n; ++i) {
        result += x[i] * y[i];
    }
    return result;
}

void divMat(double **A, double x, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; ++j) {
            A[i][j] = A[i][j] / x;
        }
    }
}

double **prodMatMat(double **A, double **B, int n, int m) {
    int i, j, k;
    double sum;
    double **O = (double **) malloc(sizeof(double) * n);
    for (i = 0; i < n; i++) {
        O[i] = (double *) malloc(n * sizeof(double));
        if (O[i] == NULL) {
            exit(EXIT_FAILURE);
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            for (k = 0; k < m; k++) {
                sum = sum + A[i][k] * B[k][j];
            }
            O[i][j] = sum;
            sum = 0;
        }
    }
    return O;
}

double *prodMatVect(double **M, double *x, int n) {
    double *v = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        v[i] = 0;
        for (int j = 0; j < n; j++) {
            v[i] += M[i][j] * x[j];
        }
    }
    return v;
}

double **transposar(double **a, int m, int n) {
    double **temp = (double **) malloc(sizeof(double) * m);
    for (int i = 0; i < n; i++) {
        temp[i] = (double *) malloc(n * sizeof(double));
        if (temp[i] == NULL) {
            exit(EXIT_FAILURE);
        }
    }
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            temp[i][j] = a[j][i];
        }
    }
    return temp;
}

double determinant(double **A, int k) {
    double s = 1, det = 0;
    int i, j, m, n, c;
    double **b = (double **) malloc(sizeof(double) * k);
    for (i = 0; i < k; i++) {
        b[i] = (double *) malloc(k * sizeof(double));
        if (b[i] == NULL) {
            exit(EXIT_FAILURE);
        }
    }
    if (k == 1) {
        return (A[0][0]);
    } else {
        det = 0;
        for (c = 0; c < k; c++) {
            m = 0;
            n = 0;
            for (i = 0; i < k; i++) {
                for (j = 0; j < k; j++) {
                    b[i][j] = 0;
                    if (i != 0 && j != c) {
                        b[m][n] = A[i][j];
                        if (n < (k - 2))
                            n++;
                        else {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (A[0][c] * determinant(b, k - 1));
            s = -1 * s;
        }
    }
    return (det);
}


double **matCofact(double **A, int f) {
    int p, q, m, n, i, j;
    double **b = (double **) malloc(sizeof(double) * f);
    double **fac = (double **) malloc(sizeof(double) * f);
    double **adj = (double **) malloc(sizeof(double) * f);
    for (i = 0; i < f; i++) {
        b[i] = (double *) malloc(f * sizeof(double));
        fac[i] = (double *) malloc(f * sizeof(double));
        adj[i] = (double *) malloc(f * sizeof(double));
        if (b[i] == NULL || fac[i] == NULL || adj[i] == NULL) {
            exit(EXIT_FAILURE);
        }
    }
    for (q = 0; q < f; q++) {
        for (p = 0; p < f; p++) {
            m = 0;
            n = 0;
            for (i = 0; i < f; i++) {
                for (j = 0; j < f; j++) {
                    if (i != q && j != p) {
                        b[m][n] = A[i][j];
                        if (n < (f - 2))
                            n++;
                        else {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
        }
    }
    for (i = 0; i < f; i++) {
        for (j = 0; j < f; j++) {
            adj[i][j] = fac[j][i];
        }
    }

    adj = transposar(adj, f, f);
    return adj;
}

int gauss(int n, double **A, double *b, double tol) {
    int i, j, k, is_solved;
    double *temp = (double *) malloc(sizeof(double) * n);
    for (i = 0; i < n; i++) {
        if (A[i][i] == 0) {
            fprintf(stderr, "Divisió per 0!\n");
            return 1;
        }
        for (j = i + 1; j < n; j++) {
            double term = A[j][i] / A[i][i];
            if (b != NULL) {
                b[j] = b[j] - term * b[i];
            }
            for (k = 0; k <= n; k++) {
                A[j][k] = A[j][k] - term * A[i][k];
            }
        }
    }
    is_solved = resoltrisup(n, A, b, temp, tol);
    if (is_solved == 0) {
        for (i = 0; i < n; ++i) {
            b[i] = temp[i];
        }
    }
    free(temp);
    return is_solved;
}

int gausspiv(int n, double **A, double *b, double tol) {
    int i, j, k, pivot, is_solved;
    double *tempB = (double *) malloc(sizeof(double) * n);
    for (k = 0; k < n - 1; k++) {
        for (i = k + 1; i < n; i++) {
            if (fabs(A[k][k]) < fabs(A[i][k])) {
                for (j = 0; j < n; j++) {
                    double temp;
                    temp = A[k][j];
                    A[k][j] = A[i][j];
                    A[i][j] = temp;
                }
            }
        }
        for (i = k + 1; i < n; i++) {
            double term = A[i][k] / A[k][k];
            for (j = 0; j < n; j++) {
                A[i][j] = A[i][j] - term * A[k][j];
            }
        }
    }


    is_solved = resoltrisup(n, A, b, tempB, tol);
    if (is_solved == 0) {
        for (i = 0; i < n; ++i) {
            b[i] = tempB[i];
        }
    }
    free(tempB);
    return is_solved;
}


double **invertMat(double **A, int n) {
    double det = determinant(A, n);
    double **inv = matCofact(A, n);
    inv = transposar(inv, n, n);
    divMat(inv, det, n);
    return inv;
}


double **checkLU(int n, double **a, double **acp) {
    double **uInv = invertMat(acp, n);
    double **L = prodMatMat(a, uInv, n, n);
    double **A = prodMatMat(L, acp, n, n);
    return A;
}


