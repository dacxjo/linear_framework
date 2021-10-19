//
// Created by David Blandon on 13/10/21.
//

#ifndef BLANDONTORREZDAVID_PRAC1_FUNS_LINALG_H
#define BLANDONTORREZDAVID_PRAC1_FUNS_LINALG_H

/**
 * Resol una matriu triangular superior donada
 * utilitzant substitució enrere (Backwards Substitution)
 * @param n - Dimensió de la matriu
 * @param A - Matriu
 * @param b - Vector terme independent
 * @param x - Vector solució
 * @param tol - Tolerància
 * @return - 0 Si s'ha pogut resoldre el sistema , 1 en cas contrari
 */
int resoltrisup(int n, double **A, double *b, double *x, double tol);

/**
 * Comprova si una matriu donada és triangular superior
 * @param A - Matriu a comprovar
 * @param n - Dimension de la matriu
 * @return - 0 si la matriu no és triangular superior i 1 al contrari
 */
_Bool checktrisup(double **A, int n);

/**
 * Resol una matriu triangular inferior donada
 * utilitzant substitució cap endavant (Forward Substitution)
 * @param n - Dimensió de la matriu
 * @param A - Matriu
 * @param b - Vector terme independent
 * @param x - Vector solució
 * @param tol - Tolerància
 * @return - 0 Si s'ha pogut resoldre el sistema , 1 en cas contrari
 */
int resoltriinf(int n, double **A, double *b, double *x, double tol);

/**
 * Comprova si una matriu donada és triangular inferior
 * @param A - Matriu a comprovar
 * @param n - Dimension de la matriu
 * @return - 0 si la matriu no és triangular inferior i 1 al contrari
 */
_Bool checktriinf(double **A, int n);

/**
* Calc el producte escalar entre dos vectors donats
* @param n - Dimensió del vector
* @param x - Vector X
* @param y - Vector Y
* @return - Producte escalar
*/
double prod_esc(int n, double *x, double *y);

/**
 * Producte entre una matriu i un vector
 * @param M - Matriu
 * @param x - Vector
 * @param n - Dimensió
 * @return
 */
double *prodMatVect(double **M, double *x,int n);

/**
* Calcula la solució d'un matriu per mitjà del. mètode de gauss sense pivots
* @param n - Dimensió del vector
* @param A - Matriu
* @param b - Vector terme independent
* @param tol - Tolerància
* @return - 0 Si s'ha pogut resoldre el sistema , 1 en cas contrari
*/
int gauss(int n, double **A, double *b, double tol);

#endif //BLANDONTORREZDAVID_PRAC1_FUNS_LINALG_H
