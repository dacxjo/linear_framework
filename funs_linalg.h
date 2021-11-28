#ifndef BLANDONTORREZDAVID_PRAC1_FUNS_LINALG_H
#define BLANDONTORREZDAVID_PRAC1_FUNS_LINALG_H

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
void prodMatVect(double **M, double *x,double *destination,int n);

/**
 * Producte entre dues matrius
 * @param A - Matriu A
 * @param B - Matriu B
 * @param n - Files
 * @param n - Columnes
 */
void prodMatMat(double **A, double **B,double **destination,int n,int m);

/**
 * Suma entre dues matrius
 * @param A - Matriu A
 * @param B - Matriu B
 * @param n - Files
 * @param n - Columnes

 */
void sumMatMat(double **A, double **B,double **destination,int n);

/**
 * Resta entre dues matrius
 * @param A - Matriu A
 * @param B - Matriu B
 * @param n - Files
 * @param n - Columnes

 */
void restMatMat(double **A, double **B,double **destination,int n);

/**
 * Resta entre dues vectors
 * @param V - Vector A
 * @param W - Vector B
 * @param destination - Vector on es guarda el resultat
 * @param n - Dimensio

 */
void restVectVect(double *V, double *W,double *destination,int n);

/**
 * Donada una matriu, retorna la seva transposada
 * @param M - Matriu
 * @param m - Nombre de files
 * @param n - Nombre de columnes
 * @return
 */
double **transposar(double **M, int m,int n);

/**
 * Genera una matriu Identitat
 * @param n - Ordre de la matriu
 * @param M - Matriu
 */
void genMatId(int n,double **M);

/**
 * Genera una matriu nula
 * @param n - Ordre de la matriu
 * @param M - Matriu
 */
void genMatNul(int n,double **M);

/**
 * Genera una matriu de Hessenberg
 * @param n - Ordre de la matriu
 * @param M - Matriu
 */
void genMatFrankHessenberg(int n,double **M);

/**
 * Genera un vector nul
 * @param n - Dimensio
 * @param V - Vector
 */
void genVectNul(int n,double *V);

/**
 * Genera un vector identitat
 * @param n - Dimensio
 * @param V - Vector
 */
void genVectId(int n,double *V);

/**
 * Calcula la norma euclidiana d'un vector
 * @param n - Dimensio
 * @param V - Vector
 */
double calcNormEucl(int n,double *V);

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
 * @param n - Dimensio de la matriu
 * @return - 0 si la matriu no és triangular superior i 1 al contrari
 */
int checktrisup(double **A, int n);

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
 * @param n - Dimensio de la matriu
 * @return - 0 si la matriu no és triangular inferior i 1 al contrari
 */
int checktriinf(double **A, int n);

/**
* Calcula la solució d'un matriu per mitjà del. mètode de gauss sense pivots
* @param n - Dimensió de la matriu
* @param A - Matriu
* @param b - Vector terme independent
* @param tol - Tolerància
* @return - 0 Si s'ha pogut resoldre el sistema , 1 en cas contrari
*/
int gauss(int n, double **A, double *b, double tol);

/**
* Calcula la solució d'un matriu per mitjà del mètode de gauss amb pivots
* @param n - Dimensió de la matriu
* @param A - Matriu
* @param v - Vector terme independent
* @param tol - Tolerància
* @return - 0 Si s'ha pogut resoldre el sistema , 1 en cas contrari
*/
int gausspiv(int n, double **A, double *v, double tol);

/**
* Versio modificada de Gauss per obtenir una matriu alterada
* @param n - Dimensió de la matriu
* @param A - Matriu
*/
void gaussLU(int n, double **A);


/**
* Factoritza una matriu alterada per a obtenir les matrius L i U
* @param n - Ordre
* @param A - Matriu alterada
* @param L - Matriu Lower
* @param U - Upper Matriu
*/
void luDecompose(int n, double **acp, double **L, double **U);

/**
* @param n - Dimensió de la matriu
* @param A - Matriu original
* @param acp - Matriu alterada
*/
double checkLU(int n, double **a, double **acp);
#endif
