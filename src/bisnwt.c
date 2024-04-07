#include <stdio.h>
#include <math.h>
#include <stdbool.h>

int bisnwt (double a, double b, double *arr,double *dlt, double tol, int maxit, double (*f)(double,void*), double (*df)(double,void*), void *prm){

    double C, a_im, b_im, c_im, xn, xn1, fx, dfx;
    

    while (true) {

        while( fabs(b-a) > *dlt ){ //Apliquem el metode de la biseccio sempre que la llargada de l'interval [a,b] sigui mes gran que δ

            C = (a+b)/2;
            a_im = (*f)(a, prm);
            c_im = (*f)(C,prm);

            if (a_im * c_im <= 0) {
                b = C;
            } else {
                a = C;
            }
        }

        if (*dlt <= tol){
            *arr = C;
            return -1;
        } else {
            xn = C;
            for (int i = 1; i <= maxit; i++) {
                fx = (*f)(xn, prm);
                dfx = (*df)(xn, prm);
                xn1 = xn - (fx / dfx);
                if (fabs(xn - xn1) < tol) {
                    *arr = xn1;
                    return i;
                }
            xn=xn1;
        }
        *dlt /= 2;
        }
    }
}