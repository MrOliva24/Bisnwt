#include <stdio.h>
#include <math.h>
#include <stdbool.h>


int bisnwt (double a, double b, double *arr,double *dlt, double tol, int maxit, double (*f)(double,void*), double (*df)(double,void*), void *prm){

    double C, a_im, b_im, c_im, xn, xn1, fx, dfx;
    
    printf("\ndelta %.16g\n", *dlt);
    while (true) {

        while( fabs(b-a) > *dlt ){ //Apliquem el metode de la biseccio sempre que la llargada de l'interval [a,b] sigui mes gran que δ

            C = (a+b)/2;
            a_im = (*f)(a, prm);
            c_im = (*f)(C,prm);

            if (a_im * c_im < 0) {
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


double f (double x, void *prm) {
   return exp(x)-2;

}

double df (double x, void *prm) {
   return exp(x);
}

int main () {

    double dlt, a=-9, b=1, arr, tol=1e-12;
    int maxit=10, iter;
    double dlts[3] = {10, 2.5, .01};

    for (int i = 0; i<3; i++){
        dlt = dlts[i];
        iter = bisnwt(a,b,&arr,&dlt,tol,maxit,&f,&df,NULL);
        printf("Utilitzant dlt inicial = %.16g, f(x) = 0 -> x = %.16g\n Resultat assolit en %d iteracions\n\n",dlts[i], arr, iter);
    }

    return 0;    
}
