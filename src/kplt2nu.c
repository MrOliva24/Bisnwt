#include <stdio.h>
#include <math.h>
#include "bisnwt.h"
#define PI 3.14159265358979323846


double kepler_eq (double , void *);
double kepler_der (double , void *); 


int main (int argc, char *argv[]) {

    //double PI = acos(-1.0);


   double e, T, M0, tf;
   int nt;
   if (argc<6
         || sscanf(argv[1], "%lf", &e)!=1
         || sscanf(argv[2], "%lf", &T)!=1
         || sscanf(argv[3], "%lf", &M0)!=1
         || sscanf(argv[4], "%lf", &tf)!=1
         || sscanf(argv[5], "%d", &nt)!=1
      ) {
      fprintf(stderr,"%s e T M0 tf nt\n", argv[0]);
      return -1;
   }

   double M, ti, a, b, E, v, iter, tp;
   int n;
   M = M0;
   double prm[2] = {e, M};

   double dlt = 2.5, tol=1e-12, arr;
   int maxit = 10;


   tp = -(M*T)/(2*PI);
   //prm[0] = e;
   //prm[1] = M;
   double interval = tf/nt;

   for (int i=0; i<nt; i++){
    
        ti = i* interval; 
        prm[1] = 2*PI*((ti - tp)/ T); //M
        a = prm[1] - PI;
        b = prm[1] + PI; // [Μ-π,Μ+π]

        iter = bisnwt(a,b,&arr,&dlt,tol,maxit,&kepler_eq,&kepler_der,prm);

        E = arr;
        n = floor(E / (2 * PI));
        v = acos((prm[0]-cos(E))/(prm[0]*cos(E)-1));


        if (E >= PI + (2 * PI * n) && E <= 2*PI + (2 * PI * n)) v = 2*PI - v;

        
        v += 2 * PI * n;
        printf("%.16g %.16g %.16g\n", ti, prm[1], v);

    }

   return 0;
}

double kepler_eq (double x, void * prm){

   double e = ((double *)prm)[0];
   double M = ((double *)prm)[1];
   return x - e * sin(x)-M;
}

double kepler_der (double x, void * prm){
   
   double e = ((double *)prm)[0];
   return 1 - e * cos(x);
}