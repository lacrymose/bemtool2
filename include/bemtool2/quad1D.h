#ifndef QUAD_H
#define QUAD_H

#include "calculus.h"
#include "triangle_dunavant_rule.hpp"

/*=========================================
||  Quadrature dans le tetrahedre unite  ||
=========================================*/
/*
int keast_order_num ( int rule );
void keast_rule (int rule, int order_num, 
		 double xyz[], double w[] );

void quad3D(const int& N, vectR3& x, vectR& w){
  
  int order = keast_order_num(N);  
  resize(x,order);
  resize(w,order);
  
  double qw[order];
  double qp[3*order];
  keast_rule(N,order,qp,qw);

  for(int j=0; j<order; j++){
    x[j][0] = qp[3*j];
    x[j][1] = qp[3*j+1];
    x[j][2] = qp[3*j+2];;
    w[j] = qw[j];}
  
}
*/

/*=========================================
||  Quadrature dans le triangle unite    ||
=========================================*/

// int dunavant_order_num ( int rule );
// void dunavant_rule ( int rule, int order_num, 
// 		     double xy[], double w[] );

void quad2D(const int& N, std::vector<R2>& x, std::vector<Real>& w);


/*=======================================
||  Quadrature sur le segment unite    ||
=======================================*/

void quad1D(const int& order, std::vector<Real>& x, std::vector<Real>& w);


/*===========================================
||  Composite Gauss-Legendre quadrature    ||
||  sur le segment unite (strategie HP)    ||
||  adaptee a une integrande singuliere    || 
||  en 0                                   ||
===========================================*/

//cgauleg_redux(double *x, double *w, int q){
void hp_quad1D(const int& order, std::vector<Real>& x, std::vector<Real>& w);



#endif
