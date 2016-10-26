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

void quad2D(const int& N, std::vector<R2>& x, std::vector<Real>& w){
  
  int order = dunavant_order_num(N);
  x.resize(order);
  w.resize(order);
  
  double qw[order];
  double qp[2*order];
  dunavant_rule(N,order,qp,qw);
  
  for(int j=0; j<order; j++){
    x[j][0] = qp[2*j];
    x[j][1] = qp[2*j+1];
    w[j] = 0.5*qw[j];}
  
}


/*=======================================
||  Quadrature sur le segment unite    ||
=======================================*/

void quad1D(const int& order, std::vector<Real>& x, std::vector<Real>& w) {
  
  x.resize(order); 
  w.resize(order);
  
  double d1,d2pn,d3pn,d4pn,dp,dpn,e1,fx,h,p,pk,pkm1,pkp1,t,u,v,x0,xtemp;
  int i,iback,k,m,mp1mi,ncopy,nmove;
  double pi = 3.141592653589793;
  
  if(order<1){
    std::cerr << "quadrature doit etre d'ordre >1";
    exit(1);}
  
  e1 = ( double ) ( order * ( order + 1 ) );
  m = ( order + 1 ) / 2;
  
  for ( i = 1; i <= m; i++ )
    {
      mp1mi = m + 1 - i;
      t = ( double ) ( 4 * i - 1 ) * pi / ( double ) ( 4 * order + 2 );
      x0 = cos(t)*( 1.0 - ( 1.0 - 1.0 / ( double ) ( order ) ) 
		    / ( double ) ( 8 * order * order ) );
      pkm1 = 1.0;
      pk = x0;
      
      for ( k = 2; k <= order; k++ ){
	pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
	pkm1 = pk;
	pk = pkp1;}
      
      d1 = ( double ) ( order ) * ( pkm1 - x0 * pk );
      dpn = d1 / ( 1.0 - x0 * x0 );
      d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );
      d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );
      d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );
      u = pk / dpn;
      v = d2pn / dpn;
      
      //  Initial approximation H:
      h = -u*( 1.0 + 0.5 * u * ( v + u*( v*v - d3pn / (3.0*dpn) ) ) );
      
      //  Refine H using one step of Newton's method:
      p = pk + h*( dpn + 0.5 * h * ( d2pn + h / 3.0 
				     * ( d3pn + 0.25 * h * d4pn ) ) );
      
      dp = dpn + h*( d2pn + 0.5 * h*( d3pn + h * d4pn / 3.0 ) );
      h = h - p / dp;
      xtemp = x0 + h;
      x[mp1mi-1] = xtemp;
      
      fx = d1-h*e1*(pk + 0.5*h*( dpn + h / 3.0 
				 *(d2pn + 0.25*h*( d3pn + 0.2*h*d4pn ))));
      w[mp1mi-1] = 2.0 * ( 1.0 - xtemp * xtemp ) / ( fx * fx );
    }
  
  if( (order%2)==1 ){x[0] = 0.0;}
  
  //  Shift the data up.
  nmove = ( order + 1 ) / 2;
  ncopy = order - nmove;
  
  for ( i = 1; i <= nmove; i++ ){
    iback = order + 1 - i;
    x[iback-1] = x[iback-ncopy-1];
    w[iback-1] = w[iback-ncopy-1];}

  //  Reflect values for the negative abscissas.
  for ( i = 1; i <= order - nmove; i++ ){
    x[i-1] = - x[order-i];
    w[i-1] = w[order-i];}
  
  for(i=0; i<order; i++){
    x[i]   = (1+x[i])*0.5;
    w[i] = w[i]*0.5;}
  

}


/*===========================================
||  Composite Gauss-Legendre quadrature    ||
||  sur le segment unite (strategie HP)    ||
||  adaptee a une integrande singuliere    || 
||  en 0                                   ||
===========================================*/

//cgauleg_redux(double *x, double *w, int q){
void hp_quad1D(const int& order, std::vector<Real>& x, std::vector<Real>& w)
{
    
  int i, j, ii, qq, n;
  Real sigma, xl, xr;
  std::vector<Real> _w, _x;
  n = order*(order+1)/2-1;
  x.resize(n); w.resize(n);
  
  /* Compute composite Gauss-Legendre quadrature rule */
  sigma = (sqrt(2.0)-1.0)*(sqrt(2.0)-1.0);
  
  xl = sigma; 
  xr = 1.0; 
  ii = n-1; 
  qq = order;
  for(i=1; i<order; i++){    
    
    /* Initialize Gauss-Legendre quadrature */
    quad1D(qq,_x,_w);  
    
    /* Transform quadrature points */
    for(j=qq-1; j>=0; j--){
      w[ii] = (xr-xl)*_w[j];
      x[ii] = (xr-xl)*_x[j]+xl;
      ii--;}    
    
    qq--;
    xr = xl;
    xl = xl*sigma;    
  }
  
}


#endif
