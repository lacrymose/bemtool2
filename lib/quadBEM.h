#ifndef QUADBEM_H
#define QUADBEM_H

#include "quad1D.h"
#include <iostream>
#include <fstream>
using namespace std;


template <int dim>
class quadBEM;


template <>
class quadBEM<2>{

 public: 
  typedef R2   qp_t;
  
 private: 
  vector<R2>   x_[4];
  vector<R2>   y_[4];
  vector<Real> w_[4];
  
 public:  
  const vector<R2>&   x(const int& rule) const {return x_[rule];}
  const vector<R2>&   y(const int& rule) const {return y_[rule];}
  const vector<Real>& w(const int& rule) const {return w_[rule];}
  
  quadBEM<2>(const int& order){
    
    R2 x,y; Real w,dw;
    vector<Real> t,dt;
    quad1D(order,t,dt);
    int nt = t.size();
    
    vector<R2>   qp;
    vector<Real> qw;
    quad2D(order,qp,qw);
    int nq = qp.size();
    
    //================================//
    //   Cas 0: triangles disjoints   //
    //================================//
    
    
    // Dans le cas ou les triangles sont 
    // disjoints, on utilise deux regles de 
    // quadrature de Dunavant tensorisees 
    for(int j=0; j<nq; j++){
      for(int k=0; k<nq; k++){
	
 	x[0] = 1.-qp[j][0];
	x[1] = qp[j][1];
	y[0] = 1.-qp[k][0];
	y[1] = qp[k][1];
	w    = qw[j]*qw[k];
	
	x_[0].push_back(x);
	y_[0].push_back(y);
	w_[0].push_back(w);
	
      }
    }
    

    /*
    for(int j0=0; j0<nt; j0++){
      for(int j1=0; j1<nt; j1++){
	for(int j2=0; j2<nt; j2++){
	  for(int j3=0; j3<nt; j3++){
	    
	    x[0] = t[j0];
	    x[1] = t[j0]*t[j1];
	    y[0] = t[j2];
	    y[1] = t[j2]*t[j3];
	    w    = t[j0]*t[j2]*dt[j0]*dt[j1]*dt[j2]*dt[j3];
	    
	    x_[0].push_back(x);
	    y_[0].push_back(y);
	    w_[0].push_back(w);
	    
	  }
	}	
      }
    }
    */    
    
    //==================================//
    
    // Dans les cas de triangles adjacents
    // on utlises les changements de variable 
    // de Sauter-Schwab avec 4 regles de 
    // quadrature de Gauss-Legendre tensorisees.
    
    for(int j0=0; j0<nt; j0++){
      for(int j1=0; j1<nt; j1++){
	for(int j2=0; j2<nt; j2++){
	  for(int j3=0; j3<nt; j3++){
	    
	    const Real& s  = t[j0];
	    const Real& e1 = t[j1];
	    const Real& e2 = t[j2];
	    const Real& e3 = t[j3];
	    
	    dw = dt[j0]*dt[j1]*dt[j2]*dt[j3];	    
	    
	    //================================//
	    //   Cas 1: un sommet commun      //
	    //================================//
	    
	    x[0] = s;
	    x[1] = s*e1;
	    y[0] = s*e2;
	    y[1] = s*e2*e3;
	    w    = s*s*s*e2*dw;
	    
	    x_[1].push_back(x);
	    y_[1].push_back(y);
	    w_[1].push_back(w);

	    //------------//
	    
	    x[0] = s*e2;
	    x[1] = s*e2*e3;
	    y[0] = s;
	    y[1] = s*e1;
	    w    = s*s*s*e2*dw;

	    x_[1].push_back(x);
	    y_[1].push_back(y);
	    w_[1].push_back(w);

	    //================================//
	    //   Cas 2: deux sommets communs  //
	    //================================//	    
	    
	    x[0] = s;
	    x[1] = s*e1*e3;
	    y[0] = s*(1-e1*e2);
	    y[1] = s*e1*(1-e2);
	    w    = s*s*s*e1*e1*dw;
	    
	    x_[2].push_back(x);
	    y_[2].push_back(y);
	    w_[2].push_back(w);

	    //---------------//

	    x[0] = s;
	    x[1] = s*e1;
	    y[0] = s*(1.-e1*e2*e3);
	    y[1] = s*e1*e2*(1-e3);
	    w    = s*s*s*e1*e1*e2*dw;
	    
	    x_[2].push_back(x);
	    y_[2].push_back(y);
	    w_[2].push_back(w);

	    //---------------//

	    x[0] = s*(1.-e1*e2);
	    x[1] = s*e1*(1.-e2);
	    y[0] = s;
	    y[1] = s*e1*e2*e3;
	    w    = s*s*s*e1*e1*e2*dw;
	    
	    x_[2].push_back(x);
	    y_[2].push_back(y);
	    w_[2].push_back(w);

	    //---------------//

	    x[0] = s*(1-e1*e2*e3);
	    x[1] = s*e1*e2*(1.-e3);
	    y[0] = s;
	    y[1] = s*e1;
	    w    = s*s*s*e1*e1*e2*dw;
	    
	    x_[2].push_back(x);
	    y_[2].push_back(y);
	    w_[2].push_back(w);

	    //---------------//

	    x[0] = s*(1.-e1*e2*e3);
	    x[1] = s*e1*(1.-e2*e3);
	    y[0] = s;
	    y[1] = s*e1*e2;
	    w    = s*s*s*e1*e1*e2*dw;
	    
	    x_[2].push_back(x);
	    y_[2].push_back(y);
	    w_[2].push_back(w);
	    
	    //================================//
	    //   Cas 3: trois sommets communs //
	    //================================//	    

	    x[0] = s;
	    x[1] = s*(1.-e1+e1*e2 );
	    y[0] = s*(1.-e1*e2*e3);
	    y[1] = s*(1.-e1);
	    w    = s*s*s*e1*e1*e2*dw;
	    
	    x_[3].push_back(x);
	    y_[3].push_back(y);
	    w_[3].push_back(w);

	    //------------//

	    x[0] = s*(1.-e1*e2*e3);
	    x[1] = s*(1.-e1);
	    y[0] = s;
	    y[1] = s*(1.-e1+e1*e2);
	    w    = s*s*s*e1*e1*e2*dw;
	    
	    x_[3].push_back(x);
	    y_[3].push_back(y);
	    w_[3].push_back(w); 
	    
	    //------------//

	    x[0] = s;
	    x[1] = s*e1*(1.-e2+e2*e3);
	    y[0] = s*(1.-e1*e2);
	    y[1] = s*e1*(1.-e2);
	    w    = s*s*s*e1*e1*e2*dw;
	    
	    x_[3].push_back(x);
	    y_[3].push_back(y);
	    w_[3].push_back(w);

	    //------------//

	    x[0] = s*(1.-e1*e2);
	    x[1] = s*e1*(1.-e2);
	    y[0] = s;
	    y[1] = s*e1*(1.-e2+e2*e3);
	    w    = s*s*s*e1*e1*e2*dw;
	    
	    x_[3].push_back(x);
	    y_[3].push_back(y);
	    w_[3].push_back(w);

	    //------------//

	    x[0] = s*(1.-e1*e2*e3);
	    x[1] = s*e1*(1-e2*e3);
	    y[0] = s;
	    y[1] = s*e1*(1.-e2);
	    w    = s*s*s*e1*e1*e2*dw;
	    
	    x_[3].push_back(x);
	    y_[3].push_back(y);
	    w_[3].push_back(w);

	    //------------//

	    x[0] = s;
	    x[1] = s*e1*(1.-e2);
	    y[0] = s*(1.-e1*e2*e3);
	    y[1] = s*e1*(1.-e2*e3);
	    w    = s*s*s*e1*e1*e2*dw;
	    
	    x_[3].push_back(x);
	    y_[3].push_back(y);
	    w_[3].push_back(w);

	  }
	}	
      }
    }


    //==================================//
    // Transformation pour se ramener a
    // l'element de reference usuel
    //==================================//
    R2x2 B; B(0,0)=1.;B(1,1)=1.;B(0,1)=-1.;
    for(int q=0; q<4; q++){
      for(int j=0; j<w_[q].size(); j++){
	x_[q][j] = B*x_[q][j];
	y_[q][j] = B*y_[q][j];
      }    
    }

    
    
  }  
  
};






template <>
class quadBEM<1>{

 public:
  typedef Real qp_t;
  
 private: 
  vector<Real> x_[3];
  vector<Real> y_[3];
  vector<Real> w_[3];
  
 public:  
  const vector<Real>& x(const int& rule) const {return x_[rule];}
  const vector<Real>& y(const int& rule) const {return y_[rule];}
  const vector<Real>& w(const int& rule) const {return w_[rule];}
  
  quadBEM<1>(const int& order){
    
    Real wxi, xxi, weta, xeta;
    vector<Real> w1, x1, w2, x2;
    int n, q, q1, q2;
    
    n = (int) ceil((order+1)/2.0); // mefiance...
    quad1D(n,x1,w1); q1 = x1.size();
    hp_quad1D(n,x2,w2); q2 = x2.size();

    
    //////////////////////////
    //   Elements disjoints  
    for(int i=0; i<q1; i++){
      wxi = w1[i];
      xxi = x1[i];
      for(int j=0; j<q1; j++){
	weta = w1[j];
	xeta = x1[j];
	
	w_[0].push_back(weta*wxi);
	x_[0].push_back(xxi);
	y_[0].push_back(xeta);
      }
    }

    
    ///////////////////////////
    //  Un seul noeud commun
    for(int i=0; i<q2; i++){
      wxi = w2[i];
      xxi = x2[i];
      for(int j=0; j<q1; j++){
	weta = w1[j];
	xeta = x1[j];

	w_[1].push_back(xxi*weta*wxi);
	x_[1].push_back(xxi);
	y_[1].push_back(xeta*xxi);
	
	w_[1].push_back(xxi*weta*wxi);
	x_[1].push_back(xxi*xeta);
	y_[1].push_back(xxi);
		
      }
    }

    
    //////////////////////////
    //  Elements identiques
    for(int i=0; i<q1; i++){
      wxi = w1[i];
      xxi = x1[i];
      for(int j=0; j<q2; j++){
	weta = w2[j];
	xeta = x2[j];

	w_[2].push_back( xxi*weta*wxi  );
	x_[2].push_back( xxi*(1.-xeta) );
	y_[2].push_back( xxi );
	
	w_[2].push_back( xxi*weta*wxi );
	x_[2].push_back( xxi );
	y_[2].push_back( xxi*(1.-xeta) );
	
      }
    }
    
    
  }


};










#endif
