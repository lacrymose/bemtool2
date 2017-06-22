#ifndef QUADPOT_H
#define QUADPOT_H

#include "quad.h"
#include <iostream>
#include <fstream>


template <int dim>
class quadPOT;

template <>
class quadPOT<2>{

 public:
  typedef R2   qp_t;

 private:
  std::vector<R2>   x_;
  std::vector<Real> w_;

 public:
  const std::vector<R2>&   x() const {return x_;}
  const std::vector<Real>& w() const {return w_;}

  quadPOT<2>(const int& order){

    R2 x,y; Real w,dw;
    std::vector<Real> t,dt;
    quad1D(order,t,dt);
    int nt = t.size();

    std::vector<R2>   qp;
    std::vector<Real> qw;
    quad2D(order,qp,qw);
    int nq = qp.size();

    //================================//
    //   Cas 0: triangles disjoints   //
    //================================//

    for(int j=0; j<nq; j++){
      x[0] = 1.-qp[j][0];
	    x[1] = qp[j][1];
	    w   = qw[j];
      x_.push_back(x);
       w_.push_back(w);
    }


    //==================================//

    // Dans les cas de triangles adjacents
    // on utlises les changements de variable
    // de Sauter-Schwab avec 4 regles de
    // quadrature de Gauss-Legendre tensorisees.

  //   for(int j0=0; j0<nt; j0++){
  //     for(int j1=0; j1<nt; j1++){
	// for(int j2=0; j2<nt; j2++){
	//   for(int j3=0; j3<nt; j3++){
  //
	//     const Real& s  = t[j0];
	//     const Real& e1 = t[j1];
	//     const Real& e2 = t[j2];
	//     const Real& e3 = t[j3];
  //
	//     dw = dt[j0]*dt[j1]*dt[j2]*dt[j3];
  //
	//     //================================//
	//     //   Cas 1: un sommet commun      //
	//     //================================//
  //
	//     x[0] = s;
	//     x[1] = s*e1;
	//     y[0] = s*e2;
	//     y[1] = s*e2*e3;
	//     w    = s*s*s*e2*dw;
  //
	//     x_[1].push_back(x);
	//     y_[1].push_back(y);
	//     w_[1].push_back(w);
  //
	//     //------------//
  //
	//     x[0] = s*e2;
	//     x[1] = s*e2*e3;
	//     y[0] = s;
	//     y[1] = s*e1;
	//     w    = s*s*s*e2*dw;
  //
	//     x_[1].push_back(x);
	//     y_[1].push_back(y);
	//     w_[1].push_back(w);
  //
	//     //================================//
	//     //   Cas 2: deux sommets communs  //
	//     //================================//
  //
	//     x[0] = s;
	//     x[1] = s*e1*e3;
	//     y[0] = s*(1-e1*e2);
	//     y[1] = s*e1*(1-e2);
	//     w    = s*s*s*e1*e1*dw;
  //
	//     x_[2].push_back(x);
	//     y_[2].push_back(y);
	//     w_[2].push_back(w);
  //
	//     //---------------//
  //
	//     x[0] = s;
	//     x[1] = s*e1;
	//     y[0] = s*(1.-e1*e2*e3);
	//     y[1] = s*e1*e2*(1-e3);
	//     w    = s*s*s*e1*e1*e2*dw;
  //
	//     x_[2].push_back(x);
	//     y_[2].push_back(y);
	//     w_[2].push_back(w);
  //
	//     //---------------//
  //
	//     x[0] = s*(1.-e1*e2);
	//     x[1] = s*e1*(1.-e2);
	//     y[0] = s;
	//     y[1] = s*e1*e2*e3;
	//     w    = s*s*s*e1*e1*e2*dw;
  //
	//     x_[2].push_back(x);
	//     y_[2].push_back(y);
	//     w_[2].push_back(w);
  //
	//     //---------------//
  //
	//     x[0] = s*(1-e1*e2*e3);
	//     x[1] = s*e1*e2*(1.-e3);
	//     y[0] = s;
	//     y[1] = s*e1;
	//     w    = s*s*s*e1*e1*e2*dw;
  //
	//     x_[2].push_back(x);
	//     y_[2].push_back(y);
	//     w_[2].push_back(w);
  //
	//     //---------------//
  //
	//     x[0] = s*(1.-e1*e2*e3);
	//     x[1] = s*e1*(1.-e2*e3);
	//     y[0] = s;
	//     y[1] = s*e1*e2;
	//     w    = s*s*s*e1*e1*e2*dw;
  //
	//     x_[2].push_back(x);
	//     y_[2].push_back(y);
	//     w_[2].push_back(w);
  //
	//     //================================//
	//     //   Cas 3: trois sommets communs //
	//     //================================//
  //
	//     x[0] = s;
	//     x[1] = s*(1.-e1+e1*e2 );
	//     y[0] = s*(1.-e1*e2*e3);
	//     y[1] = s*(1.-e1);
	//     w    = s*s*s*e1*e1*e2*dw;
  //
	//     x_[3].push_back(x);
	//     y_[3].push_back(y);
	//     w_[3].push_back(w);
  //
	//     //------------//
  //
	//     x[0] = s*(1.-e1*e2*e3);
	//     x[1] = s*(1.-e1);
	//     y[0] = s;
	//     y[1] = s*(1.-e1+e1*e2);
	//     w    = s*s*s*e1*e1*e2*dw;
  //
	//     x_[3].push_back(x);
	//     y_[3].push_back(y);
	//     w_[3].push_back(w);
  //
	//     //------------//
  //
	//     x[0] = s;
	//     x[1] = s*e1*(1.-e2+e2*e3);
	//     y[0] = s*(1.-e1*e2);
	//     y[1] = s*e1*(1.-e2);
	//     w    = s*s*s*e1*e1*e2*dw;
  //
	//     x_[3].push_back(x);
	//     y_[3].push_back(y);
	//     w_[3].push_back(w);
  //
	//     //------------//
  //
	//     x[0] = s*(1.-e1*e2);
	//     x[1] = s*e1*(1.-e2);
	//     y[0] = s;
	//     y[1] = s*e1*(1.-e2+e2*e3);
	//     w    = s*s*s*e1*e1*e2*dw;
  //
	//     x_[3].push_back(x);
	//     y_[3].push_back(y);
	//     w_[3].push_back(w);
  //
	//     //------------//
  //
	//     x[0] = s*(1.-e1*e2*e3);
	//     x[1] = s*e1*(1-e2*e3);
	//     y[0] = s;
	//     y[1] = s*e1*(1.-e2);
	//     w    = s*s*s*e1*e1*e2*dw;
  //
	//     x_[3].push_back(x);
	//     y_[3].push_back(y);
	//     w_[3].push_back(w);
  //
	//     //------------//
  //
	//     x[0] = s;
	//     x[1] = s*e1*(1.-e2);
	//     y[0] = s*(1.-e1*e2*e3);
	//     y[1] = s*e1*(1.-e2*e3);
	//     w    = s*s*s*e1*e1*e2*dw;
  //
	//     x_[3].push_back(x);
	//     y_[3].push_back(y);
	//     w_[3].push_back(w);
  //
	//   }
	// }
  //     }
  //   }
  //
  //
    //==================================//
    // Transformation pour se ramener a
    // l'element de reference usuel
    //==================================//
    R2x2 B; B(0,0)=1.;B(1,1)=1.;B(0,1)=-1.;

    for(int j=0; j<w_.size(); j++){
		x_[j] = B*x_[j];
      }




  }

};




template <>
class quadPOT<1>{

 public:
  typedef Real qp_t;

 private:
  std::vector<Real> x_;
  std::vector<Real> w_;

 public:
  const std::vector<Real>& x() const {return x_;}
  const std::vector<Real>& w() const {return w_;}

  quadPOT<1>(const int& order){

    Real wxi, xxi, weta, xeta;
    std::vector<Real> w1, x1, w2, x2;
    int n, q, q1, q2;

    n = (int) ceil((order+1)/2.0); // mefiance...
    quad1D(n,x1,w1); q1 = x1.size();
    hp_quad1D(n,x2,w2); q2 = x2.size();


    //////////////////////////
    //   Elements disjoints
    for(int i=0; i<q1; i++){
      wxi = w1[i];
      xxi = x1[i];

      w_.push_back(wxi);
      x_.push_back(xxi);
    }


  //   ///////////////////////////
  //   //  Un seul noeud commun
  //   for(int i=0; i<q2; i++){
  //     wxi = w2[i];
  //     xxi = x2[i];
  //     for(int j=0; j<q1; j++){
	// weta = w1[j];
	// xeta = x1[j];
  //
	// w_[1].push_back(xxi*weta*wxi);
	// x_[1].push_back(xxi);
	// y_[1].push_back(xeta*xxi);
  //
	// w_[1].push_back(xxi*weta*wxi);
	// x_[1].push_back(xxi*xeta);
	// y_[1].push_back(xxi);
  //
  //     }
  //   }
  //
  //
  //   //////////////////////////
  //   //  Elements identiques
  //   for(int i=0; i<q1; i++){
  //     wxi = w1[i];
  //     xxi = x1[i];
  //     for(int j=0; j<q2; j++){
	// weta = w2[j];
	// xeta = x2[j];
  //
	// w_[2].push_back( xxi*weta*wxi  );
	// x_[2].push_back( xxi*(1.-xeta) );
	// y_[2].push_back( xxi );
  //
	// w_[2].push_back( xxi*weta*wxi );
	// x_[2].push_back( xxi );
	// y_[2].push_back( xxi*(1.-xeta) );
  //
  //     }
  //   }


  }


};










#endif
