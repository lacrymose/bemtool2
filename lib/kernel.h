#ifndef KERNEL_H
#define KERNEL_H

#include "normal.h"
#include "dof.h"
#include <bessel.hpp>

using namespace boost::math;


Cplx H_0(Real x){return cyl_bessel_j(0,x)+ iu*cyl_neumann(0,x);}
Cplx H_1(Real x){return cyl_bessel_j(1,x)+ iu*cyl_neumann(1,x);}

/*=================================
||  NOYAUX INTEGRAUX PROPAGATIFS ||
=================================*/

//____________________________________
// Noyau test constant 
class CST_3D{  

  typedef R2 qp_t;
  Cplx val;

 public:  
  static const int dim = 2;
  CST_3D(const Real& k0) {};
  inline Cplx& ker(const R3& nx, const R3& ny, const R3& u){return val = Cplx(1.,0.);}

  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx, const R3& ny, const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
  {return val = Cplx(1.,0.);} 

  
};

//___________________________
// Single Layer Potential 2-D

class SLP_2D{  

  typedef Real qp_t;
  const   Real k;
  Real r;
  Cplx val;  

 public:  
  static const int dim = 1;
  
 SLP_2D(const Real& k0): k(k0) {};
  
  inline Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); return val = 0.25*iu*H_0(k*r); }
  
  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx, const R3& ny,   const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
  {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;} 
  
};

//___________________________
// Single Layer Potential 3-D

class SLP_3D{  

  typedef R2 qp_t;
  const Real k;
  Real r;
  Cplx val;  

 public:  
  static const int dim = 2;
  
 SLP_3D(const Real& k0): k(k0) {};
  
  inline Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); return val = exp(iu*k*r)/(4*pi*r); }

  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx, const R3& ny,   const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
  {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;} 
  
};

//___________________________
// Double Layer Potential 2-D

class DLP_2D{  
  
  typedef Real qp_t;
  const Real k;
  Real r;
  Cplx val;  
  
 public:  
  static const int dim = 1;
  
 DLP_2D(const Real& k0): k(k0) {};
  
  inline Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); return val = -(ny,x_y)*(1./r)*0.25*iu*k*H_1(k*r); }
  
  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx,  const R3& ny,  const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
  {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;} 
    
};


//___________________________
// Double Layer Potential 3-D

class DLP_3D{  

  typedef R2 qp_t;
  const Real k;
  Real r,r3;
  Cplx val;  

 public:  
  static const int dim = 2;

 DLP_3D(const Real& k0): k(k0) {};
  
  inline Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); r3 = r*r*r; 
    return val = (ny,x_y)*(iu*k*r-1.)*exp(iu*k*r)/(4*pi*r3); }
  
  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx,  const R3& ny,  const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
  {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;} 
    
};

//____________________________________
// Adjoint Double Layer Potential 2-D

class TDLP_2D{  
  
  typedef Real qp_t;
  const Real k;
  Real r;
  Cplx val;  
  
 public:  
  static const int dim = 1;
  
 TDLP_2D(const Real& k0): k(k0) {};
  
  inline Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); return val = -(nx,x_y)*(1./r)*0.25*iu*k*H_1(k*r); }
  
  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx,  const R3& ny,  const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
  {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;} 
    
};

//____________________________________
// Adjoint Double Layer Potential 3-D
class TDLP_3D{  

  typedef R2 qp_t;
  const Real k;
  Real r,r3;
  Cplx val;  

 public:  
  static const int dim = 2;

 TDLP_3D(const Real& k0): k(k0) {};

  Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); r3 = r*r*r;
    return val = (nx,x_y)*(iu*k*r-1.)*exp(iu*k*r)/(4*pi*r3); }

  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx,  const R3& ny,  const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
  {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;} 
  
};

//____________________________________
// Hypersingular Layer Potential 2-D

class HSP_2D{  
  
  typedef Real qp_t;
  const Real k,k2;
  Real r;
  Cplx val;  
  
 public:  
  static const int dim = 1;
  
 HSP_2D(const Real& k0): k(k0), k2(k0*k0) {};
  
  Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); return val = 0.25*iu*H_0(k*r);}
  
  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx,  const R3& ny,  const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py){
    val = ( vprod( nx, grad(phix, jx, px[kx]) ), vprod( ny, grad(phiy, jy, py[ky]) ) )*z;
    val+= -k2*(nx,ny)*phix(s,jx,kx)*phiy(t,jy,ky)*z;            
    return val;}   
  
};


//____________________________________
// Hypersingular Layer Potential 3-D
class HSP_3D{  
  
  typedef R2 qp_t;
  const Real k,k2;
  Real r;
  Cplx val;  
  
 public:  
  static const int dim = 2;
  
 HSP_3D(const Real& k0): k(k0), k2(k0*k0) {};

  Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); return val = exp(iu*k*r)/(4*pi*r);}
  
  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx,  const R3& ny,  const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py){
    val = ( vprod( nx, grad(phix, jx, px[kx]) ), vprod( ny, grad(phiy, jy, py[ky]) ) )*z;
    val+= -k2*(nx,ny)*phix(s,jx,kx)*phiy(t,jy,ky)*z;            
    return val;}   
  
};




/*===============================
||  CLASSE OPERATEUR INTEGRAL  ||
===============================*/


template <class space_x, class space_y, class kernel_t> class bem{

 public:
  
  static const int dim        = kernel_t::dim;
  static const int dim_loc_x  = space_x::dim_loc;
  static const int dim_loc_y  = space_y::dim_loc;
  
  typedef normal<dim>                     normal_t;   
  typedef mesh_<dim>                      mesh_t;
  typedef loc_<dim>                       loc_t;
  typedef elt_<dim>                       elt_t;
  typedef mat<dim_loc_x,dim_loc_y, Cplx>  mat_t;
  typedef typename quadBEM<dim>::qp_t     qp_t;
  typedef typename jac_<dim>::type        jac_t;  
  typedef const quadBEM<dim>              quad_t;
  typedef space_x                         phix_t;
  typedef space_y                         phiy_t;
  
 private:
  
  const vector<loc_t>&   loc;  
  const vect<elt_t>&     elt;
  const normal_t&        nx; 
  const normal_t&        ny;
  const mesh_t&          meshx;
  const mesh_t&          meshy;

  const Real    k2;
  quad_t        qr;
  kernel_t      kernel;
  phix_t        phix;
  phiy_t        phiy;
  mat_t         Melt;
  
  //=========================//
  // Variables intermediaires
  
  elt_t   x,y;       // sommets des elts permutes    
  jac_t   dx,dy;     // jacobien paire elt-ref -> paire elt permutes
  Real    h;         // 2^d * volume paire elt
  int     jx,jy;     // numero des elts (no. local au maillage)
  int     rule;      // regle de quadrature (gestion singularite)
  R3      x_y,x0_y0; // x-y et x0-y0 
  N3      px,py;     // permutations des indices dans les triangles
  Cplx    z;
    
  //====================================//
  //  Choix de la regle de quadrature
  
  void choose_quad(const elt_t& ex, const elt_t& ey){    
    rule = 0; Melt=0.; x=ex; y=ey;
    px[0]=0; px[1]=1; px[2]=2;
    py[0]=0; py[1]=1; py[2]=2;    
    for(int p=0; p<dim+1; p++){
      for(int q=rule; q<dim+1; q++){
	if( &x[p]==&y[q] ){
	  swap(x,rule,p); swap(px,rule,p); 
	  swap(y,rule,q); swap(py,rule,q); 
	  rule++; break;
	}
      }
    }
    
  }
  
 public:  
  
  //=========================//
  //      Constructeur 
 bem(const Real& k, const normal_t& nx0, const normal_t& ny0, int order=15): 
  k2(k*k), kernel(k), qr(order), loc(get_loc_<dim>::apply()), elt(get_elt_<dim>::apply()),
    nx(nx0), meshx(mesh_of(nx0)), ny(ny0), meshy(mesh_of(ny0)) {
    phix.attach_to(mesh_of(nx0)); phiy.attach_to(mesh_of(ny0)); };  
  
  //=====================================//
  // Calcul des interactions elementaires
  const mat_t& operator()(const elt_t& ex, const elt_t& ey){
    
    // Calcul de la regle de quadrature
    choose_quad(ex,ey);
    h     = det_jac(x)*det_jac(y);
    x0_y0 = x[0]-y[0];      
    dx    = mat_jac(x);
    dy    = mat_jac(y);
    const vector<qp_t>& s = qr.x(rule);
    const vector<qp_t>& t = qr.y(rule);
    const vector<Real>& w = qr.w(rule);
    
    // numeros locaux des triangles
    jx    = loc[ &ex-&elt[0] ][meshx];
    jy    = loc[ &ey-&elt[0] ][meshy];
    
    // Boucle sur les points de quadrature    
    for(int j=0; j<s.size(); j++) {
      
      x_y = x0_y0 + dx*s[j] - dy*t[j];      
      z = h*w[j]*kernel.ker(nx[jx],ny[jy],x_y);
      
      for(int kx=0; kx<dim_loc_x; kx++){
	for(int ky=0; ky<dim_loc_y; ky++){
	  Melt(px[kx],py[ky]) += kernel(phix, s[j], jx, kx,
					phiy, t[j], jy, ky,
					nx[jx], ny[jy], x_y,
					h, w[j], z, px, py);
	  
	}
      }
    }

    return Melt;
    
  }
  
  
};




#endif
