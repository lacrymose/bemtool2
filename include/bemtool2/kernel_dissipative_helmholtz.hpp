#ifndef KERNEL_DH_HPP
#define KERNEL_DH_HPP


#include <boost/math/special_functions/bessel.hpp>

namespace bemtool{
/*=================================
||  NOYAUX INTEGRAUX DISSIPATIFS ||
=================================*/
//___________________________
// Single Layer Potential 2-D

// class SLP_DH_2D{
//
//   typedef Real qp_t;
//   const   Real k;
//   Real r;
//   Cplx val;
//
//  public:
//   static const int dim = 1;
//
//  SLP_DH_2D(const Real& k0): k(k0) {};
//
//   inline Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
//     r = norm2(x_y); return val = (0.5/pi)*K_0(k*r); }
//
//   inline Cplx& ker(const R3& ny, const R3& x_y){
//     r = norm2(x_y); return val = (0.5/pi)*K_0(k*r); }
//
//   template <class phix_t, class phiy_t>
//     inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
// 			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
// 			    const R3& nx, const R3& ny,   const R3& x_y,
// 			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
//   {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;}
//
//   template <class phi_t>
//     inline Cplx& operator()(const phi_t& phi, const qp_t& t, const int& jy, const int& ky,
// 			    const R3& ny,   const R3& x_y,
// 			    const Real& h, const Real& w)
//   {return val = phi(t,jy,ky)*w*h*ker(ny,x_y);}
//
// };

//___________________________
// Single Layer Potential 3-D

class SLP_DH_3D{

  typedef R2 qp_t;
  const Real k;
  Real r;
  Cplx val;

 public:
  static const int dim = 2;

 SLP_DH_3D(const Real& k0): k(k0) {};

  inline Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); return val = exp(-k*r)/(4*pi*r); }

  inline Cplx& ker(const R3& ny, const R3& x_y){
    r = norm2(x_y); return val = exp(-k*r)/(4*pi*r); }

  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx, const R3& ny,   const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
  {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;}

  template <class phi_t>
    inline Cplx& operator()(const phi_t& phi, const qp_t& t, const int& jy, const int& ky,
			    const R3& ny,   const R3& x_y,
			    const Real& h, const Real& w)
  {return val = phi(t,jy,ky)*w*h*ker(ny,x_y);}

};

//___________________________
// Double Layer Potential 2-D
// class DLP_DH_2D{
//
//   typedef Real qp_t;
//   const Real k;
//   Real r;
//   Cplx val;
//
//  public:
//   static const int dim = 1;
//   DLP_2D_D(const Real& k0): k(k0) {};
//
//   inline Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
//      r = norm2(x_y); return val = (-0.5/pi)*k*K_1(k*r)*(ny,u)*(1/r); }
//
//   inline Cplx& ker(const R3& ny, const R3& x_y){
//      r = norm2(x_y); return val = (-0.5/pi)*k*K_1(k*r)*(ny,u)*(1/r); }
//
//   template <class phix_t, class phiy_t>
//      inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
//            const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
//            const R3& nx,  const R3& ny,  const R3& x_y,
//            const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
//    {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;}
//
//    template <class phi_t>
//      inline Cplx& operator()(const phi_t& phi, const qp_t& t, const int& jy, const int& ky,
//            const R3& ny,  const R3& x_y,
//            const Real& h, const Real& w)
//    {return val = phi(t,jy,ky)*h*w*ker(ny,x_y);}
//
// };

//___________________________
// Double Layer Potential 3-D


class DLP_DH_3D{

  typedef R2 qp_t;
  const Real k;
  Real r,r3;
  Cplx val;

 public:
  static const int dim = 2;

 DLP_DH_3D(const Real& k0): k(k0) {};

  inline Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); r3 = r*r*r;
    return val = (ny,x_y)*(-k*r-1.)*exp(-k*r)/(4*pi*r3); }

  inline Cplx& ker(const R3& ny, const R3& x_y){
    r = norm2(x_y); r3 = r*r*r;
    return val = (ny,x_y)*(-k*r-1.)*exp(-k*r)/(4*pi*r3); }

  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx,  const R3& ny,  const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
  {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;}

  template <class phi_t>
    inline Cplx& operator()(const phi_t& phi, const qp_t& t, const int& jy, const int& ky,
			    const R3& ny,   const R3& x_y,
			    const Real& h, const Real& w)
  {return val = phi(t,jy,ky)*w*h*ker(ny,x_y);}
};

//____________________________________
// Adjoint Double Layer Potential 2-D

// class TDLP_DH_2D{
//
//   typedef Real qp_t;
//   const Real k;
//   Real r;
//   Cplx val;
//
//  public:
//   static const int dim = 1;
//
//  TDLP_2D_D(const Real& k0): k(k0) {};
//
//   inline Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
//     r = norm2(x_y); val = (-0.5/pi)*k*K_1(k*r);
//     val = (nx,u)*(1/r)*val; return val;
//   }
//
//   template <class phix_t, class phiy_t>
//     inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
// 			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
// 			    const R3& nx,  const R3& ny,  const R3& x_y,
// 			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
//   {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;}
//
// };

//____________________________________
// Adjoint Double Layer Potential 3-D
class TDLP_DH_3D{

  typedef R2 qp_t;
  const Real k;
  Real r,r3;
  Cplx val;

 public:
  static const int dim = 2;

 TDLP_DH_3D(const Real& k0): k(k0) {};

  Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); r3 = r*r*r;
    return val = (nx,x_y)*(-k*r-1.)*exp(-k*r)/(4*pi*r3); }

  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx,  const R3& ny,  const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py)
  {return val = phix(s,jx,kx)*phiy(t,jy,ky)*z;}

};


//____________________________________
// Hypersingular Layer Potential 2-D

// class HSP_DH_2D{
//
//   typedef Real qp_t;
//   const Real k,k2;
//   Real r;
//   Cplx val;
//
//  public:
//   static const int dim = 1;
//
//  HSP_2D_D(const Real& k0): k(k0), k2(k0*k0) {};
//
//   Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
//     r = norm2(x_y); return val = 0.25*iu*H_0(k*r);}
//
//   template <class phix_t, class phiy_t>
//     inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
// 			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
// 			    const R3& nx,  const R3& ny,  const R3& x_y,
// 			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py){
//     val = ( vprod( nx, grad(phix, jx, px[kx]) ), vprod( ny, grad(phiy, jy, py[ky]) ) )*z;
//     val+= -k2*(nx,ny)*phix(s,jx,kx)*phiy(t,jy,ky)*z;
//     return val;}
//
// };




//____________________________________
// Hypersingular Layer Potential 3-D
class HSP_DH_3D{

  typedef R2 qp_t;
  const Real k,k2;
  Real r;
  Cplx val;

 public:
  static const int dim = 2;

 HSP_DH_3D(const Real& k0): k(k0), k2(k0*k0) {};

  Cplx& ker(const R3& nx, const R3& ny, const R3& x_y){
    r = norm2(x_y); return val = exp(-k*r)/(4*pi*r);}

  template <class phix_t, class phiy_t>
    inline Cplx& operator()(const phix_t& phix, const qp_t& s, const int& jx, const int& kx,
			    const phiy_t& phiy, const qp_t& t, const int& jy, const int& ky,
			    const R3& nx,  const R3& ny,  const R3& x_y,
			    const Real& h, const Real& w, const Cplx& z, const N3& px, const N3& py){
    val = ( vprod( nx, grad(phix, jx, px[kx]) ), vprod( ny, grad(phiy, jy, py[ky]) ) )*z;
    val+= -k2*(nx,ny)*phix(s,jx,kx)*phiy(t,jy,ky)*z;
    return val;}

};

}

#endif
