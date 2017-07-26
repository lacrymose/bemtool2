#ifndef KERNEL_L_HPP
#define KERNEL_L_HPP


namespace bemtool{
/*==================================
||       NOYAUX DE LAPLACE        ||
||  rem: le "L" signifie Laplace  ||
==================================*/
//
// // Single Layer Potential
// class SLP_L_2D{
//  private: const Real k; Real r; Cplx val;
//  public:
//  SLP_L(const Real& k0): k(k0) {};
//   Cplx& operator()(const R2& nx, const R2& ny, const R2& u){
//     r = norm2(u); val = (0.5/pi)*log(r); return val;} };
//
// // Double Layer Potential
// class DLP_L_2D{
//  private: const Real k; Real r; Cplx val;
//  public:  DLP_L(const Real& k0): k(k0) {};
//   Cplx& operator()(const R2& nx, const R2& ny, const R2& u){
//     r = norm2(u); val = 0.5/(pi*r);
//     val = (ny,u)*(1/r)*val; return val;} };
//
// // Transpose Double Layer Potential
// class TDLP_L_2D{
//  private: const Real k; Real r; Cplx val;
//  public:  TDLP_L(const Real& k0): k(k0) {};
//   Cplx& operator()(const R2& nx, const R2& ny, const R2& u){
//     r = norm2(u); val = 0.5/(pi*r);
//     val = (nx,u)*(1/r)*val; return val;} };
//
// // Hyper Singular Potential
// class HSP_L_2D{
//  private: const Real k; Real r; Cplx val;
//  public:  HSP_L(const Real& k0): k(k0) {};
//   Cplx& operator()(const R2& nx, const R2& ny, const R2& u){
//     r = norm2(u); val = (0.5/pi)*log(r); return val;} };

}
#endif
