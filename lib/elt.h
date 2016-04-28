#ifndef ELT_H
#define ELT_H

#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include "calculus.h"

using namespace std;


//==========================//
//        Element           //
//==========================//

// Liste des indices constituant
// les faces de l'element
template <int dim> array< dim+1, array<dim, int> > face_idx(){
  array< dim+1, array<dim, int> > I;
  for(int j=0; j<dim+1; j++)
    for(int k=0; k<dim; k++){
      I[j][k] = (j+1+k)%(dim+1);}
  return I; }


template <int dim, class val_t = R3>
  class elt_{
  
 public:
 typedef val_t        v_t;
 typedef elt_<dim> this_t;
 
 private:
 const v_t*   v_[dim+1];
 
 public:
 elt_<dim,val_t>(){};
 
 template <class r_t> elt_<dim,val_t>(const r_t& r_){
   for(int j=0; j<dim+1;j++){v_[j] = &r_[j];} }
 
 template <class r_t> void operator=(const r_t& r_){
   for(int j=0; j<dim+1;j++){v_[j] = &r_[j];} }  
 
 const v_t& operator[](const int& j) const {return *v_[j];}  
 
 template <class i_t> subarray<const this_t,i_t> operator[] (const i_t& i_) const {
   return subarray<const this_t,i_t>(*this,i_);}  
 
 friend ostream& operator<<(ostream& os, const this_t& e){
   for(int j=0; j<dim+1; j++){os << e[j] << endl;} return os;}
 
 friend bool operator==(const this_t& l_, const this_t& r_){
   for(int j=0; j<dim+1; j++){if( &l_[j] != &r_[j]  ){return false;}}
   return true;}
 
 friend void swap(this_t& e, const int& j, const int& k){
   const v_t* p = &e[j]; e.v_[j] = e.v_[k]; e.v_[k] = p;}
 
 friend array< dim+1 ,elt_<dim-1> > faces_of(const this_t& e){
   array< dim+1, elt_<dim-1> > f;
   array<dim, int> I;
   for(int k=0; k<dim+1; k++){
     for(int p=0; p<dim; p++){
       I[p] = (k+1+p)%(dim+1);}
     f[k] = e[I]; order(f[k]);}
   return f;
 }
 
 
};

typedef elt_<dim0> elt_0D; // noeud
typedef elt_<dim1> elt_1D; // arrete
typedef elt_<dim2> elt_2D; // triangle
typedef elt_<dim3> elt_3D; // tetrahedre

//____________________________
// Numeros des noeuds des elts
// mis par ordre croissant  
template <int dim>
inline void order(elt_<dim>& e){
  for(int j=1; j<dim+1; j++){
    for(int k=0; k<dim+1-j; k++){
      if(&e[k]-&e[k+1]>0){swap(e,k,k+1);}
    }
  }
  
}


//___________________
// Calculs de volumes
inline Real vol(const elt_1D& e){return norm2(e[1]-e[0]);}
inline Real vol(const elt_2D& e){return 0.5*norm2(vprod(e[1]-e[0],e[2]-e[0])); }
inline Real vol(const elt_3D& e){return (e[3]-e[0],vprod(e[2]-e[0],e[1]-e[0]))/6.;}

//_________________________
// Transfo sur l'elt de ref
template <int dim> struct jac_   {typedef mat<3,dim,Real> type;};
template <>        struct jac_<1>{typedef R3              type;};

inline R3   mat_jac(const elt_1D& e){return e[1]-e[0];}
inline R3x2 mat_jac(const elt_2D& e){return mat_(e[1]-e[0], e[2]-e[0]);}
inline R3x3 mat_jac(const elt_3D& e){return mat_(e[1]-e[0], e[2]-e[0], e[3]-e[0]);}

inline Real det_jac(const elt_1D& e){return norm2(e[1]-e[0]);}
inline Real det_jac(const elt_2D& e){return norm2(vprod(e[1]-e[0],e[2]-e[0]));}
inline Real det_jac(const elt_3D& e){return abs( (vprod(e[1]-e[0],e[2]-e[0]), e[3]-e[0]) );}

//___________________________
// Barycentre de l'elt de ref
inline R3 center(const elt_1D& e){return (1./2.)*(e[0]+e[1]);}
inline R3 center(const elt_2D& e){return (1./3.)*(e[0]+e[1]+e[2]);}
inline R3 center(const elt_3D& e){return (1./4.)*(e[0]+e[1]+e[2]+e[3]);}


#endif
