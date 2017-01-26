#ifndef ELT_H
#define ELT_H

#include <cassert>
#include "calculus.h"

namespace bemtool{



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

 friend std::ostream& operator<<(std::ostream& os, const this_t& e){
   for(int j=0; j<dim+1; j++){os << e[j] << std::endl;} return os;}

 friend bool operator==(const this_t& l_, const this_t& r_){
   for(int j=0; j<dim+1; j++){if( &l_[j] != &r_[j]  ){return false;}}
   return true;}

 friend void swap(this_t& e, const int& j, const int& k){
   const v_t* p = &e[j]; e.v_[j] = e.v_[k]; e.v_[k] = p;}
 
 friend int get_dim(const this_t& elt){return dim;}
 
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

//______________________
// Declarations de type
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
inline Real det_jac(const elt_3D& e){return std::abs( (vprod(e[1]-e[0],e[2]-e[0]), e[3]-e[0]) );}

//___________________________
// Barycentre de l'elt de ref
inline R3 center(const elt_1D& e){return (1./2.)*(e[0]+e[1]);}
inline R3 center(const elt_2D& e){return (1./3.)*(e[0]+e[1]+e[2]);}
inline R3 center(const elt_3D& e){return (1./4.)*(e[0]+e[1]+e[2]+e[3]);}

//______________________________________________________
// Comparaison de l'orientation de deux elements VOISINS
// 1 = orientation identique
// 0 = orientation differente
inline bool comp(const elt_1D& e0, const elt_1D& e1){
  if( e0 == e1 ){ return true;  }
  if( (&e0[0]==&e1[0]) || (&e0[1]==&e1[1]) ){ return false; }
  if( (&e0[0]==&e1[1]) || (&e0[1]==&e1[0]) ){ return true;  }
  std::cout << "\nelt.h: comparaison d'elements non voisins" << std::endl;
  std::abort();
}

inline bool comp(const elt_2D& e0, const elt_2D& e1){
  if( e0 == e1 ){ return true;  }
  int jj[2],kk[2], n=0;
  for(int j=0; j<3; j++){
    for(int k=0; k<3; k++){
      if( &e0[j]==&e1[k] ){
	jj[n]=j; kk[n]=k; n++;}
    }
  }
  if(n!=2){std::cout << "\nelt.h: comparaison d'elements non voisins" << std::endl; abort();}
  if( (3+jj[1]-jj[0])%3 != (3+kk[1]-kk[0])%3 ){ return true;}
  return false;

}


//=====================================//
// Calcul du vecteur normal a un element
//=====================================//

inline R3 normal_to(const elt_1D& e){
  R3 e2; e2[2]=1.; return vprod(e2, e[1]-e[0]);}

inline R3 normal_to(const elt_2D& e){
  return vprod(e[1]-e[0],e[2]-e[1]);}


//=============================//
// Calcul des vecteurs normaux
// aux faces d'un element
//=============================//

template <int dim>
inline void compute_normal_to_faces(const elt_<dim>& e, array<dim+1,R3>& n_){};


template<>
inline void  compute_normal_to_faces<dim1>(const elt_1D& e, array<2,R3>& n_){
  n_[0] = e[1]-e[0];
  normalize(n_[0]);
  n_[1] = (-1.)*n_[0];
}


template<>
inline void  compute_normal_to_faces<dim2>(const elt_2D& e, array<3,R3>& n_){

  R3 u1,u2;
  for(int k=0; k<3; k++){

    N2 I;
    I[0]=(k+1)%3;
    I[1]=(k+2)%3;
    elt_1D f = e[I];

    u1 = f[1]-f[0];
    normalize(u1);

    n_[k] = f[0]-e[k];
    n_[k] = n_[k] - (n_[k],u1)*u1;
    normalize(n_[k]);

  }

}


template<>
inline void  compute_normal_to_faces<dim3>(const elt_3D& e, array<4,R3>& n_){

  R3 u1,u2;
  for(int k=0; k<4; k++){

    N3 I;
    I[0]=(k+1)%4;
    I[1]=(k+2)%4;
    I[2]=(k+3)%4;
    elt_2D f = e[I];

    u1 = f[1]-f[0];
    normalize(u1);
    u2 = f[2]-f[0];
    u2 = u2 - (u2,u1)*u1;
    normalize(u2);

    n_[k] = f[0]-e[k];
    n_[k] = n_[k] -(n_[k],u1)*u1 -(n_[k],u2)*u2;
    normalize(n_[k]);

  }

}

//===========================//
// Routine auxiliaire
// calculer la puissance
// d'une surface/courbe
//===========================//

inline Real solid_angle(const R3& p, const elt_1D& e){
  R2x2 M;
  for(int j=0; j<2; j++){
    for(int k=0;k<2;k++){
      M(j,k) = e[k][j] - p[j];
    }
  }
  return det(M);
}

inline Real solid_angle(const R3& p, const elt_2D& e){
  R3x3 M;
  for(int j=0; j<3; j++){
    for(int k=0;k<3;k++){
      M(j,k) = e[k][j] - p[j];
    }
  }
  return (-1.)*det(M);
}


}
#endif
