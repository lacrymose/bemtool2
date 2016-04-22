#ifndef FEM_P1
#define FEM_P1

#include "mesh.h"


R2x2 StiffP1(const elt_1D& e){
  Real l = vol(e);
  R2x2 K;
  K(0,0) =  1./l;   K(1,1) =  1./l;
  K(0,1) = -1./l;   K(1,0) = -1./l;  
  return K;
}


R3x3 StiffP1(const elt_2D& e){
  R3 U[3];
  U[0] = e[2]-e[1];
  U[1] = e[0]-e[2];
  U[2] = e[1]-e[0];
  Real area = vol(e);
  
  R3x3 K;
  for(int j=0; j<3; j++){
    for(int k=0; k<3; k++){
      K(j,k) = (U[j],U[k])/(4.*area);
    }
  }
  return K;
}


R2x2 MassP1(const elt_1D& e){
  Real l = vol(e);
  R2x2 M;
  M(0,0) = l/3.; M(1,1) = l/3.;
  M(0,1) = l/6.; M(1,0) = l/6.;  
  return M;
}


R3x3 MassP1(const elt_2D& e){
  Real a = vol(e);
  R3x3 M;
  for(int j=0; j<3; j++){
    for(int k=0; k<3; k++){    
      if(j==k){ M(j,k) = a/6. ; }
      else    { M(j,k) = a/12.; }
    }
  }
  return M;
}


R4x4 MassP1(const elt_3D& e){
  Real v = vol(e);
  R4x4 M;
  for(int j=0; j<4; j++){
    for(int k=0; k<4; k++){
      if(j==k){ M(j,k) = v/10.; }
      else    { M(j,k) = v/20.; }
    }
  }
  return M;
}


#endif
