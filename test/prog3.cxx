#include <iostream>
#include <fstream>
#include "bemtool2/tools.h"
#include "bemtool2/algo.h"

using namespace std;

int main(){

  load_node("mesh/circle2.msh");  
  mesh_1D Gamma; load(Gamma,0);
  int nbelt = nb_elt(Gamma);  
  nrml_1D n_(Gamma);    
  
  P1_1D dof; dof.attach_to(Gamma);
  int nbdof = nb_dof(dof);
  gmm_dense A(nbdof,nbdof);
  
  Real kappa = 1.;
  bem<P1_1D,P1_1D, SLP_2D>   Vop(kappa,n_,n_);

  for(int j=0; j<nbelt; j++){
    const elt_1D& tj = Gamma[j];
    const N2&     jj = dof[j];
    
    for(int k=0; k<nbelt; k++){
      const elt_1D& tk = Gamma[k];
      const N2&     kk = dof[k];
      
      A(jj,kk) += Vop (tj,tk);      
      
    }
  }


  int n1=1, n2=1;
  vect<Cplx> Un,Vk;
  resize(Un,nbdof);
  resize(Vk,nbdof);  
  for(int j=0; j<nbelt; j++){
    const elt_1D& tj = Gamma[j];
    const N2&     jj = dof[j];
    
    Un[ jj[0] ] = pow(tj[0][0]+iu*tj[0][1], n1);
    Vk[ jj[0] ] = pow(tj[0][0]+iu*tj[0][1], n2);

    Un[ jj[1] ] = pow(tj[1][0]+iu*tj[1][1], n1);
    Vk[ jj[1] ] = pow(tj[1][0]+iu*tj[1][1], n2);
  }

  vect<Cplx> W;
  resize(W,nbdof);
  mv_prod(W,A,Un);
  Cplx sum = 0.;
  for(int j=0; j<nbdof; j++){
    sum += conj(Vk[j])*W[j]; }
  cout << "sum:\t " << sum << endl; 
  

  
  
  
  
}



