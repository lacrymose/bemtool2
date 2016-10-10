#include <iostream>
#include <fstream>
#include "lib/tools.h"

using namespace std;

class planewave{

private:
  Cplx kappa;
  R3   dir;
  
public:

  planewave(const Cplx& k, const R3& d): kappa(k), dir(d) {};  

  Cplx operator()(const R3& x){
    return exp(iu*kappa*(dir,x));}
  
  C3 operator()(const elt_2D& t){
    C3 U;
    U[0] = exp(iu*kappa*(dir,t[0]));
    U[1] = exp(iu*kappa*(dir,t[1]));
    U[2] = exp(iu*kappa*(dir,t[2]));
    return U;}
    
};


int main(){
  
  //#########################//
  Real kappa = 1.;
  load_node("mesh/double-sphere.msh");
  const vect<R3>& node = get_node();  
  int nb_node = size(node);
  
  //#########################//  
  mesh_2D Gamma[3];
  load(Gamma[0],0);
  load(Gamma[1],1);  
  load(Gamma[2],2);    
  
  mesh_2D Omega;  
  Omega += Gamma[0];
  Omega += Gamma[1];

  /*  
  //#########################//
  Real kappa = 1.;
  load_node("mesh/sphere.msh");
  const vect<R3>& node = get_node();  
  int nb_node = size(node);
  
  //#########################//  
  mesh_2D Omega;  
  load(Omega,1);
  */  
  
  int nbtri = nb_elt(Omega);
  P1_2D dof; dof.attach_to(Omega);
  int nbdof = nb_dof(dof);
  
  nrml_2D n_(Omega); swap(n_);
  
  //#########################//
  gmm_dense A(2*nbdof,2*nbdof),M(2*nbdof,2*nbdof);
  bem<P1_2D,P1_2D, SLP_3D>   Vop(kappa,n_,n_);
  bem<P1_2D,P1_2D, DLP_3D>   Kop(kappa,n_,n_);
  bem<P1_2D,P1_2D, TDLP_3D>  TKop(kappa,n_,n_);
  bem<P1_2D,P1_2D, HSP_3D>   Wop(kappa,n_,n_);

  progress bar("assembly", nbtri*nbtri);
  for(int j=0; j<nbtri; j++){
    const elt_2D& tj = Omega[j];
    const N3&     jj = dof[j];
    
    for(int k=0; k<nbtri; k++,bar++){
      const elt_2D& tk = Omega[k];
      const N3&     kk = dof[k];
      
      A(      jj,      kk) += Kop (tj,tk);      
      A(      jj,nbdof+kk) += Vop (tj,tk);            
      A(nbdof+jj,      kk) += Wop (tj,tk);
      A(nbdof+jj,nbdof+kk) += TKop(tj,tk);
      
    }
    
    A(jj,jj)             += (-0.5)*MassP1(tj);
    A(nbdof+jj,nbdof+jj) += (-0.5)*MassP1(tj);
    
    M(jj,jj)             += MassP1(tj);
    M(nbdof+jj,nbdof+jj) += MassP1(tj);
  }
  bar.end();

  
  //#########################//  
  R3 dir; dir[0] = 1.; planewave uinc(kappa,dir);
  vect<Cplx> U; resize(U,2*nbdof); fill(U,0.);
  for(int j=0; j<nbtri; j++){
    
    const elt_2D& t = Omega[j];
    const N3&     I = dof[j];    
    U[I] = uinc(t);
    
    C3 V;
    R3 G = center(t);
    if(abs(G[0])>1e-7){
      V[0] = iu*kappa*(dir,t[0])*uinc(t[0]);
      V[1] = iu*kappa*(dir,t[1])*uinc(t[1]);
      V[2] = iu*kappa*(dir,t[2])*uinc(t[2]);
    } else{V = iu*kappa*dir[0]*uinc(t);}
    U[nbdof+I] = V;
    
  }
  
  vect<Cplx> V; resize(V,2*nbdof); fill(V,0.);  
  mv_prod(V,A,U); Cplx sum1 = 0.;
  vect<Cplx> W; resize(W,2*nbdof); fill(W,0.);   
  mv_prod(W,M,U); Cplx sum2 = 0.;

  for(int j=0; j<2*nbdof; j++){
    sum1 += V[j]*conj(U[j]);
    sum2 += W[j]*conj(U[j]); }
  cout << "sum1/sum2 :\t" << sqrt(abs(sum1/sum2)) << endl;

  //###########################//
  fill(W,0.); cg_solve(M,W,V);
  write(Omega,"mesh/visu.mesh");
  ofstream file; file.open("mesh/visu.bb");
  file << "3	 1\t" << 2*nbdof << "\t 2" << endl;
  int count=0;
  for(int j=0; j<nbdof; j++){
    file << W[j].real() << "\t";
    count++; if(!count%50){
      count=0; file << endl;}
  }
  file.close();
  
}
