#include <iostream>
#include <fstream>
#include "lib/tools.h"
#include "lib/algo.h"

using namespace std;

int main(){

  //#########################//  
  /*
  array<5,R3> M;
  M[0][0]= 0; M[0][1]=0;
  M[1][0]= 1; M[1][1]=0;
  M[2][0]= 1; M[2][1]=1;
  M[3][0]=-1; M[3][1]=1;

  N2 I0, I1;
  I0[0]=1; I0[1]=2;
  I1[0]=3; I1[1]=2;

  elt_1D e0(M[I0]);
  elt_1D e1(M[I1]);
  cout << normal_to(e0) << endl;   
  */
  
  /*
  cout << "\ne0:\n" << e0;
  cout << "\ne1:\n" << e1;  
  cout << "\ncomp(e0,e1):  "<< comp(e0,e1);
  cout << endl << endl;   
  */
  
  /*
  N3 I0, I1;
  I0[0]=1; I0[1]=2; I0[2]=0;
  I1[0]=0; I1[1]=3; I1[2]=4;
  
  elt_2D e0(M[I0]);
  elt_2D e1(M[I1]);  

  cout << "\ne0:\n" << e0;
  cout << "\ne1:\n" << e1;  
  cout << "\ncomp(e0,e1):  "<< comp(e0,e1);
  cout << endl << endl; 
  */  
  //#########################//  
  /*
  load_node("mesh/carre2.msh");
  mesh_2D Omega;
  load(Omega,0);
  int nbtri = nb_elt(Omega);
  cout << "nbtri = " << nbtri << endl;
  
  nrml_2D n_(Omega);
  for(int j=0; j<nbtri; j++){
    cout << n_[j] << endl;}
  */  
  //#########################//  
  /*
  load_node("mesh/sphere.msh");
  mesh_2D Omega;
  load(Omega,1);
  int nbtri = nb_elt(Omega);
  cout << "nbtri = " << nbtri << endl;
  
  nrml_2D n_(Omega);
  for(int j=0; j<nbtri; j++){
    const elt_2D& t = Omega[j];
    cout << (center(t),n_[j]) << endl;
  }
  */    
  //#########################//
  
  load_node("mesh/rectangle_troue.msh");  
  mesh_1D Gamma[3];  
  load(Gamma[0],0);
  load(Gamma[1],1);
  load(Gamma[2],2);  
  
  mesh_1D Omega;  
  Omega+=Gamma[0];
  Omega+=Gamma[1];
  Omega+=Gamma[2];  

  connected1D component(Omega);
  cout << "nb component: "<< nb_(component) << endl;
  cout << "nb elt 0:" << component[0].size() << endl;
  cout << "nb elt 1:" << component[1].size() << endl;  
  cout << "nb elt 2:" << component[2].size() << endl;    
  
  nrml_1D n_(Omega);  

  R3 p; p[0]=1.; p[1]=1.;
  for(int j=0; j<8; j++){
    int jj = component[1][j];
    R3 u = center(Omega[jj])-p;
    cout << (u,n_[jj]) << endl;
  }
  
  //###############################//
  /*  
  load_node("mesh/pave_troue.msh");  
  mesh_2D Gamma[2];  
  load(Gamma[0],0);
  load(Gamma[1],1);

  mesh_2D Omega;  
  Omega+=Gamma[0];
  Omega+=Gamma[1];
  Omega = unbounded;
  
  
  nrml_2D n_(Omega);  
  
  connected2D component(Omega);
  cout << "nb component: "<< nb_(component) << endl;
  cout << "nb elt 0:" << component[0].size() << endl;
  cout << "nb elt 1:" << component[1].size() << endl;  
  
  R3 p; p[0]=0.5; p[1]=0.5; p[2]=0.5;
  for(int j=0; j<component[1].size(); j++){
    int jj = component[1][j];
    R3 u = center(Omega[jj])-p;
    cout << (u,n_[jj]) << endl;
  }  
  */  
  
  
  
  
  
}



