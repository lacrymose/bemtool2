#include <iostream>
#include <fstream>
#include "lib/tools.h"
#include "lib/algo.h"

using namespace std;

int main(){

  //#########################//  
  /*
  load_node("mesh/carre2.msh");
  mesh_2D Omega;
  load(Omega,0);
  int nbtri = nb_elt(Omega);
  cout << "nbtri = " << nbtri << endl;
  P1_2D dof; dof.attach_to(Omega);
  int nbdof = nb_dof(dof);
  adjacency2D adj(Omega);
  */
  //#########################//
  /*
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
  */
  //###############################//
  
  load_node("mesh/pave_troue.msh");  
  mesh_2D Gamma[2];  
  load(Gamma[0],0);
  load(Gamma[1],1);

  mesh_2D Omega;  
  Omega+=Gamma[0];
  Omega+=Gamma[1];
  
  connected2D component(Omega);
  cout << "nb component: "<< nb_(component) << endl;
  cout << "nb elt 0:" << component[0].size() << endl;
  cout << "nb elt 1:" << component[1].size() << endl;  
  
}



