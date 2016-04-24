#include <iostream>
#include <fstream>
#include "lib/tools.h"
#include "lib/algo.h"

using namespace std;

int main(){
  
  //#########################//
  Real kappa = 1.;
  load_node("mesh/carre2.msh");
  const vect<R3>& node = get_node();  
  int nb_node = size(node);
  
  //#########################//  
  mesh_2D Omega;
  load(Omega,0);
   
  int nbtri = nb_elt(Omega);
  P1_2D dof; dof.attach_to(Omega);
  int nbdof = nb_dof(dof);

  nrml_2D n_(Omega); swap(n_);
  adjacency2D adj(Omega);
  
  

  
  
  

  
  
  
}
