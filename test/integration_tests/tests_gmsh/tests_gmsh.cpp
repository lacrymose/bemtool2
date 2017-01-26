#include <bemtool2/mesh.h>
#include <bemtool2/gmsh_calls.h>
#include <math.h>

using namespace bemtool;
////=============================================================================////
////========================= plot_2D_data ======================================////
////=============================================================================////
  
void plot_2D_data(std::string filename) {


  ////================ Construction d'un maillage =================////
	
  gmsh_disc("disc",1.,0.05);

  
  ////================ Chargement d'un maillage ===================////
  
  geometry geom;
  load_node_gmsh(geom,"disc");

  mesh_2D Omega(geom);
  load_elt_gmsh(Omega,0);

  ////================ Creation d'une fonction ====================////
  
  const vect<R3>& node = get_node(geom);
  int nb_nodes  = size(node);
  vect<Real> f; resize(f,nb_nodes);  fill(f,0.);
  for (int i=0;i<nb_nodes;i++){
    f[i] = sin(sqrt(pow(node[i][0],2)+pow(node[i][1],2)+pow(node[i][2],2))*M_PI*3);
  }

  ////=================== Creation d'un plot ======================////
  
  write_gmsh(Omega,f,filename);
}


////=============================================================================////
////========================= plot_3D_data ======================================////
////=============================================================================////

void plot_3D_data(std::string filename) {


  ////================ Construction d'un maillage =================////

  gmsh_sphere("sphere",1.,0.05);


  ////================ Chargement d'un maillage ===================////

  geometry geom;
  load_node_gmsh(geom,"sphere");

  mesh_2D Omega(geom);
  load_elt_gmsh(Omega,0);
  gmsh_clean("sphere");

  ////================ Creation d'une fonction ====================////

  const vect<R3>& node = get_node(geom);
  int nb_nodes  = size(node);
  vect<Real> f; resize(f,nb_nodes);  fill(f,0.);
  for (int i=0;i<nb_nodes;i++){
    f[i] = sin(node[i][1]*M_PI*3);
  }


  ////=================== Creation d'un plot ======================////

  write_gmsh(Omega,f,filename);
}


////=============================================================================////
////============================= Main ==========================================////
////=============================================================================////

int main(int argc, char const *argv[]) {
  plot_2D_data("2D_plot_test");
  plot_3D_data("3D_plot_test");

  return 0;
}
