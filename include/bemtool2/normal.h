#ifndef NORMAL_H
#define NORMAL_H

#include <vector>
#include <queue>
#include <cassert>
#include "mesh.h"
#include "algo.h"

//==========================//
//       Normale            //
//==========================//

template<int dim>
class normal{

 public:
  //______________
  // Alias de type
  typedef normal<dim>    this_t;
  typedef mesh_<dim>     mesh_t;
  typedef elt_<dim>      elt_t;
  typedef loc_<dim>      loc_t;
  
 private:
  //_________________________
  // Instances pre-existantes
  const mesh_t&             mesh;
  const vect<elt_t>&        elt;
  const std::vector<loc_t>&      loc;
  
  //_______________
  // Donnee membres
  std::vector<bool>              orientation;
  std::vector<R3>                nor;
  static const R3           none;
  
  //_______________
  // Données auxilaires 
  get_elt_<dim> temp_elt;
  get_loc_<dim> temp_loc;
  
 public:
  
  normal(const mesh_t&);
  
  const R3& operator[](const int& j) const {
    return nor[j];}
  
  const R3& operator[](const elt_t& e) const {
    const int& j = loc[ &e-&elt[0] ][mesh];
    if(j==-1){return none;}
    return nor[j];}
  
  inline friend void swap(this_t& N){
    for(int j=0; j<N.nor.size(); j++){
      N.nor[j] = (-1.)*N.nor[j];} }
  
  inline friend void swap(this_t& N, const int& j){
    N.nor[j] = (-1.)*N.nor[j]; }
  
  inline friend const mesh_t& mesh_of(const this_t& N) {
    return N.mesh;}
  
};

template<int dim>
const R3 normal<dim>::none = 0.;

//===========================//
// Constructeur de la normale
template<int dim> normal<dim>::normal(const mesh_t& m):
mesh(m), elt(temp_elt.apply(get_geometry(m))), loc(temp_loc.apply(get_geometry(m))),temp_elt(), temp_loc() {
  
  int nbelt  = nb_elt(mesh);
  bool ok = true;
  orientation.assign(nbelt,ok);
  std::vector<bool> visited(nbelt,false);
  
  //====================================
  // Calcul de l'adjacence entre elements
  adjacency<mesh_t> adj(mesh);
  
  //====================================
  // Recherche des composantes connexes
  connected<mesh_t> component(mesh);
  int nb_component = nb_(component);
  Real global_orientation[nb_component];

  //====================================
  // initialisation  de la recherche
  // d'un point extremal du maillage
  int  Iext = 0;
  Real Ext = norm2(center( mesh[ component[Iext][0] ] ));

  //===============================//
  //   Breadth First Search sur    //
  //   chaque composante connexe   //
  //===============================//
  for(int I=0; I<nb_component; I++){

    int nbe = component[I].size();
    int nb_visited = 0;
    std::queue<int> visit;

    int j0 = component[I][0];
    visit.push(j0);
    visited[j0]=true;
    
    while(nb_visited<nbe){
      
      j0 = visit.front();
      const elt_t& e0 = mesh[j0];
      visit.pop();
      nb_visited++;
      
      if( norm2(center(e0)) > Ext){
	Ext = norm2(center(e0));
	Iext = I; }
      
      for(int k0=0; k0<dim+1; k0++){
	const int& j1 = adj[j0][k0];
	const elt_t& e1 = mesh[j1];
	
	if(!visited[j1]){
	  bool same = comp(e1,e0);
	  if(same){orientation[j1]=orientation[j0];}
	  else{orientation[j1]=!orientation[j0];}
	  visited[j1]=true;
	  visit.push(j1);
	}
      }
    }
    
    global_orientation[I] = 0.;
    const R3& p = mesh[ component[I][0] ][0];
    for(int j=0; j<nbe; j++){
      j0 = component[I][j];
      const elt_t& e = mesh[j0];
      Real r = solid_angle(p,e);
      if(!orientation[j0]){r = -r;}
      global_orientation[I] += r;
    }
    
    if(global_orientation[I]>0){      
      for(int j=0; j<nbe; j++){
	j0 = component[I][j];
	orientation[j0] = !orientation[j0];
      }
    }
        
  }

  //=====================================
  // Calcul effectif des vecteurs normaux
  nor.resize(nbelt);
  for(int j=0; j<nbelt; j++){
    nor[j]=normal_to(mesh[j]);
    normalize(nor[j]);
    if(orientation[j]){nor[j] = (-1.)*nor[j];}
  }
  
  //=====================================
  // Si le domaine est borne la
  // composante exterieure du bord
  // doit etre orientee dans l'autre sens
  if(mesh.is_bounded()){
    for(int j=0; j<component[Iext].size(); j++){
      int jj = component[Iext][j];
      nor[jj] = (-1.)*nor[jj];
    }
  }
  
}


typedef normal<1> nrml_1D;
typedef normal<2> nrml_2D;



#endif
