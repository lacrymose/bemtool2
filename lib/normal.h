#ifndef NORMAL_H
#define NORMAL_H

#include <vector>
#include <queue>
#include <cassert>
#include "mesh.h"


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
  const vector<loc_t>&      loc;
  
  //_______________
  // Donnee membres
  vector<R3>                nor;
  static const R3           none;
  
  
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

  
 private:
  
  void orienting(const int&, const int&,
		 const int&, const int&,
		 vector<int>&);
  
};

template<int dim>
const R3 normal<dim>::none = 0.;

//===========================//
// Orientation des elements
// de proche en proche

template<>
void normal<dim1>::orienting(const int& j0, const int& k0,
			     const int& j1, const int& k1, vector<int>& aux){
  
  const elt_1D& e = mesh[j1];

  nor[ j1 ][0] = - e[1][1] + e[0][1];
  nor[ j1 ][1] =   e[1][0] - e[0][0];
  nor[ j1 ][2] =   0.;
  normalize(nor[j1]);
  aux [ j1 ] = 1;
  
  if( ( aux[ j0 ]== 1 && k0==k1 ) ||
      ( aux[ j0 ]==-1 && k0!=k1 ) ){
    
    aux[ j1 ] = -1;
    nor[ j1 ] = (-1.)*nor[ j1 ];
    
  }
  
  
}


template<>
void normal<dim2>::orienting(const int& j0, const int& k0,
			     const int& j1, const int& k1, vector<int>& aux){
  N2 I;
  
  I[0] = (k0+1)%3; I[1] = (k0+2)%3;
  elt_1D f0 = mesh[ j0 ][I];
  
  I[0] = (k1+1)%3; I[1] = (k1+2)%3;
  elt_1D f1 = mesh[ j1 ][I];
  
  const elt_2D& e = mesh[j1];
  nor[j1] = vprod( e[1]-e[0], e[2]-e[0]);
  normalize(nor[j1]);
  aux [ j1 ] = 1;  
  
  if( ( aux[ j0 ]== 1 && f0==f1 ) ||
      ( aux[ j0 ]==-1 && f0!=f1 ) ){
    
    aux[ j1 ] = -1;
    nor[ j1 ] = (-1.)*nor[ j1 ];	
    
  }
  
}







//===========================//
// Constructeur de la normale
// par Breadth First Search
template<int dim> normal<dim>::normal(const mesh_t& m):
mesh(m), elt(get_elt_<dim>::apply()), loc(get_loc_<dim>::apply()) {
  
  typedef elt_<dim-1>      face_t;
  typedef array<dim+1,N2>  num_t;
    
  //===============================//
  //  Calcul du graphe d'adjacence
  //  entre elements du maillage
  //===============================//
  
  const R3& n0 = get_node(0);
  int nbnode = nb_node();
  int nbelt  = nb_elt(m);
  int end = -1; 
  
  vector<int>      first(nbnode,end);
  vector<int>      next;
  vector<face_t>   face;
  vector<num_t>    adj(nbelt);
  vector<N2>       num;  
  
  for(int j=0; j<nbelt; j++){
    bool exist;
    array< dim+1 ,face_t> aux = faces_of(m[j]);    
    
    for(int k=0; k<dim+1; k++){
      face_t f = aux[k]; exist = false;
      int& head = first[&f[0]-&n0];
      
      //==================================//
      for(int q = head; q!=end; q=next[q]){
	if(f==face[q]){
	  exist=true;
	  adj[j][k][0] = num[q][0];
	  adj[j][k][1] = num[q][1];
	  adj[ num[q][0] ][ num[q][1] ][0] = j;
	  adj[ num[q][0] ][ num[q][1] ][1] = k;
	  break;}
      }
      
      //==================================//
      if(!exist){
	int q = face.size();
	face.push_back(f);
	next.push_back(head);
	head = q;
	N2 J; J[0] = j; J[1] = k;
	num.push_back(J);
      }
      
    }
    
    
  }
  
  
  //===============================//
  //     Breadth First Search      //
  //===============================//
  
  nor.resize(nbelt,R3());  
  vector<int>      aux;  
  aux.resize(nbelt,0);
  
  int nb_visited = 0;
  queue<int> visit;
  vector<bool> visited(nbelt,false);
  
  // Initialisation de l'algo
  int j0 = 0;
  visit.push(j0);
  orienting(j0+1,0,j0,0,aux);
  visited[j0]=true;
  
  // Lancement de l'algo
  while(nb_visited<nbelt){
    
    // Reinitialisation dans le cas
    // de plusieurs composantes connexes
    if(visit.empty()){
      while(visited[j0]){j0++;}
      orienting(j0+1,0,j0,0,aux); }
    
    else{ j0 = visit.front(); visit.pop(); }
    nb_visited++;
    
    // Boucle sur les voisins de
    // l'element courant
    for(int k0=0; k0<dim+1; k0++){
      
      const int& j1 = adj[j0][k0][0];
      const int& k1 = adj[j0][k0][1];      
      
      // Si voisin pas deja visite:
      // propagation de l'algo
      if(!visited[j1]){
		
	visited[j1]=true;	
	visit.push(j1);
	orienting(j0,k0,j1,k1,aux);
	
      }
      
    }
    
  }

}


typedef normal<1> nrml_1D;
typedef normal<2> nrml_2D;





#endif
