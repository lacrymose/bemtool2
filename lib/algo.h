#ifndef ALGO_H
#define ALGO_H

#include <vector>
#include <iostream>
#include <queue>
#include "calculus.h"
#include "mesh.h"


using namespace std;

template <typename m_t>
class adjacency{
  
 public:
  static const int Dim =         m_t::Dim;
  static const int dim =         m_t::Dim;
  typedef adjacency<m_t>         this_t;
  typedef m_t                    mesh_t;
  typedef elt_<dim>              elt_t;
  typedef elt_<dim-1>            face_t;
  typedef array<dim+1,face_t>    face_array_t;
  typedef array<dim+1,int>       num_array_t;
  
 private:
  const mesh_t&           mesh;
  const vect<elt_t>&      elt;
  vector<face_t>          face;
  vector<num_array_t>     neig;
  vector<num_array_t>     back;
  
 public:
  // Constructeur
  adjacency(const mesh_t&);
  
  // Acces aux donnees
  const num_array_t&
    operator[](const int& j){
    return neig[j];}
  
  const num_array_t&
    operator()(const int& j){
    return back[j];}
  
  friend const vector<face_t>&
    faces_of(const this_t& adj){
    return adj.face;}
  
};


template <typename m_t>
adjacency<m_t>::adjacency(const mesh_t&  m):
mesh(m), elt(get_elt_<dim>::apply()) {
  
  const R3& n0 = get_node(0);
  int nbnode = nb_node();
  int nbelt  = nb_elt(m);
  int end    = -1; 
  
  vector<int>  first; 
  vector<int>  next;
  vector<N2>   num;  

  first.resize(nbnode,end);
  neig.resize(nbelt);
  back.resize(nbelt);  
    
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
	  neig[j][k] = num[q][0];
	  back[j][k] = num[q][1];
	  neig[ num[q][0] ][ num[q][1] ] = j;
	  back[ num[q][0] ][ num[q][1] ] = k;
	  break;
	}
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
      //==================================//      
      
    }
  }
  
}

typedef adjacency<mesh_1D> adjacency1D;
typedef adjacency<mesh_2D> adjacency2D;
typedef adjacency<mesh_3D> adjacency3D;




template <typename m_t>
class connected{
  
 public:
  static const int Dim =    m_t::Dim;
  static const int dim =    m_t::Dim;  
  typedef connected<m_t>    this_t;
  typedef m_t               mesh_t;
  typedef elt_<dim>         elt_t;
  typedef elt_<dim-1>       face_t;

 private:
  // tableau  avec les no. des elements
  // de chaque composante.
  // num[j][k] est le no du k ieme elt du
  // la composante no. j
  vector< vector<int> >     num;

  // nbre de composantes
  int                       nbc; 

  // reference vers le maillage 
  const mesh_t&             mesh;
  
 public:
  connected(const mesh_t&);
  
  const vector<int>& operator[](const int& j){
    return num[j];}
  
  friend const int& nb_(const this_t& component){
    return component.nbc;}
  
};

  
//===============================//
//     Breadth First Search      //
//===============================//
template <typename m_t>
connected<m_t>::connected(const mesh_t& m): mesh(m) {  
  
  nbc = 1;
  num.resize(nbc);
  
  int nbelt  = nb_elt(mesh);
  adjacency<mesh_t> adj(mesh);

  int nb_visited = 0;
  queue<int> visit;
  vector<bool> visited(nbelt,false);
  
  // Initialisation de l'algo
  int j0 = 0;

  visit.push(j0);
  visited[j0]=true;
  num[nbc-1].push_back(j0);
  
  // Lancement de l'algo
  while(nb_visited<nbelt){
    
    // Reinitialisation dans le cas
    // de plusieurs composantes connexes
    if(visit.empty()){
      nbc++; num.resize(nbc);
      j0=0; while(visited[j0]){j0++;}
      visit.push(j0);
      visited[j0]=true;
      num[nbc-1].push_back(j0);      
    }
    else{
      j0 = visit.front();
      visit.pop();
    }
    nb_visited++;
    
    // Boucle sur les voisins de
    // l'element courant
    for(int k0=0; k0<dim+1; k0++){
      const int& j1 = adj[j0][k0];
      
      // Si voisin pas deja visite:
      // propagation de l'algo
      if(!visited[j1]){
	visited[j1]=true;	
	visit.push(j1);
	num[nbc-1].push_back(j1);      	
      }
    }
  }
  
  
  
  
  
  
  
  
}

typedef connected<mesh_1D> connected1D;
typedef connected<mesh_2D> connected2D;
typedef connected<mesh_3D> connected3D;



#endif
