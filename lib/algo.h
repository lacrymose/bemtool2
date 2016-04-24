#ifndef ALGO_H
#define ALGO_H

#include <vector>
#include <iostream>
#include "calculus.h"
#include "mesh.h"

using namespace std;


template <int dim>
class adjacency{
  
 public:
  typedef adjacency<dim>         this_t;
  typedef mesh_<dim>             mesh_t;
  typedef elt_<dim>              elt_t;
  typedef elt_<dim-1>            face_t;
  typedef array<dim+1,face_t>    face_array_t;
  typedef array<dim+1,N2>        num_array_t;
  
 private:
  const mesh_<dim>&       mesh;
  const vect<elt_t>&      elt;
  vector<face_t>          face;
  vector<num_array_t>     neig;

 public:
  // Constructeur
  adjacency(const mesh_t&);
  
  // Acces aux donnees
  const num_array_t&
    operator[](const int& j){
    return neig[j];}
  
  friend const vector<face_t>&
    faces_of(const this_t& adj){
    return adj.face;}
  
};


template <int dim>
adjacency<dim>::adjacency(const mesh_t&  m):
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
	  neig[j][k][0] = num[q][0];
	  neig[j][k][1] = num[q][1];
	  neig[ num[q][0] ][ num[q][1] ][0] = j;
	  neig[ num[q][0] ][ num[q][1] ][1] = k;
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
      //==================================//      
      
    }
  }
  

  
  
}


typedef adjacency<1> adjacency1D;
typedef adjacency<2> adjacency2D;
typedef adjacency<3> adjacency3D;

#endif
