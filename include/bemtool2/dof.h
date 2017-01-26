#ifndef DOF_H
#define DOF_H

#include <vector>
#include "mesh.h"


namespace bemtool {




// //#############################//
// //   Raviart-Thomas ordre 0    //
// //#############################//
// 
// class RT0{
//   
//  public:
// 
//   static const int dim_loc = 3;
// 
//   //______________
//   // Alias de type
//   typedef mesh_2D        mesh_t;
//   typedef elt_2D         elt_t;
//   typedef loc_2D         loc_t;
//   
//  private:
//   //_______________________
//   // Instance pre-existante
//   const mesh_t*          mesh;  
//   const vect<elt_t>&     elt;
//   const std::vector<loc_t>&   loc;
//   
//   //________________
//   // Donnees membres
//   vect<R3>         div;    // divergences de fct de forme
//   vect<N3>         dof;    // numerotation des dofs  
//   vect<R3x2>       jac;    // jacobien elt-ref -> elt-courant
//   vect<_3xR3>      nor;    // normale aux aretes
//   int              nb_dof; // nombre degres liberte
//   _3xR2            b;      // sommets de l'elt-ref
//   static const N3  none;
//   
//   
//  public:
//   
//   RT0();
//   
//   void attach_to(const mesh_t&);
//   
//   R3 operator()(const R2& x, const int& j, const int& k){
//     return div[j][k]*( jac[j]*(x-b[k]) );}
//   
//   friend const Real& div(const RT0& phi, const int& j, const int& k){
//     return phi.div[j][k];}
//   
//   //__________
//   // Acces ddl
//   const N3& operator[](const int& j){ return dof[j];}
//   
//   const N3& operator[](const elt_2D& e){
//     const int& j = loc[ &e-&elt[0] ][*mesh];
//     if(j==-1){return none;}
//     return dof[j];}
//   
//   friend const int& nb_dof(const RT0& phi){return phi.nb_dof;}
//   
//   //_______________
//   // Calcul de flux
//   template <class fct> Cplx proj(const fct&, const int&, const int&);    
//   
// };
// 
// const N3 RT0::none = -1;
// 
// RT0::RT0(): elt(get_elt_<2>::apply()), loc(get_loc_<2>::apply()) { };
// 
// void RT0::attach_to(const mesh_2D& m){
//   
//   mesh = &m;
//   
//   b[1][0] = 1.;
//   b[2][1] = 1.;
//   
//   int nbn = nb_node();
//   int nbt = nb_elt(m);
//   int end = -1;
//   const R3& n0 = get_node(0);
//   
//   resize(dof,nbt);
//   resize(jac,nbt);
//   resize(div,nbt);
//   resize(nor,nbt);
//   
//   std::vector<int>     first(nbn,-1);
//   std::vector<int>     next;
//   std::vector<elt_1D>  edge;
//   
//   for(int j=0; j<nbt; j++){
//     const elt_2D& t = (*mesh)[j];
//     
//     R3 u1 = t[1]-t[0], u2 = t[2]-t[0]; 
//     jac[j] = (1./norm2(vprod(u1,u2)))*mat_(u1,u2);
//     compute_normal_to_faces<2>( (*mesh)[j], nor[j]);
//     
//     for(int k=0; k<3; k++){      
//       N2 I; I[0] = (k+1)%3; I[1] = (k+2)%3;
//       elt_1D e = t[I]; order(e);
//       int& head = first[ &e[0]-&n0 ];
//       
//       //================================//
//       bool exist = false;
//       for(int q=head; q!=end; q=next[q]){
// 	if(e==edge[q]){
// 	  exist=true;
// 	  dof[j][k]=q;
// 	  div[j][k]=-1.;
// 	  break;}
//       }
//       
//       //=================//
//       if(!exist){
// 	int q = edge.size();
// 	dof[j][k]=q;
// 	div[j][k]=+1.;
// 	edge.push_back(e);
// 	next.push_back(head);
// 	head = q;
//       }
//       
//       //=================//
//       nor[j][k] = div[j][k]*nor[j][k];
//             
//     }
//   }
//   
//   // 1 edge = 1 dof
//   nb_dof = edge.size();
//   
// }
// 
// 
// template <class fct>
// Cplx RT0::proj(const fct& f, const int& j, const int& k){    
//   const R3& a = (*mesh)[j][ (k+1)%3 ];
//   const R3& b = (*mesh)[j][ (k+2)%3 ];
//   R3 m = 0.5*(a+b);
//   // Quadrature sur une arete par une regle de Simpson
//   return (norm2(b-a)/6.)*( nor[j][k],(f(a) + 4.*f(m) + f(b) ) );
// }




//#############################//
//    Elements P1 Lagrange     //
//#############################//

template <int dim>
class P1_{
  
 public:
  static const int dim_loc = dim+1;
  //______________
  // Alias de type  
  typedef mesh_<dim>         mesh_t;
  typedef elt_<dim>          elt_t;
  typedef loc_<dim>          loc_t;
  typedef bemtool::array<dim,Real>    Rd;
  typedef bemtool::array<dim+1,int>   Nloc;
  typedef mat<3,dim,Real>    R3xd;
  typedef bemtool::array<dim+1,R3>    locxR3;
  typedef bemtool::array<dim+1,Rd>    locxRd;  
  
 private:
  //_______________________
  // Instance pre-existante
  const mesh_t*          mesh;  
  const vect<elt_t>&     elt;
  const std::vector<loc_t>&   loc;
  
  //_______________
  // Donnees membres
  vect<locxR3>        nabla;   // gradient des fct de forme
  vect<Nloc>          dof;     // numerotation de dofs
  std::vector< std::vector<std::pair<const elt_t*,int> > >  elts_of_dofs;// liste des triangles qui contiennent le support du dof ainsi que leur indice local
  locxRd              b;       // sommets de l'elt-ref  
  int                 nb_dof;  // nombre degres liberte
  static const Nloc   none; 
  
  //_______________
  // Données auxilaires 
  get_elt_<dim> temp_elt;
  get_loc_<dim> temp_loc;
  
  P1_(const P1_&); 
  
  
 public:
  
//   P1_();

  P1_(const mesh_t&);
  
  void attach_to(const mesh_t&);
  
  Real operator()(const Rd& x, const int& j, const int& k) const {
    if(k==0){return  1.-(Rd(1.),x);} return x[k-1];}
  
  friend const R3& grad(const P1_<dim>& phi, const int& j, const int& k) {
    return phi.nabla[j][k]; }
  
  //__________
  // Acces ddl
  const Nloc& operator[](const int& j) const { return dof[j];}
  
  const Nloc& operator[](const elt_t& e) const {
    const int& j = loc[ &e-&elt[0] ][*mesh];
    if(j==-1){return none;}
    return dof[j];}
  
  friend const int& nb_dof(const P1_<dim>& phi) {return phi.nb_dof;}
  
  friend const std::vector<std::pair<const elt_t*,int> >& get_elts_of_dof(const P1_<dim>& phi, const int& j) {return phi.elts_of_dofs[j];}
  
  //______________
  // Valeur nodale
  template <class fct> Cplx proj(const fct&, const int&, const int&);    
  
};

template <int dim>
const typename P1_<dim>::Nloc P1_<dim>::none = -1;

// template <int dim> P1_<dim>::P1_(const geometry& geom_):
// elt(temp_elt.apply(geom_)), loc(temp_loc.apply(geom_)),temp_elt(), temp_loc() {};

 template <int dim> P1_<dim>::P1_(const mesh_t& m):
  elt(temp_elt.apply(get_geometry(m))), loc(temp_loc.apply(get_geometry(m))), temp_elt(), temp_loc() {
   mesh = &m;
   
   
   
   for(int j=0; j<dim; j++){
     b[j+1][j]=1.;}
   
   int nbnod    = nb_node(get_geometry(m));
   int nbelt    = nb_elt(m);
   const R3& n0 = get_node(get_geometry(m),0);
   
   elts_of_dofs.resize(nbnod);
   resize(dof,nbelt);
   resize(nabla,nbelt);
   nb_dof = 0;
   
   std::vector<int>  num(nbnod,-1);
   for(int j=0; j<nbelt; j++){
     
     compute_normal_to_faces<dim>( (*mesh)[j],nabla[j]);
     for(int k=0; k<dim+1; k++){
       
       const R3& ak = (*mesh)[j][ k ];
       const R3& aj = (*mesh)[j][ (k+1)%(dim+1) ];
       
       const int vert = &ak - &n0;
       if(num[vert]==-1){ num[vert]=nb_dof; nb_dof++;}
       dof[j][k] = num[vert];
       elts_of_dofs[num[vert]].push_back(std::pair<const elt_t*,int> (&((*mesh)[j]),k));
       nabla[j][k] = ( 1./(ak-aj,nabla[j][k]) )*nabla[j][k];
     }
     
   }

 }

//template <int dim>
//void P1_<dim>::attach_to(const mesh_t& m){ 
//  mesh = &m;
//  
//  for(int j=0; j<dim; j++){
//    b[j+1][j]=1.;}
//  
//  int nbnod    = nb_node(get_geometry(m));
//  int nbelt    = nb_elt(m);
//  const R3& n0 = get_node(get_geometry(m),0);
//
//  resize(dof,nbelt);
//  resize(nabla,nbelt);
//  nb_dof = 0;
//  
//  std::vector<int>  num(nbnod,-1);
//  for(int j=0; j<nbelt; j++){
//
//    compute_normal_to_faces<dim>( (*mesh)[j],nabla[j]);
//    for(int k=0; k<dim+1; k++){
//
//      const R3& ak = (*mesh)[j][ k ];
//      const R3& aj = (*mesh)[j][ (k+1)%(dim+1) ];
//      
//      const int vert = &ak - &n0;
//      if(num[vert]==-1){ num[vert]=nb_dof; nb_dof++;}
//      dof[j][k] = num[vert];      
//      nabla[j][k] = ( 1./(ak-aj,nabla[j][k]) )*nabla[j][k];
//    }
//    
//  }
//}

template <int dim> template <class fct>
  Cplx P1_<dim>::proj(const fct& f, const int& j, const int& k){    
  return f( (*mesh)[j][k]); }

typedef P1_<1> P1_1D;
typedef P1_<2> P1_2D;

}
#endif
