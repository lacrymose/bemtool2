#ifndef MESH_H
#define MESH_H

#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include "calculus.h"
#include "elt.h"


//==========================//
//   Numerotation locale    //  
//   d'un element           //
//==========================//

// Declarations prealables 
// sur les maillages
template <int dim> class mesh_;
template <class m_t> int nb_elt (const m_t&);
typedef mesh_<dim1> mesh_1D;
typedef mesh_<dim2> mesh_2D;  
typedef mesh_<dim3> mesh_3D;  


template <int dim>
class loc_{
  
 private:
  
  typedef const void*        admesh_t; 
  std::vector<admesh_t>           admesh;
  std::vector<int>                num;
  static const int           none; 
  
 public:
  
  loc_<dim>(const loc_<dim>& nl):
  admesh(nl.admesh), num(nl.num) {};
  
  template <class m_t>
    loc_<dim>(const m_t& m_, const int& num0):
  admesh(1,&m_), num(1,num0) {};
  
  template <class m_t>
    void push(const m_t& m_, const int& num0){
    admesh.push_back(&m_); num.push_back(num0);}
  
  template <class m_t>
    const int& operator[](const m_t& m_) const {
    for(int j=0; j<admesh.size(); j++){
      if(&m_==admesh[j]){return num[j];}}
    return none;}
  
};


template<int dim>
const int loc_<dim>::none = -1;

typedef loc_<dim0>  loc_0D;
typedef loc_<dim1>  loc_1D;
typedef loc_<dim2>  loc_2D;

//==========================//
//    Registre d'elements   //  
//==========================//

template <int dim>
class list_elt_{
  
 public:
  typedef list_elt_<dim>  this_t;
  typedef mesh_<dim>      mesh_t;
  typedef elt_<dim>       elt_t;
  typedef loc_<dim>       loc_t;  
  
 private:  
  // Liste des elements enregistres
  vect<elt_t>       elt;
  
  // Numerotations locales
  std::vector<loc_t>     loc;
  
  // Premier noeud de la scene
  const R3*         ad0;
  
  // Gestion des doublons
  std::vector<int>       first;
  std::vector<int>       next;
  
 public:  
  list_elt_(){};  
  int  push(const mesh_t&, elt_t);
  void init(const int& nb_node, const R3& p0){
    ad0 = &p0;
    assert(first.size()==0);
    first.resize(nb_node,-1);}
  
  friend const int                    size(const this_t& l)    {return size(l.elt);}
  friend const vect  <elt_t>&         get_elt(const this_t& l) {return l.elt;}
  friend const std::vector<loc_t>&    get_loc(const this_t& l) {return l.loc;}  
  
};

// Enregistrement d'un element
template <int dim>
int list_elt_<dim>::push(const mesh_<dim>& m, elt_<dim> e){
    
  order(e);
  int end = -1;
  bool exist = false;
  int& head = first[ &e[0]-ad0 ];
  int p = head;
  while(p!=end){
    if(e==elt[p]){exist=true; break;}
    p = next[p];}
  if(exist){
    loc[p].push(m,nb_elt(m));
    return p;}
  
  p = size(elt);
  elt.push_back(e);
  next.push_back(head);
  head = p;
  
  loc.push_back(loc_<dim>(m,nb_elt(m)));
  return p;
  
}


typedef list_elt_<dim0>  list_elt_0D;
typedef list_elt_<dim1>  list_elt_1D;
typedef list_elt_<dim2>  list_elt_2D;


//====================================//
//  Structures auxiliaires pour       //
//  l'acces aux registres d'elements  //
//====================================//
class geometry;

template <int dim> struct get_elt_{
  inline const vect< elt_<dim> >& apply(const geometry& geom);};

template <int dim> struct get_loc_{
  inline const std::vector<loc_<dim> >& apply(const geometry& geom);};

//==========================//
//  Données géométriques    //
//==========================//


class geometry{
  
 private:
  // variable  de travail
  const char*    meshfile;
  vect<R3>       node;
  list_elt_1D    elt1D;
  list_elt_2D    elt2D;
  
  // Interdit de recopier la geometrie
  geometry(const geometry&);
  
 public:
  geometry ();
  
  
  friend void load_node_gmsh(geometry& , const char*);
  friend int push(geometry&, const mesh_<1>&, elt_1D);
  friend int push(geometry&, const mesh_<2>&, elt_2D);
  
  friend const int nb_node(const geometry&);
  friend const int nb_elt1D(const geometry&);  
  friend const int nb_elt2D(const geometry&); 
  
  friend const char*              meshfile(const geometry&);
  friend const R3&                get_node(const geometry&, const int&);
  friend const vect<R3>&          get_node(const geometry&);
  template<class i_t> friend const subarray<const vect<R3>,i_t> get_node(const geometry&, const i_t&);
  
  friend const vect<elt_1D>&       get_elt1D(const geometry&);
  friend const vect<elt_2D>&       get_elt2D(const geometry&);
  friend struct get_elt_<1>;
  friend struct get_elt_<2>;
  
  friend const std::vector<loc_1D>&     get_loc1D(const geometry&);
  friend const std::vector<loc_2D>&     get_loc2D(const geometry&);
  friend struct get_loc_<1>;
  friend struct get_loc_<2>;
  
};
  
  
// Initialisation et instanciation
inline geometry::geometry(){
  meshfile = 0;
//   vect<R3>      geometry::node;
//   list_elt_1D   geometry::elt1D;
//   list_elt_2D   geometry::elt2D;
// geometry      geom;
}

// Definition des methodes
inline const char*         meshfile(const geometry& geom) {return geom.meshfile;}
inline const int           nb_node (const geometry& geom) {return size(geom.node ); }
inline const int           nb_elt1D(const geometry& geom) {return size(geom.elt1D); }
inline const int           nb_elt2D(const geometry& geom) {return size(geom.elt2D); }


// Acces aux noeuds
inline const vect<R3>&    get_node(const geometry& geom)            {return geom.node;  }
inline const R3&          get_node(const geometry& geom, const int& j){return geom.node[j];}
template <class i_t> const subarray<const vect<R3>,i_t> get_node(const geometry& geom, const i_t& i_) {
  return subarray<const vect<R3>,i_t>(geom.node,i_);}

// Acces aux elements
inline            const vect<elt_1D>& get_elt1D (const geometry& geom)        {return get_elt(geom.elt1D);}
inline            const vect<elt_2D>& get_elt2D (const geometry& geom)        {return get_elt(geom.elt2D);}
template<> inline const vect<elt_1D>& get_elt_<1>::apply(const geometry& geom){return get_elt1D(geom);}
template<> inline const vect<elt_2D>& get_elt_<2>::apply(const geometry& geom){return get_elt2D(geom);}

// Acces aux numerotations locales
inline             const std::vector<loc_1D>& get_loc1D (const geometry& geom)        {return get_loc(geom.elt1D);}
inline             const std::vector<loc_2D>& get_loc2D (const geometry& geom)        {return get_loc(geom.elt2D);}
template <> inline const std::vector<loc_1D>& get_loc_<1>::apply(const geometry& geom){return get_loc1D(geom);}
template <> inline const std::vector<loc_2D>& get_loc_<2>::apply(const geometry& geom){return get_loc2D(geom);}


// Ajout d'elt sans doublonnage
inline int push(geometry& geom, const mesh_1D& m, elt_1D e){return geom.elt1D.push(m,e);}
inline int push(geometry& geom, const mesh_2D& m, elt_2D e){return geom.elt2D.push(m,e);}


// Chargement des noeuds du maillage
inline void load_node_gmsh(geometry& geom,const char* filename){
  
 vect<R3>& node = geom.node;
  assert(!size(node));
  geom.meshfile = filename;
  std::string filename_string = filename;
  filename_string += ".msh";
  
  // lecture du fichier
  std::ifstream file; file.open(filename_string.c_str());
 
  if( file.fail() ){
    std::cout << "Probleme lors du chargement des noeuds" << std::endl;
    exit(EXIT_FAILURE);}
  
  std::string line;
  std::istringstream iss;
  int poubelle;
  R3 p; 
  
  while( line != "$Nodes" ){
    getline(file,line);}      
  
  // Chargt nb de noeuds
  file >> poubelle;
  getline(file,line);       
  
  // Chargt donnees noeuds
  getline(file,line);       
  while( line != "$EndNodes" ){     
    iss.str(line);
    iss >> poubelle;
    iss>>p;
    node.push_back(p);
    iss.clear();
    getline(file,line);
  }
  
  file.close();

  geom.elt1D.init(nb_node(geom), node[0]);  
  geom.elt2D.init(nb_node(geom), node[0]);    
  
}




//============================//
//  Type pour gere si un      //
//  domaine est borne ou non  //
//============================//

enum boundedness{ yes, no };
static const boundedness unbounded = no;

  
//==========================//
//         Maillage         //
//==========================//

template <int dim>
class mesh_{

 public:
  //______________
  // Alias de type
  typedef mesh_<dim>     this_t;
  typedef elt_<dim>      elt_t;
  static const int Dim = dim;
  
 private:
  //_________________________
  // Instances pre-existantes
  geometry&       geom;
//   const vect<R3>&       node;
//   const vect<elt_t>&    elt;
  
  //_______________
  // Donnee membres
  std::vector<int>      num_elt;
  bool                  bounded;    
  
  
  template <class r_t> void operator<<(const r_t& r_){
    int J = push(geom, *this,elt_t(r_));
    num_elt.push_back(J);} 
  
 public:
  
  mesh_<dim>(geometry& geom_): geom(geom_), bounded(true){};
  // Operateur de recopie quand le maillage est vide (pour faire des vecteurs de maillage)
  mesh_<dim>(const mesh_<dim>& m): geom(m.geom), bounded(true){assert((m.num_elt).size()==0);};
  
  template <class m_t> friend void load_elt_gmsh(m_t&, int);
  template <class m_t> friend int  nb_elt (const m_t&);
  const elt_t& operator[](const int& j) const {get_elt_<dim> temp; return (temp.apply(geom))[num_elt[j]];};  
  
  friend const geometry& get_geometry(const mesh_<dim>& m) {return m.geom;};
  
  void operator+=(const this_t& m){
    for(int j=0; j<nb_elt(m); j++){
      int J = push(geom, *this,m[j]);
      num_elt.push_back(J);} }
  
  void operator=(const boundedness& b){
    if(b==no){bounded=false;} }
  
  bool is_bounded() const {return bounded;}
  
  friend void write_medit(const mesh_<dim>& m, char const * const name){
    
    const vect<R3>& node = get_node(m.geom);  
    int nb_node = size(m.node);
    std::ofstream file; file.open(name);    
    file << "MeshVersionFormatted 1\n";
    file << "Dimension 3\n";
    file << "Vertices\n";
    file << nb_node << std::endl;
    for(int j=0; j<nb_node; j++){
      file << m.node[j] << "\t 0 \n";}
    file << std::endl;
    file << "Triangles\n";
    file << nb_elt(m) << std::endl;
    for(int j=0; j<nb_elt(m); j++){  
      file << &m[j][0]-&m.node[0] +1 << "\t";
      file << &m[j][1]-&m.node[0] +1 << "\t";
      file << &m[j][2]-&m.node[0] +1 << "\t";
      file << "0 \n";
    }
    file.close();
    
  }
  
};

template <class m_t>
int nb_elt (const m_t& m){return m.num_elt.size();}

// Routine de chargement de maillage
template <class m_t>
void load_elt_gmsh(m_t& m, int ref = -1){
  
  const int dim = m_t::Dim;
  bemtool::array<dim+1,int> I;  
  const char* filename = meshfile(m.geom);
  std::string filename_string = filename;
  filename_string += ".msh";
  // Variables  auxiliaires
  int poubelle, elt_type;
  int tag, nb_tags;
  
  // Traitement chaines caractere
  std::ifstream file;
  std::string line;
  std::istringstream iss;
  
  file.open(filename_string.c_str(), std::ifstream::in);
  if( file.fail() ){
    std::cout << "Probleme lors du chargement des elements" << std::endl;
    exit(EXIT_FAILURE);}
  
  // Deplacement rubrique elements
  while( line != "$Elements" ){
    getline(file,line);}      
  // Ligne nb elements
  getline(file,line);
  // Chrgt donnees elements
  getline(file,line);       

  int nb_elt = 0;
  while( line != "$EndElements" ){    
    iss.str(line);    
    iss >> poubelle;
    iss >> elt_type;      
    iss >> nb_tags;
    iss >> tag;      
    
    if(elt_type==dim && tag == ref){
      for(int j=0; j<nb_tags-1; j++){
	iss >> poubelle;}      
      
      // acquisition des numeros              
      iss >> I; I--;
      // ajout de l'elt dans le mesh
      m << get_node(m.geom,I);
      
    }
    
    iss.clear();
    getline(file,line);
  }
  
  file.close();
  
}






#endif
