#ifndef MESH_H
#define MESH_H

#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include "calculus.h"
#include "elt.h"


using namespace std;





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
  vector<admesh_t>           admesh;
  vector<int>                num;
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
  vector<loc_t>     loc;
  
  // Premier noeud de la scene
  const R3*         ad0;
  
  // Gestion des doublons
  vector<int>       first;
  vector<int>       next;
  
 public:  
  list_elt_(){};  
  int  push(const mesh_t&, elt_t);
  void init(const int& nb_node, const R3& p0){
    ad0 = &p0;
    assert(first.size()==0);
    first.resize(nb_node,-1);}
  
  friend int                     size(this_t& l)          {return size(l.elt);}
  friend const vect  <elt_t>&    get_elt(const this_t& l) {return l.elt;}
  friend const vector<loc_t>&    get_loc(const this_t& l) {return l.loc;}  
  
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

template <int dim> struct get_elt_{
  static inline const vect< elt_<dim> >& apply();};

template <int dim> struct get_loc_{
  static inline const vector<loc_<dim> >& apply();};


//==========================//
//  Données géométriques    //
//==========================//

class geometry{
  
 private:
  // variable  de travail
  static const char*    meshfile;
  static vect<R3>       node;
  static list_elt_1D    elt1D;
  static list_elt_2D    elt2D;
  
 public:
  friend void load_node(const char*);
  static int push(const mesh_<1>&, elt_1D);
  static int push(const mesh_<2>&, elt_2D);
  
  friend int nb_node();
  friend int nb_elt1D();  
  friend int nb_elt2D(); 
  
  friend const char*              meshfile();
  friend const R3&                get_node(const int&);
  friend const vect<R3>&          get_node();
  template<class i_t> friend const subarray<const vect<R3>,i_t> get_node(const i_t&);
  
  friend const vect<elt_1D>&       get_elt1D();
  friend const vect<elt_2D>&       get_elt2D();
  friend struct get_elt_<1>;
  friend struct get_elt_<2>;
  
  friend const vector<loc_1D>&     get_loc1D();
  friend const vector<loc_2D>&     get_loc2D();
  friend struct get_loc_<1>;
  friend struct get_loc_<2>;
  
};


// Initialisation et instanciation
const char*   geometry::meshfile = 0;
vect<R3>      geometry::node;
list_elt_1D   geometry::elt1D;
list_elt_2D   geometry::elt2D;
geometry      geom;


// Definition des methodes
const char*   meshfile() {return geometry::meshfile;    }
int           nb_node () {return size(geometry::node ); }
int           nb_elt1D() {return size(geometry::elt1D); }
int           nb_elt2D() {return size(geometry::elt2D); }


// Acces aux noeuds
const vect<R3>&    get_node()            {return geometry::node;  }
const R3&          get_node(const int& j){return geometry::node[j];}
template <class i_t> const subarray<const vect<R3>,i_t> get_node(const i_t& i_) {
  return subarray<const vect<R3>,i_t>(geometry::node,i_);}

// Acces aux elements
           const vect<elt_1D>& get_elt1D ()        {return get_elt(geometry::elt1D);}
           const vect<elt_2D>& get_elt2D ()        {return get_elt(geometry::elt2D);}
template<> const vect<elt_1D>& get_elt_<1>::apply(){return get_elt(geometry::elt1D);}
template<> const vect<elt_2D>& get_elt_<2>::apply(){return get_elt(geometry::elt2D);}

// Acces aux numerotations locales
            const vector<loc_1D>& get_loc1D ()        {return get_loc(geometry::elt1D);}
            const vector<loc_2D>& get_loc2D ()        {return get_loc(geometry::elt2D);}
template <> const vector<loc_1D>& get_loc_<1>::apply(){return get_loc(geometry::elt1D);}
template <> const vector<loc_2D>& get_loc_<2>::apply(){return get_loc(geometry::elt2D);}


// Ajout d'elt sans doublonnage
int geometry::push(const mesh_1D& m, elt_1D e){return elt1D.push(m,e);}
int geometry::push(const mesh_2D& m, elt_2D e){return elt2D.push(m,e);}


// Chargement des noeuds du maillage
void load_node(const char* filename){
  
   vect<R3>& node = geometry::node;
  assert(!size(node));
  geometry::meshfile = filename;
  
  // lecture du fichier
  ifstream file; file.open(filename);
  if( file.fail() ){
    cout << "fichier de maillage non trouve" << endl;
    exit(EXIT_FAILURE);}
  
  string line;
  istringstream iss;
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

  geometry::elt1D.init(nb_node(), node[0]);  
  geometry::elt2D.init(nb_node(), node[0]);    
  
}


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
  const vect<R3>&       node;
  const vect<elt_t>&    elt;
  
  //_______________
  // Donnee membres
  vector<int>           num_elt;
  
  
  // Pas de constructeur par recopie
  mesh_<dim>(const mesh_<dim>&);
  
  template <class r_t> void operator<<(const r_t& r_){
    int J = geometry::push(*this,elt_t(r_));
    num_elt.push_back(J);} 
  
 public:
  
  mesh_<dim>(): node(get_node()), elt(get_elt_<dim>::apply()){};
  template <class m_t> friend void load(m_t&, int);
  template <class m_t> friend int  nb_elt (const m_t&);
  const elt_t& operator[](const int& j) const { return elt[num_elt[j]];};  
  
  void operator+=(const this_t& m){
    for(int j=0; j<nb_elt(m); j++){
      int J = geometry::push(*this,m[j]);
      num_elt.push_back(J);} }
  
  friend void write(const mesh_<dim>& m, char const * const name){
        
    const vect<R3>& node = get_node();  
    int nb_node = size(node);
    ofstream file; file.open(name);    
    file << "MeshVersionFormatted 1\n";
    file << "Dimension 3\n";
    file << "Vertices\n";
    file << nb_node << endl;
    for(int j=0; j<nb_node; j++){
      file << node[j] << "\t 0 \n";}
    file << endl;
    file << "Triangles\n";
    file << nb_elt(m) << endl;
    for(int j=0; j<nb_elt(m); j++){  
      file << &m[j][0]-&node[0] +1 << "\t";
      file << &m[j][1]-&node[0] +1 << "\t";
      file << &m[j][2]-&node[0] +1 << "\t";
      file << "0 \n";
    }
    file.close();
    
  }
  
  
  
};

template <class m_t>
int nb_elt (const m_t& m){return m.num_elt.size();}

// Routine de chargement de maillage
template <class m_t>
void load(m_t& m, int ref = -1){
  
  const int dim = m_t::Dim;
  array<dim+1,int> I;  
  const char* filename = meshfile();
  
  // Variables  auxiliaires
  int poubelle, elt_type;
  int tag, nb_tags;
  
  // Traitement chaines caractere
  ifstream file;
  string line;
  istringstream iss;
  
  file.open(filename, ifstream::in);
  if( file.fail() ){
    cout << "fichier de maillage non trouve" << endl;
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
      m << get_node(I);
      
    }
    
    iss.clear();
    getline(file,line);
  }
  
  file.close();
  
}






#endif
