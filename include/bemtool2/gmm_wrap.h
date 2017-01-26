#ifndef GMM_WRAP_H
#define GMM_WRAP_H

#include <gmm/gmm.h>
#include "calculus.h"

namespace bemtool{



/*===========================
||   Interface avec gmm    ||
===========================*/

template <class gmm_mat_t,class val_type>
  class gmm_mat{
  
 private:
  gmm_mat_t    mat;
  std::vector<std::size_t> ipvt;
  
 public:
  typedef gmm_mat<gmm_mat_t,val_type>   mat_t;
  typedef vect<val_type>               vect_t;
  typedef val_type                        v_t;
  typedef typename gmm::linalg_traits<gmm_mat_t>::reference  ref_t;
  
  gmm_mat(){};
 gmm_mat(const int n, const int m): mat(n,m), ipvt(n){};
  
  ref_t operator()(const int i, const int j){    
    return mat(i,j);}
  
  const ref_t operator()(const int i, const int j) const{
    return mat(i,j);}
  
  template <class row_t, class col_t>
    submat<mat_t,row_t,col_t> 
    operator()(const row_t& I0,const col_t& J0 ){
    return submat<mat_t,row_t,col_t>(*this,I0,J0);}
  
  template <class row_t, class col_t>
    submat<const mat_t,row_t,col_t> 
    operator()(const row_t& I0,const col_t& J0 ) const{
    return submat<const mat_t,row_t,col_t>(*this,I0,J0);}
  
  friend int nb_rows(const mat_t& m){
    return gmm::mat_nrows(m.mat);}
  
  friend int nb_cols(const mat_t& m){
    return gmm::mat_ncols(m.mat);}
  
  friend void clear(mat_t& m){
    gmm::clear(m.mat);}
  
  ///////////////////////////////////
  //    produit matrice-vecteur    //
  ///////////////////////////////////
  
  friend void mv_prod(vect_t& b, const mat_t& A, const vect_t& x){
    assert(size(b)==nb_rows(A));        
    gmm::mult(A.mat,vector_of(x),vector_of(b));}
  
  friend void add_mv_prod(vect_t& b, const mat_t& A, const vect_t& x){
    assert(size(b)==nbrows(A));        
    gmm::mult_add(A.mat,vector_of(x),vector_of(b));}
  
  ///////////////////////////////////
  //            solveurs           //
  ///////////////////////////////////
  
  friend void lu_factor(mat_t& m, mat_t& lu){
    assert(nbrows(lu)==nbcols(lu));
    lu.mat = m.mat;
    gmm::lu_factor(lu.mat,lu.ipvt);}
  
  friend void lu_solve(mat_t& lu, vect<v_t>& x, vect<v_t>& b){
    assert(x.size()==b.size());
    gmm::lu_solve(lu.mat,lu.ipvt,vector_of(x),vector_of(b));} 
  
  friend v_t cond(mat_t& m){
    gmm::condition_number(m.mat);}

  friend void invert(mat_t& m){
    gmm::lu_inverse(m.mat);}
  
  friend void mult(mat_t& m1, mat_t& m2,mat_t& m3){
    gmm::mult(m1.mat,m2.mat,m3.mat);}
  
  friend void cg_solve(mat_t& m, vect<v_t>& x, vect<v_t>& b){
    assert(size(x)==size(b));    
    gmm::identity_matrix Id;
    gmm::iteration iter(10E-9);
    iter.set_noisy(1);
    gmm::cg(m.mat,vector_of(x),vector_of(b),Id,Id,iter);
  }
  
  friend void gmres_solve(mat_t& m, vect<v_t>& x, vect<v_t>& b, const int restart, int verbose=1){
    assert(size(x)==size(b));    
    gmm::identity_matrix PR;
    //    gmm::diagonal_precond<gmm_mat_t> PR(m.mat); 
    //    gmm::ildlt_precond<gmm_mat_t> PR(m.mat);     
    gmm::iteration iter(10E-9);
    iter.set_name("OUTER");
    iter.set_noisy(verbose);
    std::size_t restart2 = restart;
    gmm::gmres(m.mat,vector_of(x),vector_of(b),PR,restart,iter);
  }
  
  
  /////////////////////////////////////
  //     Ecriture dans un fichier    //
  /////////////////////////////////////  
  
  friend void write(mat_t& m, char const * const name){

    std::ofstream file; file.open(name);
    int nb_row = nb_rows(m);
    int nb_col = nb_cols(m);
    
    file << "# name: " << name       << std::endl;
    file << "# type: complex matrix" << std::endl;
    file << "# rows:   \t" << nb_row << std::endl;
    file << "# columns:\t" << nb_col << std::endl;  
    
    int count = 0;
    for(int k=0; k<nb_col; k++){
      for(int j=0; j<nb_row; j++){
	file << m(j,k) << "  ";
	count++;
	if(count==10){
	  file << std::endl; 
	  count = 0;}
      }
    }
    file.close();
  }


  
  
};


/*===========================
||   Definitions de type   ||
===========================*/

typedef Cplx Field;
typedef gmm_mat< gmm::row_matrix< gmm::wsvector<Field> >,Field > gmm_sparse;
typedef gmm_mat< gmm::dense_matrix<Field>,Field>                 gmm_dense;

}

#endif
