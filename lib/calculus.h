#ifndef CALCULUS_H
#define CALCULUS_H

#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <cstdarg>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>





using namespace std;

template <int N> struct dim_{ };
static const int dim0 = 0;
static const int dim1 = 1;
static const int dim2 = 2;
static const int dim3 = 3;

//==========================//
//      Types de base       //
//==========================//

typedef double        Real;
typedef complex<Real> Cplx;
const   Cplx iu(0.,1.);
const   Real pi = 3.14159265358979;

inline Cplx operator+ (const int& n,  const Cplx& z){ return (double)n + z;}
inline Cplx operator+ (const Cplx& z, const int& n) { return z + (double)n;}

inline Cplx operator- (const int& n,  const Cplx& z){ return (double)n - z;}
inline Cplx operator- (const Cplx& z, const int& n) { return z - (double)n;}

inline Cplx operator* (const int& n,  const Cplx& z){ return (double)n * z;}
inline Cplx operator* (const Cplx& z, const int& n) { return z * (double)n;}

inline Real conj(const Real& x){return x;}
inline int  conj(const int& x) {return x;}

template <class U> struct isbase               {static const bool test = false; };
template <>        struct isbase<bool>         {static const bool test = true;  };
template <>        struct isbase<int>          {static const bool test = true;  };
template <>        struct isbase<double>       {static const bool test = true;  };
template <>        struct isbase<float>        {static const bool test = true;  };
template <class U> struct isbase<complex<U> >  {static const bool test = true;  };

// Operateurs d'acces
template <class r_t, bool test = isbase<r_t>::test > struct access_ {
  static       typename r_t::v_t& op(      r_t& r_, const int& j)              {return r_[j];}
  static const typename r_t::v_t& op(const r_t& r_, const int& j)              {return r_[j];}
  static       typename r_t::v_t& op(      r_t& r_, const int& j, const int& k){return r_(j,k);}
  static const typename r_t::v_t& op(const r_t& r_, const int& j, const int& k){return r_(j,k);}
};

template <class r_t> struct access_<r_t,true> {
  static       r_t& op(      r_t& r_, const int& j)              {return r_;}
  static const r_t& op(const r_t& r_, const int& j)              {return r_;}
  static       r_t& op(      r_t& r_, const int& j, const int& k){return r_;}
  static const r_t& op(const r_t& r_, const int& j, const int& k){return r_;}
};

// Type des entrees
template <class r_t, bool test = isbase<r_t>::test> struct entry           {typedef typename r_t::v_t type; };
template <class r_t>                                struct entry<r_t,true> {typedef r_t               type; };

// Gestion des operations entre type
template <class l_t, class r_t> struct res                     {typedef l_t          type;};
template <>                     struct res< int, double >      {typedef double       type;};
template <>                     struct res< double, int >      {typedef double       type;};
template <class l_t>            struct res< l_t,complex<l_t> > {typedef complex<l_t> type;};
template <class r_t>            struct res< complex<r_t>,r_t > {typedef complex<r_t> type;};
template <class l_t>            struct res< complex<l_t>,int > {typedef complex<l_t> type;};
template <class r_t>            struct res< int,complex<r_t> > {typedef complex<r_t> type;};


template <class l_t, class r_t> struct resop{
  typedef typename res<typename entry<l_t>::type,typename entry<r_t>::type>::type  type; };


template <class l_t, class r_t, bool l_test = isbase<l_t>::test, bool r_test = isbase<r_t>::test >
  struct find_nb_rows{ static const int nr = l_t::nr;};
template <class l_t, class r_t >
  struct find_nb_rows<l_t, r_t, true, false>{ static const int nr = r_t::nr;};



//==========================//
//        Boucles           //
//==========================//

//---------------------------//
template <class r_t, int D = r_t::dim> struct ascending_loop{
  static inline void apply(r_t& r_){
    r_[D-1]=D-1; ascending_loop<r_t,D-1>::apply(r_);} };

template <class r_t> struct ascending_loop<r_t,1>{
  static inline void apply(r_t& r_){r_[0]=0;} };

//---------------------------//
template <class r_t, int D = r_t::dim> struct decrement_loop{
  static inline void apply(r_t& r_){
    r_[D-1]--; decrement_loop<r_t,D-1>::apply(r_);} };

template <class r_t> struct decrement_loop<r_t,1>{
  static inline void apply(r_t& r_){r_[0]--;} };

//---------------------------//
template <class r_t, int D = r_t::dim> struct increment_loop{
  static inline void apply(r_t& r_){
    r_[D-1]++; increment_loop<r_t,D-1>::apply(r_);} };

template <class r_t> struct increment_loop<r_t,1>{
  static inline void apply(r_t& r_){r_[0]++;} };

//---------------------------//
template <class r_t, int D = r_t::dim> struct istream_loop{
  static inline void apply(istream& is, r_t& r_){
    is >> r_[r_t::dim-D]; istream_loop<r_t,D-1>::apply(is,r_);} };

template <class r_t> struct istream_loop<r_t,1>{
  static inline void apply(istream& is, r_t& r_){
    is >> r_[r_t::dim-1];} };

//---------------------------//
template <class r_t, int d = r_t::dim> struct ostream_loop{
  static inline void apply(ostream& os, const r_t& r_){
    os<<r_[r_t::dim-d] << "\t"; ostream_loop<r_t,d-1>::apply(os,r_);} };

template <class r_t> struct ostream_loop<r_t,1>{
  static inline void apply(ostream& os, const r_t& r_){
    os<<r_[r_t::dim-1];} };

//---------------------------//
template <class l_t, class r_t, int D = l_t::dim> struct assign_loop{
  static void apply(l_t& l_, const r_t& r_){
    l_[D-1] = access_<r_t>::op(r_,D-1); assign_loop<l_t,r_t,D-1>::apply(l_,r_);} };

template <class l_t, class r_t> struct assign_loop<l_t,r_t,1>{
  static void apply(l_t& l_, const r_t& r_){ l_[0] = access_<r_t>::op(r_,0);} };

//---------------------------//
template <class l_t, class r_t, int D = l_t::dim> struct plus_assign_loop{
  static void apply(l_t& l_, const r_t& r_){
    l_[D-1] += access_<r_t>::op(r_,D-1); plus_assign_loop<l_t,r_t,D-1>::apply(l_,r_);} };

template <class l_t, class r_t> struct plus_assign_loop<l_t,r_t,1>{
  static void apply(l_t& l_, const r_t& r_){ l_[0] += access_<r_t>::op(r_,0);} };

//---------------------------//
template <class l_t, int d = l_t::dim> struct construct_loop{
  static void apply(l_t& l_){
    l_[d-1] = typename l_t::v_t();
    construct_loop<l_t,d-1>::apply(l_);} };

template <class l_t> struct construct_loop<l_t,1>{
  static void apply(l_t& l_){
    l_[0] = typename l_t::v_t();} };




//---------------------------//
template <class l_t, int j = l_t::nr, int k = l_t::nc> struct construct_dbloop{
  static inline void apply(l_t& l_){
    l_(j-1,k-1) = typename l_t::v_t(); construct_dbloop<l_t,j,k-1>::apply(l_); }};

template <class l_t, int j> struct construct_dbloop<l_t,j,1>{
  static inline void apply(l_t& l_){
    l_(j-1,0) = typename l_t::v_t(); construct_dbloop<l_t,j-1,l_t::nc>::apply(l_); }};

template <class l_t> struct construct_dbloop<l_t,1,1>{
  static inline void apply(l_t& l_){
    l_(0,0) = typename l_t::v_t(); }};

//---------------------------//
template <class l_t, class r_t, int j = l_t::nr, int k = l_t::nc> struct assign_dbloop{
  static inline void apply(l_t& l_, const r_t& r_){
    l_(j-1,k-1) = access_<r_t>::op(r_,j-1,k-1); assign_dbloop<l_t,r_t,j,k-1>::apply(l_,r_); }};

template <class l_t, class r_t, int j> struct assign_dbloop<l_t,r_t,j,1>{
  static inline void apply(l_t& l_, const r_t& r_){
    l_(j-1,0) = access_<r_t>::op(r_,j-1,0); assign_dbloop<l_t,r_t,j-1,l_t::nc>::apply(l_,r_); }};

template <class l_t, class r_t> struct assign_dbloop<l_t,r_t,1,1>{
  static inline void apply(l_t& l_, const r_t& r_){
    l_(0,0) = access_<r_t>::op(r_,0,0); }};

//---------------------------//
template <class l_t, class r_t, int j = l_t::nr, int k = l_t::nc> struct plus_assign_dbloop{
  static inline void apply(l_t& l_, const r_t& r_){
    l_(j-1,k-1) += access_<r_t>::op(r_,j-1,k-1); plus_assign_dbloop<l_t,r_t,j,k-1>::apply(l_,r_); }};

template <class l_t, class r_t, int j> struct plus_assign_dbloop<l_t,r_t,j,1>{
  static inline void apply(l_t& l_, const r_t& r_){
    l_(j-1,0) += access_<r_t>::op(r_,j-1,0); plus_assign_dbloop<l_t,r_t,j-1,l_t::nc>::apply(l_,r_); }};

template <class l_t, class r_t> struct plus_assign_dbloop<l_t,r_t,1,1>{
  static inline void apply(l_t& l_, const r_t& r_){
    l_(0,0) += access_<r_t>::op(r_,0,0); }};






//---------------------------//
template <class r_t, int nr = r_t::nr, int nc = r_t::nc> struct ostream_dbloop{
  static inline void apply(ostream& os, const r_t& r_){
    os<<r_(r_t::nr-nr,r_t::nc-nc) << "\t"; ostream_dbloop<r_t,nr,nc-1>::apply(os,r_);} };

template <class r_t, int nr> struct ostream_dbloop<r_t,nr,1>{
  static inline void apply(ostream& os, const r_t& r_){
    os<<r_(r_t::nr-nr,r_t::nc-1) << endl; ostream_dbloop<r_t,nr-1,r_t::nc>::apply(os,r_);} };

template <class r_t> struct ostream_dbloop<r_t,1,1>{
  static inline void apply(ostream& os, const r_t& r_){
    os<<r_(r_t::nr-1,r_t::nc-1);} };

//---------------------------//
template <class v_t, class l_t, class r_t, int d = l_t::nc> struct mat_mult{

  static v_t apply(const l_t& l_, const r_t& r_, const int& j){
    return access_<l_t>::op(l_,j,d-1)*access_<r_t>::op(r_,d-1) + mat_mult<v_t,l_t,r_t,d-1>::apply(l_,r_,j);}

  static v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
    return access_<l_t>::op(l_,j,d-1)*access_<r_t>::op(r_,d-1,k) + mat_mult<v_t,l_t,r_t,d-1>::apply(l_,r_,j);}
  
};

template <class v_t, class l_t, class r_t> struct mat_mult<v_t,l_t,r_t,1>{

  static v_t apply(const l_t& l_, const r_t& r_, const int& j){
    return access_<l_t>::op(l_,j,0)*access_<r_t>::op(r_,0);}

  static v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
    return access_<l_t>::op(l_,j,0)*access_<r_t>::op(r_,0,k);}
  
};


//---------------------------//
template <class l_t,class r_t, int k> struct dprod{
  static typename resop<l_t,r_t>::type apply(const l_t& l_, const r_t& r_){
    return l_[k-1]*conj(r_[k-1]) +  dprod<l_t,r_t,k-1>::apply(l_,r_);}};

template <class l_t, class r_t> struct dprod<l_t,r_t,1>{
  static typename resop<l_t,r_t>::type apply(const l_t& l_, const r_t& r_){
    return l_[0]*conj(r_[0]);}};




//==========================//
//       Operations         //
//==========================//


//%%%%%%%%%%//
// Addition //
//%%%%%%%%%%//

template <class lhs_t, class rhs_t>
struct pp{
  
  typedef lhs_t                          l_t;
  typedef rhs_t                          r_t;  
  typedef typename resop<l_t,r_t>::type  v_t;


  static v_t apply(const l_t& l_, const r_t& r_, const int& j){
    return access_<l_t>::op(l_,j) + access_<r_t>::op(r_,j);}

  static v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
    return access_<l_t>::op(l_,j,k) + access_<r_t>::op(r_,j,k);}
  
};


//%%%%%%%%%%%%%%//
// Soustraction //
//%%%%%%%%%%%%%%//

template <class lhs_t, class rhs_t>
struct mm{

  typedef lhs_t                          l_t;
  typedef rhs_t                          r_t;  
  typedef typename resop<l_t,r_t>::type  v_t;
  
  static v_t apply(const l_t& l_, const r_t& r_, const int& j){
    return access_<l_t>::op(l_,j) - access_<r_t>::op(r_,j);}

  static v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
    return access_<l_t>::op(l_,j,k) - access_<r_t>::op(r_,j,k);}
};


//%%%%%%%%%%%%%%%%//
// Multiplication //
//%%%%%%%%%%%%%%%%//

template <class lhs_t, class rhs_t, bool l_test = isbase<lhs_t>::test, bool r_test = isbase<rhs_t>::test>
struct tt{
  
  typedef lhs_t                            l_t;
  typedef rhs_t                            r_t;  
  typedef typename resop<l_t,r_t>::type    v_t;
  
  static v_t apply(const l_t& l_, const r_t& r_, const int& j){
    return mat_mult<v_t,l_t,r_t>::apply(l_,r_,j);}
  
  static v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
    return mat_mult<v_t,l_t,r_t>::apply(l_,r_,j,k);}
  
};


template <class lhs_t, class rhs_t>
struct tt<lhs_t,rhs_t, true, false>{

  typedef lhs_t                            l_t;
  typedef rhs_t                            r_t;  
  typedef typename resop<l_t,r_t>::type    v_t;
  
  static v_t apply(const l_t& l_, const r_t& r_, const int& j){
    return access_<l_t>::op(l_,j)*access_<r_t>::op(r_,j);}

  static v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
    return access_<l_t>::op(l_,j,k)*access_<r_t>::op(r_,j,k);}
  
};

template <class lhs_t, class rhs_t>
struct tt<lhs_t,rhs_t, false, true>{

  typedef lhs_t                            l_t;
  typedef rhs_t                            r_t;  
  typedef typename resop<l_t,r_t>::type    v_t;
  
  static v_t apply(const l_t& l_, const r_t& r_, const int& j){
    return access_<l_t>::op(l_,j)*access_<r_t>::op(r_,j);}

  static v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
    return access_<l_t>::op(l_,j,k)*access_<r_t>::op(r_,j,k);}
  
};



//%%%%%%%%%%%%%%%%//
//  Comparaisons  //
//%%%%%%%%%%%%%%%%//

template <class l_t, class r_t> bool operator==(const l_t& l_, const r_t& r_){
  for(int j=0; j<l_t::dim; j++){ if( access_<l_t>::op(l_,j)!=access_<r_t>::op(r_,j) ){return false;}}
  return true;}

template <class l_t, class r_t> bool operator!=(const l_t& l_, const r_t& r_){
  return !(l_==r_);}



//==========================//
//       Expressions        //
//==========================//

template <class op>
class xpr{
  
 public: 
  
  typedef typename op::l_t l_t;
  typedef typename op::r_t r_t; 
  typedef typename op::v_t v_t; 
  typedef xpr<op>       this_t;
  
  static const int nr = find_nb_rows<l_t,r_t>::nr;
  static const int dim = nr;
  
 private:
  
  const l_t& l_;
  const r_t& r_;  
  
 public:
  
 xpr(const l_t& l0, const r_t& r0): l_(l0), r_(r0) {};
  
  v_t operator[](const int& j){ return op::apply(l_,r_,j); }
  const v_t operator[](const int& j) const { return op::apply(l_,r_,j); }  
  
  v_t operator()(const int& j, const int& k){ return op::apply(l_,r_,j,k); }
  const v_t operator()(const int& j, const int& k) const { return op::apply(l_,r_,j,k); }
  
  //==== Addition
  template <class r_t> xpr< pp<this_t,r_t> > operator+(const r_t& r_) const {return xpr< pp<this_t,r_t> >(*this,r_);}
  xpr< pp<this_t, int>  > operator+(const int&  r_) const {return xpr< pp<this_t, int>  >(*this,r_);}
  xpr< pp<this_t, Real> > operator+(const Real& r_) const {return xpr< pp<this_t, Real> >(*this,r_);}
  xpr< pp<this_t, Cplx> > operator+(const Cplx& r_) const {return xpr< pp<this_t, Cplx> >(*this,r_);}   

  inline friend xpr< pp<int,this_t>  > operator+(const int&  l_, const this_t& r_){return xpr< pp<int,this_t>  >(l_,r_);}
  inline friend xpr< pp<Real,this_t> > operator+(const Real& l_, const this_t& r_){return xpr< pp<Real,this_t> >(l_,r_);}
  inline friend xpr< pp<Cplx,this_t> > operator+(const Cplx& l_, const this_t& r_){return xpr< pp<Cplx,this_t> >(l_,r_);}

  //==== Soustraction
  template <class r_t> xpr< mm<this_t,r_t> > operator-(const r_t& r_) const {return xpr< mm<this_t,r_t> >(*this,r_);}
  xpr< mm<this_t, int>  > operator-(const int&  r_) const {return xpr< mm<this_t, int>  >(*this,r_);}
  xpr< mm<this_t, Real> > operator-(const Real& r_) const {return xpr< mm<this_t, Real> >(*this,r_);}
  xpr< mm<this_t, Cplx> > operator-(const Cplx& r_) const {return xpr< mm<this_t, Cplx> >(*this,r_);}   

  inline friend xpr< mm<int,this_t>  > operator-(const int&  l_, const this_t& r_){return xpr< mm<int,this_t>  >(l_,r_);}
  inline friend xpr< mm<Real,this_t> > operator-(const Real& l_, const this_t& r_){return xpr< mm<Real,this_t> >(l_,r_);}
  inline friend xpr< mm<Cplx,this_t> > operator-(const Cplx& l_, const this_t& r_){return xpr< mm<Cplx,this_t> >(l_,r_);}
  
  //==== Multiplication
  template <class r_t> xpr< tt<this_t,r_t> > operator*(const r_t& r_) const {return xpr< tt<this_t,r_t> >(*this,r_);}
  xpr< tt<this_t, int>  > operator*(const int&  r_) const {return xpr< tt<this_t, int>  >(*this,r_);}
  xpr< tt<this_t, Real> > operator*(const Real& r_) const {return xpr< tt<this_t, Real> >(*this,r_);}
  xpr< tt<this_t, Cplx> > operator*(const Cplx& r_) const {return xpr< tt<this_t, Cplx> >(*this,r_);}   
  
  inline friend xpr< tt<this_t,int>  > operator*(const int&  l_, const this_t& r_){return xpr< tt<this_t,int>  >(r_,l_);}
  inline friend xpr< tt<this_t,Real> > operator*(const Real& l_, const this_t& r_){return xpr< tt<this_t,Real> >(r_,l_);}
  inline friend xpr< tt<this_t,Cplx> > operator*(const Cplx& l_, const this_t& r_){return xpr< tt<this_t,Cplx> >(r_,l_);}
  
  //===== Produit scalaire
  template <class r_t>typename resop<this_t,r_t>::type operator,(const r_t& r_) const {
    assert(nr==r_t::nr); return dprod<this_t,r_t,nr>::apply(*this,r_); }
  
};

template <class t> struct access_<xpr<t>, false> {
  static inline typename xpr<t>::v_t op(xpr<t>& r_, const int& j){return r_[j];}
  static inline const typename xpr<t>::v_t op(const xpr<t>& r_, const int& j){return r_[j];}
  
  static inline typename xpr<t>::v_t op(xpr<t>& r_, const int& j, const int& k){return r_(j,k);}
  static inline const typename xpr<t>::v_t op(const xpr<t>& r_, const int& j, const int& k){return r_(j,k);}
};



//==========================//
//        Sub_Array         //
//==========================//

template <class a_t, class i_t>
class subarray{
  
public:
  static const int dim  =  i_t::dim;
  static const int nr   =  i_t::dim;
  typedef typename a_t::v_t     v_t;
  typedef subarray<a_t,i_t>  this_t;

private:
  a_t&         a_;
  const i_t&   i_;
  
public:
  
  subarray(a_t& a0, const i_t& i0): a_(a0), i_(i0) {};
  
  template <class r_t> void operator=(const r_t& r_){
    assign_loop<this_t,r_t>::apply(*this,r_);} 

  template <class r_t> void operator+=(const r_t& r_){
    plus_assign_loop<this_t,r_t>::apply(*this,r_);} 
  
  v_t& operator[](const int& j){return a_[i_[j]];}
  const v_t& operator[](const int& j) const {return a_[i_[j]];}
  
  friend ostream& operator<<(ostream& os, const this_t& ar){
    ostream_loop<this_t>::apply(os,ar); return os;}
  
  friend this_t& operator>>(istream& is, this_t& ar){
    istream_loop<this_t>::apply(is,ar); return ar;}
  
  void operator++(int){increment_loop<this_t>::apply(*this);}
  
  void operator--(int){decrement_loop<this_t>::apply(*this);}  
  
  inline friend int size(const this_t& ar){return dim;} 

  //==== Addition
  template <class r_t> xpr< pp<this_t,r_t> > operator+(const r_t& r_) const {return xpr< pp<this_t,r_t> >(*this,r_);}
  xpr< pp<this_t, int>  > operator+(const int&  r_) const {return xpr< pp<this_t, int>  >(*this,r_);}
  xpr< pp<this_t, Real> > operator+(const Real& r_) const {return xpr< pp<this_t, Real> >(*this,r_);}
  xpr< pp<this_t, Cplx> > operator+(const Cplx& r_) const {return xpr< pp<this_t, Cplx> >(*this,r_);}   

  inline friend xpr< pp<int,this_t>  > operator+(const int&  l_, const this_t& r_){return xpr< pp<int,this_t>  >(l_,r_);}
  inline friend xpr< pp<Real,this_t> > operator+(const Real& l_, const this_t& r_){return xpr< pp<Real,this_t> >(l_,r_);}
  inline friend xpr< pp<Cplx,this_t> > operator+(const Cplx& l_, const this_t& r_){return xpr< pp<Cplx,this_t> >(l_,r_);}
  
  //==== Soustraction
  template <class r_t> xpr< mm<this_t,r_t> > operator-(const r_t& r_) const {return xpr< mm<this_t,r_t> >(*this,r_);}
  xpr< mm<this_t, int>  > operator-(const int&  r_) const {return xpr< mm<this_t, int>  >(*this,r_);}
  xpr< mm<this_t, Real> > operator-(const Real& r_) const {return xpr< mm<this_t, Real> >(*this,r_);}
  xpr< mm<this_t, Cplx> > operator-(const Cplx& r_) const {return xpr< mm<this_t, Cplx> >(*this,r_);}   

  inline friend xpr< mm<int,this_t>  > operator-(const int&  l_, const this_t& r_){return xpr< mm<int,this_t>  >(l_,r_);}
  inline friend xpr< mm<Real,this_t> > operator-(const Real& l_, const this_t& r_){return xpr< mm<Real,this_t> >(l_,r_);}
  inline friend xpr< mm<Cplx,this_t> > operator-(const Cplx& l_, const this_t& r_){return xpr< mm<Cplx,this_t> >(l_,r_);}
  
  //==== Multiplication
  template <class r_t> xpr< tt<this_t,r_t> > operator*(const r_t& r_) const {return xpr< tt<this_t,r_t> >(*this,r_);}
  xpr< tt<this_t, int>  > operator*(const int&  r_) const {return xpr< tt<this_t, int>  >(*this,r_);}
  xpr< tt<this_t, Real> > operator*(const Real& r_) const {return xpr< tt<this_t, Real> >(*this,r_);}
  xpr< tt<this_t, Cplx> > operator*(const Cplx& r_) const {return xpr< tt<this_t, Cplx> >(*this,r_);}   
  
  inline friend xpr< tt<int,this_t>  > operator*(const int&  l_, const this_t& r_){return xpr< tt<int,this_t>  >(l_,r_);}
  inline friend xpr< tt<Real,this_t> > operator*(const Real& l_, const this_t& r_){return xpr< tt<Real,this_t> >(l_,r_);}
  inline friend xpr< tt<Cplx,this_t> > operator*(const Cplx& l_, const this_t& r_){return xpr< tt<Cplx,this_t> >(l_,r_);}

  //===== Produit scalaire
  template <class r_t>typename resop<this_t,r_t>::type operator,(const r_t& r_) const {
    assert(nr==r_t::nr); return dprod<this_t,r_t,nr>::apply(*this,r_); }

  
};


//==========================//
//           Array          //
//==========================//

template <int Dim, class V_t>
class array{
  
public:
  static const int dim =      Dim;
  static const int nr  =      Dim;
  typedef array<Dim,V_t>   this_t;
  typedef V_t                 v_t;
  
private:
  v_t  v_[dim];
  
public:
  array<Dim,V_t>(){ construct_loop<this_t>::apply(*this); }
  
  template <class r_t> array<Dim,V_t>(const r_t& r_){
    assign_loop<this_t,r_t>::apply(*this,r_);}   
  
  template <class r_t> void operator=(const r_t& r_){
    assign_loop<this_t,r_t>::apply(*this,r_);} 

  template <class r_t> void operator+=(const r_t& r_){
    plus_assign_loop<this_t,r_t>::apply(*this,r_);} 
  
  v_t& operator[](const int& j){return v_[j];}
  const v_t& operator[](const int& j) const {return v_[j];}
  
  inline friend ostream& operator<<(ostream& os, const this_t& ar){
    ostream_loop<this_t>::apply(os,ar); return os;}  
  
  inline friend this_t& operator>>(istream& is, this_t& ar){
    istream_loop<this_t>::apply(is,ar); return ar;}
  
  void operator++(int){increment_loop<this_t>::apply(*this);}
  
  void operator--(int){decrement_loop<this_t>::apply(*this);}  
  
  inline friend int size(const this_t& ar){return dim;}   

  template <class i_t> subarray<this_t,i_t> operator[] (const i_t& i_){
    return subarray<this_t,i_t>(*this,i_);}

  template <class i_t> subarray<const this_t,i_t> operator[] (const i_t& i_) const {
    return subarray<const this_t,i_t>(*this,i_);}

  //==== Addition
  template <class r_t> xpr< pp<this_t,r_t> > operator+(const r_t& r_) const {return xpr< pp<this_t,r_t> >(*this,r_);}
  xpr< pp<this_t, int>  > operator+(const int&  r_) const {return xpr< pp<this_t, int>  >(*this,r_);}
  xpr< pp<this_t, Real> > operator+(const Real& r_) const {return xpr< pp<this_t, Real> >(*this,r_);}
  xpr< pp<this_t, Cplx> > operator+(const Cplx& r_) const {return xpr< pp<this_t, Cplx> >(*this,r_);}   

  inline friend xpr< pp<int,this_t>  > operator+(const int&  l_, const this_t& r_){return xpr< pp<int,this_t>  >(l_,r_);}
  inline friend xpr< pp<Real,this_t> > operator+(const Real& l_, const this_t& r_){return xpr< pp<Real,this_t> >(l_,r_);}
  inline friend xpr< pp<Cplx,this_t> > operator+(const Cplx& l_, const this_t& r_){return xpr< pp<Cplx,this_t> >(l_,r_);}
  
  //==== Soustraction
  template <class r_t> xpr< mm<this_t,r_t> > operator-(const r_t& r_) const {return xpr< mm<this_t,r_t> >(*this,r_);}
  xpr< mm<this_t, int>  > operator-(const int&  r_) const {return xpr< mm<this_t, int>  >(*this,r_);}
  xpr< mm<this_t, Real> > operator-(const Real& r_) const {return xpr< mm<this_t, Real> >(*this,r_);}
  xpr< mm<this_t, Cplx> > operator-(const Cplx& r_) const {return xpr< mm<this_t, Cplx> >(*this,r_);}   

  inline friend xpr< mm<int,this_t>  > operator-(const int&  l_, const this_t& r_){return xpr< mm<int,this_t>  >(l_,r_);}
  inline friend xpr< mm<Real,this_t> > operator-(const Real& l_, const this_t& r_){return xpr< mm<Real,this_t> >(l_,r_);}
  inline friend xpr< mm<Cplx,this_t> > operator-(const Cplx& l_, const this_t& r_){return xpr< mm<Cplx,this_t> >(l_,r_);}
  
  //==== Multiplication
  template <class r_t> xpr< tt<this_t,r_t> > operator*(const r_t& r_) const {return xpr< tt<this_t,r_t> >(*this,r_);}
  xpr< tt<this_t, int>  > operator*(const int&  r_) const {return xpr< tt<this_t, int>  >(*this,r_);}
  xpr< tt<this_t, Real> > operator*(const Real& r_) const {return xpr< tt<this_t, Real> >(*this,r_);}
  xpr< tt<this_t, Cplx> > operator*(const Cplx& r_) const {return xpr< tt<this_t, Cplx> >(*this,r_);}   
  
  inline friend xpr< tt<int,this_t>  > operator*(const int&  l_, const this_t& r_){return xpr< tt<int,this_t>  >(l_,r_);}
  inline friend xpr< tt<Real,this_t> > operator*(const Real& l_, const this_t& r_){return xpr< tt<Real,this_t> >(l_,r_);}
  inline friend xpr< tt<Cplx,this_t> > operator*(const Cplx& l_, const this_t& r_){return xpr< tt<Cplx,this_t> >(l_,r_);}


  //===== Produit scalaire
  template <class r_t>typename resop<this_t,r_t>::type operator,(const r_t& r_) const {
    assert(nr==r_t::nr); return dprod<this_t,r_t,nr>::apply(*this,r_); }

  
};


//==========================//
//         Sub-Matrix       //
//==========================//

template <class m_t, class ro_t, class co_t>
class submat{
  
public:  
  
  // Constantes
  static const int nr  = ro_t::nr;
  static const int nc  = co_t::nr;
  static const int dim = nr*nc;
  
  typedef typename m_t::v_t      v_t;
  typedef array<nr,v_t>          a_t;
  typedef submat<m_t,ro_t,co_t>  this_t;
  
private:

  m_t& m;
  const ro_t& I; 
  const co_t& J;

public:
  
  submat<m_t,ro_t,co_t>(m_t& m0, const ro_t& I0, const co_t& J0): m(m0), I(I0), J(J0){};
  
  template <class r_t> void operator =(const r_t& r_){ assign_dbloop<this_t,r_t>::apply(*this,r_);      }
  
  template <class r_t> void operator+=(const r_t& r_){ plus_assign_dbloop<this_t,r_t>::apply(*this,r_); }  
  
  v_t& operator()(const int& j, const int& k){ return m(I[j],J[k]); }

  const v_t& operator()(const int& j, const int& k) const { return m(I[j],J[k]); }
  
  inline friend ostream& operator<<(ostream& os, const this_t& r_){ ostream_dbloop<this_t>::apply(os,r_); return os;}
  

  //==== Addition
  template <class r_t> xpr< pp<this_t,r_t> > operator+(const r_t& r_) const {return xpr< pp<this_t,r_t> >(*this,r_);}
  xpr< pp<this_t, int>  > operator+(const int&  r_) const {return xpr< pp<this_t, int>  >(*this,r_);}
  xpr< pp<this_t, Real> > operator+(const Real& r_) const {return xpr< pp<this_t, Real> >(*this,r_);}
  xpr< pp<this_t, Cplx> > operator+(const Cplx& r_) const {return xpr< pp<this_t, Cplx> >(*this,r_);}   

  inline friend xpr< pp<int,this_t>  > operator+(const int&  l_, const this_t& r_){return xpr< pp<int,this_t>  >(l_,r_);}
  inline friend xpr< pp<Real,this_t> > operator+(const Real& l_, const this_t& r_){return xpr< pp<Real,this_t> >(l_,r_);}
  inline friend xpr< pp<Cplx,this_t> > operator+(const Cplx& l_, const this_t& r_){return xpr< pp<Cplx,this_t> >(l_,r_);}
  
  //==== Soustraction
  template <class r_t> xpr< mm<this_t,r_t> > operator-(const r_t& r_) const {return xpr< mm<this_t,r_t> >(*this,r_);}
  xpr< mm<this_t, int>  > operator-(const int&  r_) const {return xpr< mm<this_t, int>  >(*this,r_);}
  xpr< mm<this_t, Real> > operator-(const Real& r_) const {return xpr< mm<this_t, Real> >(*this,r_);}
  xpr< mm<this_t, Cplx> > operator-(const Cplx& r_) const {return xpr< mm<this_t, Cplx> >(*this,r_);}   

  inline friend xpr< mm<int,this_t>  > operator-(const int&  l_, const this_t& r_){return xpr< mm<int,this_t>  >(l_,r_);}
  inline friend xpr< mm<Real,this_t> > operator-(const Real& l_, const this_t& r_){return xpr< mm<Real,this_t> >(l_,r_);}
  inline friend xpr< mm<Cplx,this_t> > operator-(const Cplx& l_, const this_t& r_){return xpr< mm<Cplx,this_t> >(l_,r_);}
  
  //==== Multiplication
  template <class r_t> xpr< tt<this_t,r_t> > operator*(const r_t& r_) const {return xpr< tt<this_t,r_t> >(*this,r_);}
  xpr< tt<this_t, int>  > operator*(const int&  r_) const {return xpr< tt<this_t, int>  >(*this,r_);}
  xpr< tt<this_t, Real> > operator*(const Real& r_) const {return xpr< tt<this_t, Real> >(*this,r_);}
  xpr< tt<this_t, Cplx> > operator*(const Cplx& r_) const {return xpr< tt<this_t, Cplx> >(*this,r_);}   
  
  inline friend xpr< tt<int,this_t>  > operator*(const int&  l_, const this_t& r_){return xpr< tt<int,this_t>  >(l_,r_);}
  inline friend xpr< tt<Real,this_t> > operator*(const Real& l_, const this_t& r_){return xpr< tt<Real,this_t> >(l_,r_);}
  inline friend xpr< tt<Cplx,this_t> > operator*(const Cplx& l_, const this_t& r_){return xpr< tt<Cplx,this_t> >(l_,r_);}

  
};


//==========================//
//           Matrix         //
//==========================//

template <int Nr, int Nc, class V_t>
  class mat{
  
 public:
  static const int dim =  Nr*Nc;
  static const int nr  =     Nr;
  static const int nc  =     Nc;
  
  typedef V_t                v_t;
  typedef array<nr,v_t>      a_t;
  typedef mat<nr,nc,v_t>  this_t;
  
private:
  
  v_t  v_[dim];
  
public:
    
  mat<Nr,Nc,V_t>(){ construct_dbloop<this_t>::apply(*this); }
  
  template <class r_t> mat<Nr,Nc,V_t>(const r_t& r_){ assign_dbloop<this_t,r_t>::apply(*this,r_); } 
  
  template <class r_t> void operator=(const r_t& r_){ assign_dbloop<this_t,r_t>::apply(*this,r_); }

  template <class r_t> void operator+=(const r_t& r_){ plus_assign_dbloop<this_t,r_t>::apply(*this,r_); } 
  
  v_t& operator()(const int& j, const int& k){ return v_[j+k*nr]; }
  
  const v_t& operator()(const int& j, const int& k) const { return v_[j+k*nr]; }

  template <class ro_t, class co_t> submat<this_t,ro_t,co_t> 
  operator()(const ro_t& I, const co_t& J){ return submat<this_t,ro_t,co_t>(*this,I,J);}
  
  template <class ro_t, class co_t> submat< const this_t,ro_t,co_t> 
  operator()(const ro_t& I, const co_t& J) const { return submat<this_t,ro_t,co_t>(*this,I,J);}
  
  inline friend ostream& operator<<(ostream& os, const this_t& m){ ostream_dbloop<this_t>::apply(os,m); return os;}

  //==== Addition
  template <class r_t> xpr< pp<this_t,r_t> > operator+(const r_t& r_) const {return xpr< pp<this_t,r_t> >(*this,r_);}
  xpr< pp<this_t, int>  > operator+(const int&  r_) const {return xpr< pp<this_t, int>  >(*this,r_);}
  xpr< pp<this_t, Real> > operator+(const Real& r_) const {return xpr< pp<this_t, Real> >(*this,r_);}
  xpr< pp<this_t, Cplx> > operator+(const Cplx& r_) const {return xpr< pp<this_t, Cplx> >(*this,r_);}   

  inline friend xpr< pp<int,this_t>  > operator+(const int&  l_, const this_t& r_){return xpr< pp<int,this_t>  >(l_,r_);}
  inline friend xpr< pp<Real,this_t> > operator+(const Real& l_, const this_t& r_){return xpr< pp<Real,this_t> >(l_,r_);}
  inline friend xpr< pp<Cplx,this_t> > operator+(const Cplx& l_, const this_t& r_){return xpr< pp<Cplx,this_t> >(l_,r_);}
  
  //==== Soustraction
  template <class r_t> xpr< mm<this_t,r_t> > operator-(const r_t& r_) const {return xpr< mm<this_t,r_t> >(*this,r_);}
  xpr< mm<this_t, int>  > operator-(const int&  r_) const {return xpr< mm<this_t, int>  >(*this,r_);}
  xpr< mm<this_t, Real> > operator-(const Real& r_) const {return xpr< mm<this_t, Real> >(*this,r_);}
  xpr< mm<this_t, Cplx> > operator-(const Cplx& r_) const {return xpr< mm<this_t, Cplx> >(*this,r_);}   

  inline friend xpr< mm<int,this_t>  > operator-(const int&  l_, const this_t& r_){return xpr< mm<int,this_t>  >(l_,r_);}
  inline friend xpr< mm<Real,this_t> > operator-(const Real& l_, const this_t& r_){return xpr< mm<Real,this_t> >(l_,r_);}
  inline friend xpr< mm<Cplx,this_t> > operator-(const Cplx& l_, const this_t& r_){return xpr< mm<Cplx,this_t> >(l_,r_);}
  
  //==== Multiplication
  template <class r_t> xpr< tt<this_t,r_t> > operator*(const r_t& r_) const {return xpr< tt<this_t,r_t> >(*this,r_);}
  xpr< tt<this_t, int>  > operator*(const int&  r_) const {return xpr< tt<this_t, int>  >(*this,r_);}
  xpr< tt<this_t, Real> > operator*(const Real& r_) const {return xpr< tt<this_t, Real> >(*this,r_);}
  xpr< tt<this_t, Cplx> > operator*(const Cplx& r_) const {return xpr< tt<this_t, Cplx> >(*this,r_);}   
  
  inline friend xpr< tt<int,this_t>  > operator*(const int&  l_, const this_t& r_){return xpr< tt<int,this_t>  >(l_,r_);}
  inline friend xpr< tt<Real,this_t> > operator*(const Real& l_, const this_t& r_){return xpr< tt<Real,this_t> >(l_,r_);}
  inline friend xpr< tt<Cplx,this_t> > operator*(const Cplx& l_, const this_t& r_){return xpr< tt<Cplx,this_t> >(l_,r_);}
    
};


//%%%%%%%%%%%%%%%%%%%%%%%%%%//
//    Produit scalaire      //
//%%%%%%%%%%%%%%%%%%%%%%%%%%//


template <class a_t> inline typename a_t::v_t norm2(const a_t& a_){return sqrt( (a_,a_) ); }
template <class a_t> inline void normalize(a_t& a_){a_ = (1./norm2(a_))*a_;}


//%%%%%%%%%%%%%%%%%%%%%%%%%%//
//  Inversion de matrices   //
//%%%%%%%%%%%%%%%%%%%%%%%%%%//

template<class m_t>
typename m_t::v_t det(const m_t& M){
    
  if(m_t::nr== 2){
    return M(0,0)*M(1,1)-M(1,0)*M(0,1);}

  if(m_t::nr== 3){
    return M(0,0)*( M(1,1)*M(2,2) - M(2,1)*M(1,2) )
      + M(0,1)*( M(1,0)*M(2,2)-M(2,0)*M(1,2) ) 
      + M(0,2)*( M(1,0)*M(2,1)-M(2,0)*M(1,1) ); }
  
  if(m_t::nr > 3){
    cout << "matrice trop grosse" << endl;
    exit(EXIT_FAILURE);}    
}


template<class m_t>
mat<m_t::nr,m_t::nc,typename m_t::v_t>  
  inv(const m_t& M){

  mat<m_t::nr,m_t::nc,typename m_t::v_t> R;
  if(m_t::nr== 2){
    typename m_t::v_t Delta = M(0,0)*M(1,1)-M(1,0)*M(0,1);
    R(0,0) =  M(1,1)/Delta;
    R(1,1) =  M(0,0)/Delta;
    R(0,1) = -M(0,1)/Delta;
    R(1,0) = -M(1,0)/Delta; 
    return R; }
  
  if(m_t::nr== 3){
    typename m_t::v_t Delta = M(0,0)*( M(1,1)*M(2,2) - M(2,1)*M(1,2) )
      + M(0,1)*( M(1,0)*M(2,2)-M(2,0)*M(1,2) ) 
      + M(0,2)*( M(1,0)*M(2,1)-M(2,0)*M(1,1) );
    R(0,0) =  ( M(1,1)*M(2,2) - M(2,1)*M(1,2) )/Delta;
    R(1,0) = -( M(1,0)*M(2,2) - M(2,0)*M(1,2) )/Delta;
    R(2,0) =  ( M(1,0)*M(2,1) - M(2,0)*M(1,1) )/Delta;
    R(0,1) = -( M(0,1)*M(2,2) - M(2,1)*M(0,2) )/Delta;
    R(1,1) =  ( M(0,0)*M(2,2) - M(2,0)*M(0,2) )/Delta;  
    R(2,1) = -( M(0,0)*M(2,1) - M(2,0)*M(0,1) )/Delta;
    R(0,2) =  ( M(0,1)*M(1,2) - M(1,1)*M(0,2) )/Delta;
    R(1,2) = -( M(0,0)*M(1,2) - M(1,0)*M(0,2) )/Delta;
    R(2,2) =  ( M(0,0)*M(1,1) - M(1,0)*M(0,1) )/Delta;
    return R; }

  if(m_t::nr > 3){
    cout << "matrice trop grosse" << endl;
    exit(EXIT_FAILURE);}

}


template<int nr, int nc, class v_t>
mat<nc,nr,v_t> tr(const mat<nr,nc,v_t>& M){
  mat<nc,nr,v_t> R;
  for(int j=0; j<nr; j++){
    for(int k=0; k<nc; k++){
      R(k,j) = M(j,k);}}
  return R;}


//%%%%%%%%%%%%%%%%%%%%%//
//    Combinatoire     //
//%%%%%%%%%%%%%%%%%%%%%//

template <class array_t>
inline void swap(array_t& ar, const int& j, const int& k){
  typename array_t::v_t val = ar[j];
  ar[j] = ar[k]; ar[k] = val;}

template <class array_t>
inline void sort(array_t& ar){
  const int N = size(ar); 
  for(int j=1;j<N;j++){
    for(int k=0;k<N-j;k++){
      if(ar[k]>ar[k+1]){swap(ar,k,k+1);}
    }
  }
}

template <int dim>
inline void init(array<dim,int>& I){
  ascending_loop< array<dim,int> >::apply(I);  }

//%%%%%%%%%%%%%%%%%%%%%%%%%%//
//    Produit vectoriel     //
//%%%%%%%%%%%%%%%%%%%%%%%%%%//

template <class l_t, class r_t>
  inline array<3,typename resop<l_t,r_t>::type> vprod(const l_t& u, const r_t& v){
  array<3,typename resop<l_t,r_t>::type> w;
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
  return w; }


//%%%%%%%%%%%%%%%%%%%%%%%%%%//
//   Definitions de type    //
//%%%%%%%%%%%%%%%%%%%%%%%%%%//

typedef array<1,   int>     N1;
typedef array<2,   int>     N2;
typedef array<3,   int>     N3;
typedef array<4,   int>     N4;
typedef array<5,   int>     N5;
typedef array<6,   int>     N6;
typedef array<7,   int>     N7;
typedef array<8,   int>     N8;
typedef array<9,   int>     N9;
typedef array<10,  int>     N10;

typedef array<2,   Real>    R2;
typedef array<3,   Real>    R3;
typedef array<4,   Real>    R4;
typedef array<5,   Real>    R5;
typedef array<6,   Real>    R6;
typedef array<7,   Real>    R7;
typedef array<8,   Real>    R8;
typedef array<9,   Real>    R9;
typedef array<10,  Real>    R10;

typedef array<2,   Cplx>    C2;
typedef array<3,   Cplx>    C3;
typedef array<4,   Cplx>    C4;
typedef array<5,   Cplx>    C5;
typedef array<6,   Cplx>    C6;
typedef array<7,   Cplx>    C7;
typedef array<8,   Cplx>    C8;
typedef array<9,   Cplx>    C9;
typedef array<10,  Cplx>    C10;

typedef mat<2,2,   Real>    R2x2;
typedef mat<3,3,   Real>    R3x3;
typedef mat<4,4,   Real>    R4x4;
typedef mat<5,5,   Real>    R5x5;
typedef mat<6,6,   Real>    R6x6;
typedef mat<10,10, Real>    R10x10;
typedef mat<3,1,   Real>    R3x1;
typedef mat<3,2,   Real>    R3x2;
typedef mat<1,3,   Real>    R1x3;
typedef mat<2,3,   Real>    R2x3;
typedef mat<5,2,   Real>    R5x2;
typedef mat<5,3,   Real>    R5x3;

typedef mat<2,2,   Cplx>    C2x2;
typedef mat<3,3,   Cplx>    C3x3;
typedef mat<4,4,   Cplx>    C4x4;
typedef mat<6,6,   Cplx>    C6x6;
typedef mat<10,10, Cplx>    C10x10;
typedef mat<3,2,   Cplx>    C3x2;
typedef mat<2,3,   Cplx>    C2x3;

typedef array<2,N2>        _2xN2;
typedef array<3,N2>        _3xN2;
typedef array<4,N2>        _4xN2;
typedef array<2,N3>        _2xN3;
typedef array<3,N3>        _3xN3;
typedef array<4,N3>        _4xN3;
typedef array<2,N4>        _2xN4;
typedef array<3,N4>        _3xN4;
typedef array<4,N4>        _4xN4;

typedef array<2,R2>        _2xR2;
typedef array<3,R2>        _3xR2;
typedef array<4,R2>        _4xR2;
typedef array<2,R3>        _2xR3;
typedef array<3,R3>        _3xR3;
typedef array<4,R3>        _4xR3;
typedef array<2,R4>        _2xR4;
typedef array<3,R4>        _3xR4;
typedef array<4,R4>        _4xR4;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//  Construction rapide de petites matrices    //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//


mat<3,2,Real> matR3x2(const R3& X0, const R3& X1){
  R3x2 M;
  for(int j=0; j<3; j++){
    M(j,0) = X0[j];
    M(j,1) = X1[j]; }
  return M;
}


template <class x_t>
mat<x_t::nr,1, typename x_t::v_t> mat_(const x_t& x0){
  mat<x_t::nr,1, typename x_t::v_t> M;
  for(int j=0; j<x_t::nr; j++){M(j,0) = x0[j];}
  return M;
}

template <class x_t>
mat<x_t::nr,2, typename x_t::v_t> mat_(const x_t& x0, const x_t& x1){
  mat<x_t::nr,2, typename x_t::v_t> M;
  for(int j=0; j<x_t::nr; j++){
    M(j,0) = x0[j];
    M(j,1) = x1[j];
  }
  return M;
}

template <class x_t>
mat<x_t::nr,3, typename x_t::v_t> mat_(const x_t& x0, const x_t& x1, const x_t& x2){
  mat<x_t::nr,3, typename x_t::v_t> M;
  for(int j=0; j<x_t::nr; j++){
    M(j,0) = x0[j];
    M(j,1) = x1[j];
    M(j,2) = x2[j];
  }
  return M;
}






//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//   Classe vecteur grande taille    //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

template <class T> class vect {
  
 public:
  typedef T        v_t;
  typedef vect<T>  this_t;
  
 private:
  vector<T> data;
  
  // Interdit de recopier un vect
  vect<T>(const vect<T>&);
  
 public:
  
  // Constructeur par defaut
  vect<T>(){};

  // Operateurs d'acces
  T& operator[](const int& j){return data[j];}
  
  const T& operator[](const int& j) const {return data[j];}

  T& back(){ return data.back();}
  
  const T& back() const { return data.back();}

  friend vector<T>& vector_of(this_t& v){return v.data;}

  friend const vector<T>& vector_of(const this_t& v){return v.data;}
  
  void push_back(const T& t){data.push_back(t);}
  
  friend int size  (const this_t& v){return v.data.size(); }

  friend void fill (this_t& v,const T& d){
    for(int j=0; j<size(v); j++){v[j]=d;}}
  
  friend void clear (this_t& v){v.data.clear();}
  
  friend void resize(this_t& v, const int& N){v.data.resize(N);}
  
  friend ostream& operator<<(ostream& os, const this_t& v){
    for(int j=0; j<size(v); j++){os << v[j] << endl;} return os;}
  
  template <class i_t> subarray<this_t,i_t> operator[] (const i_t& i_){
    return subarray<this_t,i_t>(*this,i_);}
  
  template <class i_t> subarray<const this_t,i_t> operator[] (const i_t& i_) const {
    return subarray<const this_t,i_t>(*this,i_);}
  
};




#endif
