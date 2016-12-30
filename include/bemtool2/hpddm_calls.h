#ifndef HPDDM_CALL_H
#define HPDDM_CALL_H

#define HPDDM_NUMBERING 'C'
#include <HPDDM.hpp>
#include "gmm_wrap.h"

struct HPDDMOperator_GMM : HPDDM::EmptyOperator<std::complex<double>> {
  gmm_dense& _A;
  HPDDMOperator_GMM(gmm_dense& A) : HPDDM::EmptyOperator<std::complex<double>>(nb_rows(A)), _A(A) { }
  void GMV(const std::complex<double>* const in, std::complex<double>* const out, const int& mu = 1) const {
    vect<Cplx> b; resize(b, _n);
    std::complex<double>* pt = &(b[0]);
    std::copy_n(in, _n, pt);
    vect<Cplx> x; resize(x, _n);
    mv_prod(x, _A, b);
    pt = &(x[0]);
    std::copy_n(pt, _n, out);
  }
  template<bool = true>
  void apply(const std::complex<double>* const in, std::complex<double>* const out, const unsigned short& mu = 1, std::complex<double>* = nullptr, const unsigned short& = 0) const {
    for(int i = 0; i < _n; ++i) {
      #if 1
      out[i] = in[i] / _A(i, i);
      #else
      out[i] = in[i];
      #endif
    }
  }
};

#endif
