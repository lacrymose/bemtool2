#ifndef TESTSCIRCLE_H
#define TESTSCIRCLE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "bemtool2/mesh.h"
#include "bemtool2/normal.h"
#include "bemtool2/dof.h"
#include "bemtool2/femP1.h"
#include "bemtool2/kernel.h"
#include "bemtool2/misc.h"
// #include "gmm_wrap.h"
#include "bemtool2/gmsh_calls.h"
#include "bemtool2/gmm_wrap.h"
using namespace bemtool;
void first_kind_dirichlet_2D  (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
void first_kind_neumann_2D    (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
void second_kind_dirichlet_2D (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
void second_kind_neumann_2D   (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);

void fourier_harmonic_2D      (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
void plane_wave_harmonics_2D  (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);

#endif
