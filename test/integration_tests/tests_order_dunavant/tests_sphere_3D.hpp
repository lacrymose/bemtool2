#ifndef TESTSSPHERE_H
#define TESTSSPHERE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "bemtool2/tools.h"

using namespace bemtool;

void first_kind_dirichlet_3D  (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
//a void first_kind_neumann_2D    (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
//a void second_kind_dirichlet_2D (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
//a void second_kind_neumann_2D   (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
//a 
void fourier_harmonic_3D      (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
//a void plane_wave_harmonics_2D  (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
//a 
//a void champs_rayonne_2D        (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
#endif
