#ifndef TESTSCIRCLE_H
#define TESTSCIRCLE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "bemtool2/tools.h"

void first_kind_dirichlet_2D  (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
void first_kind_neumann_2D    (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
void second_kind_dirichlet_2D (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
void second_kind_neumann_2D   (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);

void fourier_harmonic_2D      (std::vector<Real> harmonics, Real lc, Real R, std::string output_name,int verbose=0);
void plane_wave_harmonics_2D  (std::vector<Real> harmonics, Real lc, Real R, std::string output_name,int verbose=0);
#endif
