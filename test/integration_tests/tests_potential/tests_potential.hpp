#ifndef TESTSPOTENTIAL_H
#define TESTSPOTENTIAL_H

#include "bemtool2/tools.h"

void potential_elt_2D        (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);
void potential_node_2D       (std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0);


#endif
