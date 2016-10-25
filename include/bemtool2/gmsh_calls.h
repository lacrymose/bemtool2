#ifndef GMSH_CALL_H
#define GMSH_CALL_H

#include <string>
#include "calculus.h"
#include "user.h"


////=============================================================////
////===========================  Circle =========================////
////=============================================================////
void gmsh_circle(std::string mesh_name, Real R, Real lc);

////=============================================================////
////===========================  Clean ==========================////
////=============================================================////
void gmsh_clean(std::string mesh_name);

#endif

