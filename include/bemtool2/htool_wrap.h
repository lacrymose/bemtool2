#ifndef HTOOL_WRAP_H
#define HTOOL_WRAP_H

#include "htool.hpp"
#include "tools.h"

template <class space, class kernel> class MyMatrix: public htool::VirtualMatrix{
private:
	potential<space,kernel>& Op;
	const geometry& geom;
	
public:
	MyMatrix(potential<space,kernel>& op,const geometry& geome, const int& nbpts, const int& nbdof): Op(op),geom(geome) {
		nr = nbpts;
		nc = nbdof;
	}

	const Cplx get_coef(const int& j, const int& k) const{
// 		std::cout << get_node(geom,j) << std::endl;
		return Op(get_node(geom,j),k);
	}
};


#endif