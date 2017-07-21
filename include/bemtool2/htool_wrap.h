#ifndef HTOOL_WRAP_H
#define HTOOL_WRAP_H

#include <htool/htool.hpp>
#include <mpi.h>
#include "tools.h"
namespace bemtool{


// template <class space, class kernel> class MyMatrix;
//
// template <class space, class kernel> const Real NormFrob (const MyMatrix<space, kernel>& A);


template <class space, class kernel> class MyPotential: public htool::IMatrix<Cplx>{
private:
	potential<space,kernel>& Op;
	const geometry& geom;

public:
	MyPotential(potential<space,kernel>& op,const geometry& geome, const int& nbpts, const int& nbdof): IMatrix(nbpts,nbdof), Op(op),geom(geome) {
	}

	Cplx get_coef(const int& j, const int& k) const{
// 		std::cout << get_node(geom,j) << std::endl;
		return Op(get_node(geom,j),k);
	}

	// friend const Real NormFrob<space,kernel> (const MyMatrix<space,kernel>& A);
};

// template <class space, class kernel> const Real NormFrob (const MyMatrix<space, kernel>& A){
//     Real norm=0;
//     for (int j=0;j<A.nr;j++){
//         for (int k=0;k<A.nc;k++){
//             norm = norm + pow(std::abs(A.get_coef(j,k)),2);
//         }
//     }
//     return sqrt(norm);
// }

}
#endif
