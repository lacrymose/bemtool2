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

template<typename bem>
class MyBEM: public htool::IMatrix<Cplx>{
	bem& Op;

public:
	MyBEM(bem& Op0,int nbdof0):IMatrix(nbdof0,nbdof0),Op(Op0){}

	Cplx get_coef(const int& i, const int& j)const {return Op(i,j);}

	htool::SubMatrix<Cplx> get_submatrix(const std::vector<int>& J, const std::vector<int>& K ) const{
		htool::SubMatrix<Cplx> submat(J,K);
		Op(J,K,submat);
		return submat;
	}


};
}
#endif
