#include "tests_quadratures.hpp"

void quad1D_chps_rayonne(){
	
	
	R3 A;A[0]=0;A[1]=0;A[2]=0;
	R3 B;B[0]=0;B[1]=1;B[2]=0;
	R3 C;C[0]=1;C[1]=0;C[2]=0;
	
	vect<R3> pts;pts.push_back(A);pts.push_back(B);
	
	geometry geom;
	load_node_hand(geom,pts);
	
	bemtool::array<2,int> I; I[0]=0; I[1]=1;

	mesh_1D mesh(geom);
	load_elt_hand(mesh,I);
	
	nrml_1D n_(mesh);
	Real kappa=1;
	chps_rayonne<P1_1D,CST_2D> CSTop(kappa,n_);
	mat<1,2,Cplx > CST = CSTop(C,mesh[0]);
	
	std::cout << CST << std::endl;
}