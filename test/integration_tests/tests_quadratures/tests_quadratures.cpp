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
	potential<P1_1D,CST_2D> CSTop(kappa,n_);
	mat<1,2,Cplx > CST = CSTop(C,mesh[0]);

	// 0.5 attendu
	std::cout << CST << std::endl;
}

void quad2D_chps_rayonne(){


	R3 A;A[0]=1;A[1]=0;A[2]=0;
	R3 B;B[0]=1;B[1]=1;B[2]=0;
	R3 C;C[0]=0;C[1]=0;C[2]=0;
	R3 D;D[0]=0;D[1]=0;D[2]=1;

	vect<R3> pts;pts.push_back(A);pts.push_back(B);pts.push_back(C);

	geometry geom;
	load_node_hand(geom,pts);

	bemtool::array<3,int> I; I[0]=0; I[1]=1;I[2]=2;

	mesh_2D mesh(geom);
	load_elt_hand(mesh,I);

	nrml_2D n_(mesh);
	Real kappa=1;
	potential<P1_2D,CST_3D> CSTop(kappa,n_);
	mat<1,3,Cplx > CST = CSTop(D,mesh[0]);
	// std::cout << mesh[0] <<std::endl;
	std::cout << CSTop(D,mesh[0]) << std::endl;
}
