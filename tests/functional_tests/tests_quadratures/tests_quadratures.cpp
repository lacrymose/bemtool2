#include <bemtool2/mesh.h>
#include <bemtool2/kernel.h>
#include <bemtool2/dof.h>

using namespace bemtool;
////=============================================================================////
////========================= quad1D_potential ==================================////
////=============================================================================////

void quad1D_potential(){

	////================ Construction d'un maillage =================////
	R3 A;A[0]=0;A[1]=0;A[2]=0;
	R3 B;B[0]=1;B[1]=0;B[2]=0;
	R3 C;C[0]=1;C[1]=1;C[2]=0;

	vect<R3> pts;pts.push_back(A);pts.push_back(B);

	geometry geom;
	load_node_hand(geom,pts);

	bemtool::array<2,int> I; I[0]=0; I[1]=1;

	mesh_1D mesh(geom);
	load_elt_hand(mesh,I);

	nrml_1D n_(mesh);
	Real kappa=1;

	////================= Calcul de la quadrature ===================////

	potential<P1_1D,CST_2D> CSTop(kappa,n_);
	mat<1,2,Cplx > CST = CSTop(C,mesh[0]);



	////========================= Comparaison =======================////
	// 0.5 attendu
	std::cout <<"2D - calcul par elt    "<< CST << std::endl;
	std::cout <<"2D - calcul par noeuds "<< CSTop(C,0)<<"\t"<< CSTop(C,1) << std::endl;
}

////=============================================================================////
////========================= quad2D_potential ==================================////
////=============================================================================////

void quad2D_potential(){

	////================ Construction d'un maillage =================////
	R3 A;A[0]=1;A[1]=0;A[2]=0;
	R3 B;B[0]=1;B[1]=1;B[2]=0;
	R3 C;C[0]=0;C[1]=0;C[2]=0;
	R3 D;D[0]=0;D[1]=0;D[2]=1;

	vect<R3> pts;pts.push_back(A);pts.push_back(B);pts.push_back(C);

	geometry geom;
	load_node_hand(geom,pts);
	// load_node_gmsh(geom,"test");
	bemtool::array<3,int> I; I[0]=0; I[1]=1;I[2]=2;

	mesh_2D mesh(geom);
	load_elt_hand(mesh,I);
	// load_elt_gmsh(mesh,0);

	nrml_2D n_(mesh);
	Real kappa=1;
	////================= Calcul de la quadrature ===================////
	potential<P1_2D,CST_3D> CSTop(kappa,n_);
	mat<1,3,Cplx > CST = CSTop(D,mesh[0]);

	////========================= Comparaison =======================////
	// 0.166667 attendu
	std::cout <<"3D - calcul par elt    "<< CSTop(D,mesh[0]) << std::endl;
	std::cout <<"3D - calcul par noeuds "<< CSTop(D,0) <<"\t"<< CSTop(D,1)<<"\t"<<CSTop(D,2)<< std::endl;

}

////=============================================================================////
////========================= quad1D_bem ==================================////
////=============================================================================////

void quad1D_bem(){

	////================ Construction d'un maillage =================////
	R3 A;A[0]=0;A[1]=0;A[2]=0;
	R3 B;B[0]=1;B[1]=0;B[2]=0;
	R3 C;C[0]=1;C[1]=1;C[2]=0;

	vect<R3> pts;pts.push_back(A);pts.push_back(B);

	geometry geom;
	load_node_hand(geom,pts);

	bemtool::array<2,int> I; I[0]=0; I[1]=1;

	mesh_1D mesh(geom);
	load_elt_hand(mesh,I);

	nrml_1D n_(mesh);
	Real kappa=1;

	////================= Calcul de la quadrature ===================////
	std::vector<int> dofs_I {0,1};
	std::vector<int> dofs_J {0,1};

	bem<P1_1D,P1_1D,CST_2D> CSTop(kappa,n_,n_);
	mat<2,2,Cplx > CST1 = CSTop(mesh[0],mesh[0]);
	mat<2,2,Cplx > CST2;
	CSTop(dofs_I,dofs_J,CST2);



	////========================= Comparaison =======================////
	// 0.5 attendu


	std::cout <<"2D - calcul par elt    "<< CST1 << std::endl;
	std::cout <<"2D - calcul par noeuds "<< CST2 << std::endl;
}

////=============================================================================////
////========================= quad2D_bem ==================================////
////=============================================================================////

void quad2D_bem(){

	////================ Construction d'un maillage =================////
	R3 A;A[0]=1;A[1]=0;A[2]=0;
	R3 B;B[0]=1;B[1]=1;B[2]=0;
	R3 C;C[0]=0;C[1]=0;C[2]=0;

	vect<R3> pts;pts.push_back(A);pts.push_back(B);pts.push_back(C);

	geometry geom;
	load_node_hand(geom,pts);
	// load_node_gmsh(geom,"test");
	bemtool::array<3,int> I; I[0]=0; I[1]=1;I[2]=2;

	mesh_2D mesh(geom);
	load_elt_hand(mesh,I);
	// load_elt_gmsh(mesh,0);

	nrml_2D n_(mesh);
	Real kappa=1;
	////================= Calcul de la quadrature ===================////
	bem<P1_2D,P1_2D,CST_3D> CSTop(kappa,n_,n_);
	mat<3,3,Cplx > CST = CSTop(mesh[0],mesh[0]);

	////========================= Comparaison =======================////
	// 0.166667 attendu
	std::cout <<"3D - calcul par elt    "<< CST << std::endl;
	std::cout <<"3D - calcul pasr noeuds "<< CSTop(0,0) <<"\t"<< CSTop(0,1)<<"\t"<<CSTop(0,2)<< std::endl;

}

////=============================================================================////
////============================= Main ==========================================////
////=============================================================================////

int main(){

    quad1D_potential();
		quad1D_bem();

    quad2D_potential();
		quad2D_bem();

}
