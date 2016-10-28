#include <iostream>
#include <fstream>
#include "tools.h"

using namespace std;



int main(){
////=============================================================////
////=======================  Mesh loading  ======================////
////=============================================================////

	load_node("../test/data/two_circles.msh");
	const vect<R3>& node = get_node();  
	int nb_node = size(node);
	
	mesh_1D Gamma[2];
	load(Gamma[0],0);
	load(Gamma[1],1);     

	mesh_1D Omega;  
	Omega += Gamma[0];
	Omega += Gamma[1];



////=============================================================////
////================== trace de la frontiere ====================////
////=============================================================////
	
// 	for (int j=0;j<2;j++){
// 		ofstream file; file.open(("mesh_1_3_0.1_gamma_"+NbrToStr(j)+".txt").c_str());
// 		for(int k=0; k<nb_elt(Gamma[j]); k++){
// 			elt_1D seg = Gamma[j][k];
// 			file << seg[0][0] <<" "<< seg[0][1] << endl;
// 			file << endl;
// 		}
// 		file.close();
// 	}
	
	 
////=============================================================////
////================== Calcul de la normale =====================////
////=============================================================////
	
	nrml_1D n_(Omega);
	
	for (int j=0;j < nb_elt(Omega);j++){
		elt_1D seg = Omega[j];
		R3 G;
		G[0]= 0.5 * ( seg[0][0] + seg[1][0] );
		G[1]= 0.5 * ( seg[0][1] + seg[1][1] );
		G[2]= 0;
		
		if (((G,n_[j])<0 & j<nb_elt(Gamma[0])) || ((G,n_[j])>0 & j>=nb_elt(Gamma[0]))){
			exit(1);
		}

	}
////=============================================================////



}
