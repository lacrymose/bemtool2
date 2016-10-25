
#include "bemtool2/gmsh_calls.h"

using namespace std;

////=============================================================////
////===========================  Circle =========================////
////=============================================================////
void gmsh_circle(std::string mesh_name, Real R, Real lc){
	std::cout<<"Création du maillage"<<std::endl;
	
	std::ofstream circle(mesh_name.c_str());
	if (circle.is_open()) {
		circle << "lc = "+NbrToStr(lc)+";\n";
		circle << "Point(0) = { 0 , 0 , 0 , lc};\n";
		
		//// Droite
		circle << "Point(1) = { "+NbrToStr(R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Haut
		circle << "Point(2) = { "+NbrToStr(0)+" , "+NbrToStr(R)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Gauche
		circle << "Point(3) = { "+NbrToStr(-R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Bas
		circle << "Point(4) = { "+NbrToStr(0)+" , "+NbrToStr(-R)+" , "+NbrToStr(0)+" , lc}; \n";

		//// Droite
		circle << "Circle(1) = { "+NbrToStr(1)+" , "+NbrToStr(0)+" , "+NbrToStr(2)+"}; \n";
		//// Haut
		circle << "Circle(2) = { "+NbrToStr(2)+" , "+NbrToStr(0)+" , "+NbrToStr(3)+"}; \n";
		//// Gauche
		circle << "Circle(3) = { "+NbrToStr(3)+" , "+NbrToStr(0)+" , "+NbrToStr(4)+"}; \n";
		//// Bas
		circle << "Circle(4) = { "+NbrToStr(4)+" , "+NbrToStr(0)+" , "+NbrToStr(1)+"}; \n";	
			
		
		circle << "Physical Line(0) ={ "+NbrToStr(1)+" , "+NbrToStr(2)+" , "+NbrToStr(3)+" , "+NbrToStr(4)+"};\n";
		
		
		circle.close();
	}
	else std::cout << "Unable to open file \n" <<std::endl;
	
	system(("gmsh -2 "+mesh_name).c_str());
}

////=============================================================////
////===========================  Clean ==========================////
////=============================================================////
void gmsh_clean(std::string mesh_name){
	system(("rm "+mesh_name+".geo").c_str());
	system(("rm "+mesh_name+".msh").c_str());
}