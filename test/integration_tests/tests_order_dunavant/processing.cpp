#include "tests_sphere_3D.hpp"

int main(){
    
    //// Rayon du cercle
    Real R= 1.;

	
    //// Finesses
    
    //std::vector<Real> finesses;finesses.push_back(0.3);finesses.push_back(0.1);finesses.push_back(0.07);finesses.push_back(0.05);
    //std::vector<Real> finesses;finesses.push_back(1);finesses.push_back(0.1);finesses.push_back(0.07);finesses.push_back(0.05);
    std::vector<Real> finesses;finesses.push_back(1); //finesses.push_back(0.4);finesses.push_back(0.3);finesses.push_back(0.2);

    //// Harmoniques
    std::vector<Real> harmoniques;harmoniques.push_back(1);harmoniques.push_back(2);harmoniques.push_back(3);harmoniques.push_back(4);harmoniques.push_back(5);

    //// Tests de convergence
    for (int i=0;i<finesses.size();i++){
        std::cout<<"Test for lc="+NbrToStr<Real>(finesses[i])<<std::endl;
//         
//         first_kind_dirichlet_2D (harmoniques, finesses[i], R, "first_kind_dirichlet_2D");
//         first_kind_neumann_2D   (harmoniques, finesses[i], R, "first_kind_neumann_2D");
//         second_kind_dirichlet_2D(harmoniques, finesses[i], R, "second_kind_dirichlet_2D");
//         second_kind_neumann_2D  (harmoniques, finesses[i], R, "second_kind_neumann_2D");
//         
        fourier_harmonic_3D     (harmoniques, finesses[i], R, "fourier_harmonics",1);
//         plane_wave_harmonics_2D (harmoniques, finesses[i], R, "plane_wave_harmonics");
// 		
// //		champs_rayonne_2D (harmoniques, finesses[i], R, "plane_wave_harmonics");
//         
//         
//         
    }
   //first_kind_dirichlet_3D(harmoniques, 0.25, R, "first_kind_dirichlet_3D", 1);
   //first_kind_dirichlet_3D(harmoniques, 0.4, R, "first_kind_dirichlet_3D", 1);
}
