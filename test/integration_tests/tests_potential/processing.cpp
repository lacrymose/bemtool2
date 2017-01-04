#include "tests_potential.hpp"


int main(){

    //// Rayon du cercle
    Real R= 1.;

    //// Harmoniques
    std::vector<Real> harmoniques;harmoniques.push_back(1);harmoniques.push_back(2);harmoniques.push_back(3);harmoniques.push_back(4);harmoniques.push_back(5);

	//// Finesse
	Real finesse = 0.05;
	
	//// Tests
	potential_elt_2D (harmoniques, finesse, R, "potential_elt_2D");
	potential_node_2D (harmoniques, finesse, R, "potential_node_2D");

}