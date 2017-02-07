
#include <bemtool2/mesh.h>
#include <bemtool2/normal.h>
#include <bemtool2/gmsh_calls.h>
#include <htool/view.hpp>

namespace bemtool{
	
void attach_ui(htool::Scene& s) {
	
	// Variables
	static Real lc= 0.1;
	static Real R = 1 ;
	enum enumeration {Disc,Ball,Sphere};
	static enumeration enumval = Disc;
	bool enabled = true;
	
	// Scene (screen inside)
	htool::statics& gv = htool::Scene::gv;
	
	// FormeHelper to hepl the construction of a window
	nanogui::FormHelper *gui = new nanogui::FormHelper(gv.screen);
	
	// Window
	nanogui::ref<nanogui::Window> nanoguiWindow = gui->addWindow(Eigen::Vector2i(10, 10), "Bemtool - Disc");
	nanoguiWindow->setPosition(Eigen::Vector2i(500, 200));
	
	
	// Parts of the window 
	gui->addVariable("Mesh", enumval, enabled)
           ->setItems({"Disc", "Ball","Sphere"});
	gui->addGroup("Parameters");
	gui->addVariable("Mesh size", lc)->setSpinnable(true);
	gui->addVariable("Radius", R)->setSpinnable(true);
	gui->addButton("Create mesh",[&]{
   	  	if (gv.active_project == NULL)
			std::cerr << "No active project" << std::endl;
		else {
			

				geometry vol;
				
				switch(enumval){
					case Disc:
					{
						gmsh_disc  ("disc",R,lc,0);
						load_node_gmsh(vol ,"disc");
						std::cout << "Mesh of a disc computed with radius "<< R<<" and mesh size "<<lc << std::endl;
					}
					case Sphere:
					{
						gmsh_disc  ("sphere",R,lc,0);
						load_node_gmsh(vol ,"sphere");
						std::cout << "Mesh of a sphere computed with radius "<< R<<" and mesh size "<<lc << std::endl;
					}
					case Ball:
					{
						gmsh_ball  ("ball",R,lc,0);
						load_node_gmsh(vol ,"ball");
						std::cout << "Mesh of a ball computed with radius "<< R<<" and mesh size "<<lc << std::endl;
					}
				}
				
				mesh_2D Vol(vol);
				load_elt_gmsh(Vol,0);
				
				const vect<R3>& node = get_node(vol);
				int nbnode  = size(node);
				int nbelt = nb_elt(Vol);
					
				std::vector<htool::R3> X(nbnode);
				std::vector<int> NbPts(nbelt);
				std::vector<htool::N4>  Elts(nbelt);
				std::vector<htool::R3> Normals(nbelt);


				for (int i=0;i<nbnode;i++){
					htool::R3 pt;pt[0]=node[i][0];pt[1]=node[i][1];pt[2]=node[i][2];
					X[i]=pt;
				}

				switch(enumval){
					case Disc:
					{
						for (int i=0;i<nbelt;i++){
							Elts[i][0]=&Vol[i][0]-&node[0];
							Elts[i][1]=&Vol[i][1]-&node[0];
							Elts[i][2]=&Vol[i][2]-&node[0];
							Elts[i][3]=1;
					
							NbPts[i]=3;
					
							Normals[i][0] = 0;
							Normals[i][1] = 0;
							Normals[i][2] = 1;
						}
					}
					case Sphere:
					{
						for (int i=0;i<nbelt;i++){					
							Elts[i][0]=&Vol[i][0]-&node[0];
							Elts[i][1]=&Vol[i][1]-&node[0];
							Elts[i][2]=&Vol[i][2]-&node[0];
							Elts[i][3]=1;
					
							NbPts[i]=3;
					
							Normals[i][0] = 0;
							Normals[i][1] = 0;
							Normals[i][2] = 1;
						}
					}
					case Ball:
					{
						for (int i=0;i<nbelt;i++){					
							Elts[i][0]=&Vol[i][0]-&node[0];
							Elts[i][1]=&Vol[i][1]-&node[0];
							Elts[i][2]=&Vol[i][2]-&node[0];
							Elts[i][3]=1;
					
							NbPts[i]=3;
							htool::R3 v1 = X[Elts[i][1]] - X[Elts[i][0]];
							htool::R3 v2 = X[Elts[i][2]] - X[Elts[i][0]];
							htool::R3 v3 =v1^v2;
							Normals[i][0] = v3[0];
							Normals[i][1] = v3[1];
							Normals[i][2] = v3[2];
						}
					}
				}
				
// 				gv.active_project->set_ctrs(X);
// 				gv.active_project->set_rays(Rays);
		
// 				gv.screen->performLayout();
				
				htool::GLMesh m(X,Elts,NbPts,Normals);
				s.set_mesh(m);
		}
   });
	

	gv.screen->performLayout();

}

} // namespace