#include "bemtool2/tools.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>
// #include "mpi.h"
// #include "bemtool2/htool_wrap.h"

///==================================================================================////
///============================= potential_elt_2D ===================================////
///==================================================================================////

// void potential_elt_2D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0){
// 
//     ////=======================  Mesh building  =====================////
// 
//     if (verbose>0){
//         std::cout<<"Construction du maillage"<<std::endl;
//     }
//     gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc,verbose);
//     gmsh_disc  (("disc_"+NbrToStr(lc)).c_str(),R*0.9,lc,verbose);
// 
// 
//     ////=======================  Mesh loading  ======================////
// 
//     if (verbose>0){
//         std::cout<<"Chargement du maillage"<<std::endl;
//     }
// 
//     geometry geom,vol;
//     load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());
//     load_node_gmsh(vol ,("disc_"+NbrToStr(lc)).c_str());
// 
//     mesh_1D Omega(geom);
//     mesh_2D Vol(vol);
// 
//     load_elt_gmsh(Omega,0);
//     load_elt_gmsh(Vol,0);
// 
//     gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
//     gmsh_clean(("disc_"+NbrToStr(lc)).c_str());
// 
// 
//     ////================== Calcul de la normale =====================////
// 
// 
//     nrml_1D n_(Omega);
// 
// 
//     ////============ Matrice pour le champs rayonne =================////
// 
//     const vect<R3>& node = get_node(vol);
//     int nbpt  = size(node);
// 
//     int nbelt = nb_elt(Omega);
//     P1_1D dof(Omega);
//     int nbdof = nb_dof(dof);
// 
//     Real kappa=10;
// 
//     gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
//     potential<P1_1D,SLP_2D> SLPop(kappa,n_);
//     potential<P1_1D,DLP_2D> DLPop(kappa,n_);
// 
//     progress bar("assembly", nbpt*nbelt,verbose);
//     for (int j=0; j<nbpt ;j++){
//         const N1& jj = j;
//         for (int k=0;k<nbelt;k++, bar++){
//             const elt_1D& tk = Omega[k];
//             const N2&     kk = dof[k];
// 
//             SLP (jj,kk) += SLPop(node[j],tk) ;
// 			DLP (jj,kk) += DLPop(node[j],tk) ;
// 		}
// 		
//     }
//     bar.end();
// 	
// 	for (int l=0;l<harmonics.size();l++){
// 		Real p = harmonics[l];
// 
// 		////======================= Trace solution ======================////
// 
// 		vect<Cplx> TraceDirichlet;resize(TraceDirichlet,nbdof);fill(TraceDirichlet,0.);
// 		vect<Cplx> TraceNeumann;  resize(TraceNeumann,nbdof);  fill(TraceNeumann,0.);
// 
// 		for (int j=0 ; j<nbdof ; j++){
// 			const std::vector<std::pair<int,int>>& elts = get_elts_of_dof(phi, j);
// 			const elt_1D& seg = Omega[j];
// 			
// 			R3 X=elts[0].;
// 
// 
// 			Real theta0 = atan (X0[1]/X0[0]);
// 			Real theta1 = atan (X1[1]/X1[0]);
// 
// 			if (X[0]<0 & X[1]>=0){
// 				theta0 += M_PI;
// 			}
// 			if (X[0]<0 & X[1]<0){
// 				theta0 -= M_PI;
// 			}
// 
// 
// 			C2 Vinc;
// 
// 			Vinc[0] = exp( iu*p*theta0 );
// 			Vinc[1] = exp( iu*p*theta1 );
// 
// 			TraceDirichlet [I] += 0.5*Vinc;
// 
// 
// 			Vinc[0] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta0 );
// 			Vinc[1] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta1 );
// 
// 			TraceNeumann[I] += 0.5*Vinc;
// 
// 		}
// 
// 
// 		////=================== Solution de référence ===================////
// 
// 		vect<Cplx> Ref;resize(Ref,nbpt);fill(Ref,0.);
// 		for (int j=0;j<nbpt;j++){
// 			R3 X; X[0] = node[j][0];X[1]= node[j][1]; X[2]= node[j][2];;
// 
// 			Real r     = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
// 			Real theta = atan (X[1]/X[0]);
// 
// 			if (X[0]<0 & X[1]>=0){
// 				theta += M_PI;
// 			}
// 			if (X[0]<0 & X[1]<0){
// 				theta -= M_PI;
// 			}
// 
// 			Ref[j]= boost::math::cyl_bessel_j(p,kappa*r)/boost::math::cyl_bessel_j(p,kappa*R)*exp( iu*p*theta );
// 		}
// 
// 
// 		////================ Construction de la solution ================////
// 
// 		vect<Cplx> S;resize(S,nbpt);fill(S,0.);
// 		vect<Cplx> D;resize(D,nbpt);fill(D,0.);
// 		vect<Real> Out;resize(Out,nbpt);fill(Out,0.);
// 		vect<Real> Out_ref;resize(Out_ref,nbpt);fill(Out_ref,0.);
// 		
// 
// 		mv_prod(S,SLP,TraceNeumann);
// 		mv_prod(D,DLP,TraceDirichlet);
// 		for (int i=0;i<nbpt;i++){
// 			Out[i]= std::abs(S[i]+D[i]);
// 			Out_ref[i]= std::abs(Ref[i]);
// 		}
// 
// 		write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc));
// 		write_gmsh(Vol,Out_ref,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc)+"_ref");
// 
// 	}
// }

///==================================================================================////
///============================= potential_node_2D ==================================////
///==================================================================================////

void potential_node_2D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0){

    ////=======================  Mesh building  =====================////

    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc,verbose);
    gmsh_disc  (("disc_"+NbrToStr(lc)).c_str(),R*0.9,lc,verbose);


    ////=======================  Mesh loading  ======================////

    if (verbose>0){
        std::cout<<"Chargement du maillage"<<std::endl;
    }

    geometry geom,vol;
    load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());
    load_node_gmsh(vol ,("disc_"+NbrToStr(lc)).c_str());

    mesh_1D Omega(geom);
    mesh_2D Vol(vol);

    load_elt_gmsh(Omega,0);
    load_elt_gmsh(Vol,0);

//     gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
//     gmsh_clean(("disc_"+NbrToStr(lc)).c_str());


    ////================== Calcul de la normale =====================////


    nrml_1D n_(Omega);


    ////============ Matrice pour le champs rayonne =================////

    const vect<R3>& node = get_node(vol);
    int nbpt  = size(node);

    int nbelt = nb_elt(Omega);
    P1_1D dof(Omega);
    int nbdof = nb_dof(dof);

    Real kappa=10;

    gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
    potential<P1_1D,SLP_2D> SLPop(kappa,n_);
    potential<P1_1D,DLP_2D> DLPop(kappa,n_);

    progress bar("assembly", nbpt*nbelt,verbose);
    for (int j=0; j<nbpt ;j++){
        const N1& jj = j;
		
		for (int k=0;k<nbdof;k++, bar++){
			const N1& kk = k;
            SLP (jj,kk) += SLPop(node[j],k) ;
			DLP (jj,kk) += DLPop(node[j],k) ;
			
		}
		
		
    }
    bar.end();
	
	for (int l=0;l<harmonics.size();l++){
		Real p = harmonics[l];

		////======================= Trace solution ======================////

		vect<Cplx> TraceDirichlet;resize(TraceDirichlet,nbdof);fill(TraceDirichlet,0.);
		vect<Cplx> TraceNeumann;  resize(TraceNeumann,nbdof);  fill(TraceNeumann,0.);

		for (int j=0 ; j<nbelt ; j++){
			const elt_1D& seg = Omega[j];
			const N2&     I   = dof[j];

			R3 X0; X0[0] =  seg[0][0];X0[1]=seg[0][1]; X0[2]=0;
			R3 X1; X1[0] =  seg[1][0];X1[1]=seg[1][1]; X1[2]=0;

			Real theta0 = atan (X0[1]/X0[0]);
			Real theta1 = atan (X1[1]/X1[0]);

			if (X0[0]<0 & X0[1]>=0){
				theta0 += M_PI;
			}
			if (X0[0]<0 & X0[1]<0){
				theta0 -= M_PI;
			}
			if (X1[0]<0 & X1[1]>=0){
				theta1 += M_PI;
			}
			if (X1[0]<0 & X1[1]<0){
				theta1 -= M_PI;
			}

			C2 Vinc;

			Vinc[0] = exp( iu*p*theta0 );
			Vinc[1] = exp( iu*p*theta1 );

			TraceDirichlet [I] += 0.5*Vinc;


			Vinc[0] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta0 );
			Vinc[1] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta1 );

			TraceNeumann[I] += 0.5*Vinc;

		}


		////=================== Solution de référence ===================////

		vect<Cplx> Ref;resize(Ref,nbpt);fill(Ref,0.);
		for (int j=0;j<nbpt;j++){
			R3 X; X[0] = node[j][0];X[1]= node[j][1]; X[2]= node[j][2];;

			Real r     = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
			Real theta = atan (X[1]/X[0]);

			if (X[0]<0 & X[1]>=0){
				theta += M_PI;
			}
			if (X[0]<0 & X[1]<0){
				theta -= M_PI;
			}

			Ref[j]= boost::math::cyl_bessel_j(p,kappa*r)/boost::math::cyl_bessel_j(p,kappa*R)*exp( iu*p*theta );
		}


		////================ Construction de la solution ================////

		vect<Cplx> S;resize(S,nbpt);fill(S,0.);
		vect<Cplx> D;resize(D,nbpt);fill(D,0.);
		vect<Real> Out;resize(Out,nbpt);fill(Out,0.);
		vect<Real> Out_ref;resize(Out_ref,nbpt);fill(Out_ref,0.);
		

		mv_prod(S,SLP,TraceNeumann);
		mv_prod(D,DLP,TraceDirichlet);
		for (int i=0;i<nbpt;i++){
			Out[i]= std::abs(S[i]+D[i]);
			Out_ref[i]= std::abs(Ref[i]);
		}

		write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc));
		write_gmsh(Vol,Out_ref,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc)+"_ref");

	}
}

///==================================================================================////
///============================= potential_elt_3D ===================================////
///==================================================================================////
void potential_elt_3D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose=0){

    ////=======================  Mesh building  =====================////

    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_sphere(("sphere_"+NbrToStr(lc)).c_str(),R,lc,verbose);
    gmsh_ball  (("ball_"+NbrToStr(lc)).c_str(),R*0.9,lc,verbose);


    ////=======================  Mesh loading  ======================////

    if (verbose>0){
        std::cout<<"Chargement du maillage"<<std::endl;
    }

    geometry geom,vol;
    load_node_gmsh(geom,("sphere_"+NbrToStr(lc)).c_str());
    load_node_gmsh(vol ,("ball_"+NbrToStr(lc)).c_str());

    mesh_2D Omega(geom);
    mesh_3D Vol(vol);

    load_elt_gmsh(Omega,0);
    load_elt_gmsh(Vol,0);

    gmsh_clean(("sphere_"+NbrToStr(lc)).c_str());
    gmsh_clean(("ball_"+NbrToStr(lc)).c_str());

    ////================== Calcul de la normale =====================////


    nrml_2D n_(Omega);


    ////============ Matrice pour le champs rayonne =================////

    const vect<R3>& node = get_node(vol);
    int nbpt  = size(node);

    int nbelt = nb_elt(Omega);
    P1_2D dof(Omega);
    int nbdof = nb_dof(dof);

    Real kappa=1;

    gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
    potential<P1_2D,SLP_3D> SLPop(kappa,n_);
    potential<P1_2D,DLP_3D> DLPop(kappa,n_);

    progress bar("assembly", nbpt*nbelt,verbose);
    for (int j=0; j<nbpt ;j++){
        const N1& jj = j;
        for (int k=0;k<nbelt;k++, bar++){
            const elt_2D& tk = Omega[k];
            const N3&     kk = dof[k];

            SLP (jj,kk) += SLPop(node[j],tk) ;
			DLP (jj,kk) += DLPop(node[j],tk) ;
		}
		
    }
    bar.end();
	write(SLP,"SLP");
	for (int l=0;l<harmonics.size();l++){
		 Real p = harmonics[l];
		 
		////======================= Trace solution ======================////
		 
		vect<Cplx> TraceDirichlet;resize(TraceDirichlet,nbdof);fill(TraceDirichlet,0.);
		vect<Cplx> TraceNeumann;  resize(TraceNeumann,nbdof);  fill(TraceNeumann,0.);


		for (int j=0 ; j<nbelt ; j++){
			const elt_2D& seg = Omega[j];
			const N3&     I   = dof[j];
            
            R3 X0; X0[0] =  seg[0][0];X0[1]=seg[0][1]; X0[2]=seg[0][2];
            R3 X1; X1[0] =  seg[1][0];X1[1]=seg[1][1]; X1[2]=seg[1][2];
            R3 X2; X2[0] =  seg[2][0];X2[1]=seg[2][1]; X2[2]=seg[2][2];


            Real radius0 = sqrt(X0[0] * X0[0] + X0[1] * X0[1] + X0[2] * X0[2]);
			Real radius1 = sqrt(X1[0] * X1[0] + X1[1] * X1[1] + X1[2] * X1[2]);
			Real radius2 = sqrt(X0[0] * X0[0] + X0[1] * X0[1] + X0[2] * X2[2]);

            Real theta0 = acos(X0[2]/radius0);
            Real theta1 = acos(X1[2]/radius1);
            Real theta2 = acos(X2[2]/radius2);

            Real phi0 = atan2(X0[1],X0[0]);
            Real phi1 = atan2(X1[1],X1[0]);
            Real phi2 = atan2(X2[1],X2[0]);
            
//             if (X0[0]<0 & X0[1]>=0){
//                 phi0 += M_PI;
//             }
//             if (X0[0]<0 & X0[1]<0){
//                 phi0 -= M_PI;
//             }
//             if (X1[0]<0 & X1[1]>=0){
//                 phi1 += M_PI;
//             }
//             if (X1[0]<0 & X1[1]<0){
//                 phi1 -= M_PI;
//             }
//             if (X2[0]<0 & X2[1]>=0){
//                 phi2 += M_PI;
//             }
//             if (X2[0]<0 & X2[1]<0){
//                 phi2 -= M_PI;
//             }
//             std::cout << "X : "<<X0 <<" "<<X1<<" "<<X2<<std::endl;
//             std::cout << "theta : "<<theta0 <<" "<<theta1<<" "<<theta2<<std::endl;
// 			std::cout << "phi : "<< phi0 <<" "<<phi1<<" "<<phi2<<std::endl;
			
            C3 Vinc;
            
            Vinc[0] = boost::math::spherical_harmonic(p,p,theta0,phi0);
            Vinc[1] = boost::math::spherical_harmonic(p,p,theta1,phi1);
            Vinc[2] = boost::math::spherical_harmonic(p,p,theta2,phi2);

            TraceDirichlet [I] += 1./3. * Vinc;

            Cplx snn;
            snn = kappa*(((p+1.5)/(kappa*R)) - boost::math::sph_bessel(p+1,kappa*R)/boost::math::sph_bessel(p,kappa*R));

            Vinc[0] *= snn;
            Vinc[1] *= snn;
            Vinc[2] *= snn;

			TraceNeumann[I] += 1./3. * Vinc;
			

		}

		////=================== Solution de référence ===================////

		vect<Cplx> Ref;resize(Ref,nbpt);fill(Ref,0.);
		for (int j=0;j<nbpt;j++){
			R3 X; X[0] = node[j][0];X[1]= node[j][1]; X[2]= node[j][2];;

			Real r     = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
			
			Real theta;
			if (r<1e-10)
				theta=0; // to avoid nan value
			else
				theta = acos(X[2]/r);
            
			Real phi = atan2(X[1],X[0]);

			Ref[j]= boost::math::sph_bessel(p,kappa*r)/boost::math::sph_bessel(p,kappa*R)*boost::math::spherical_harmonic(p,p,theta,phi);
			if (j==0)
				std::cout << Ref[j]<<std::endl;
		}

		////================ Construction de la solution ================////

		vect<Cplx> S;resize(S,nbpt);fill(S,0.);
		vect<Cplx> D;resize(D,nbpt);fill(D,0.);
		vect<Real> Out;resize(Out,nbpt);fill(Out,0.);
		vect<Real> Out_ref;resize(Out_ref,nbpt);fill(Out_ref,0.);
		

		mv_prod(S,SLP,TraceNeumann);
		mv_prod(D,DLP,TraceDirichlet);
		for (int i=0;i<nbpt;i++){
// 			std::cout << S[i]<< std::endl; 
			Out[i]= std::abs(S[i]+D[i]);
			Out_ref[i]= std::abs(Ref[i]);
		}

		write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc));
		write_gmsh(Vol,Out_ref,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc)+"_ref");

	}
}
// class MyMatrix: public htool::VirtualMatrix{
// private:
// 	potential<P1_1D,SLP_2D>& Op;
// 	const geometry& geom;
// 	
// public:
// 	MyMatrix(potential<P1_1D,SLP_2D>& SLPop,const geometry& geome, const int& nbpts, const int& nbdof): Op(SLPop),geom(geome) {
// // 		potential<P1_1D,SLP_2D> Op=SLPop;
// 		nr = nbpts;
// 		nc = nbdof;
// 	}
// 
// 	const Cplx get_coef(const int& j, const int& k) const{
// // 		std::cout << get_node(geom,j) << std::endl;
// 		return Op(get_node(geom,j),k);
// 	}
// };


int main(int argc, char* argv[]){
	
	std::vector<Real> harmonics(3);harmonics[0]=1;harmonics[1]=2;harmonics[2]=3;
	Real lc = 0.07;
	Real R=0.5;
	
// 	potential_node_2D(harmonics,lc,  R, "potential_node_2D",1);
	potential_elt_3D(harmonics,lc,  R, "potential_elt_3D",1);
	
}
//     MPI_Init(&argc, &argv);
//     /*# Init #*/
//     int rankWorld, sizeWorld;
//     MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);
//     MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);
// 	htool::SetNdofPerElt(1);
// 	htool::SetEpsilon(0.9);
// 	htool::SetEta(0.9);
// 	
//     //// Rayon du cercle
//     Real R= 1.;
// // 
// //     //// Harmoniques
// //     std::vector<Real> harmoniques;harmoniques.push_back(1);harmoniques.push_back(2);harmoniques.push_back(3);harmoniques.push_back(4);harmoniques.push_back(5);
// // 
// // 	//// Finesse
// 	Real lc = 0.1;
// 	
// 	//// Tests
// // 	potential_elt_2D (harmoniques, finesse, R, "potential_elt_2D");
// // 	potential_node_2D (harmoniques, finesse, R, "potential_node_2D");
// 
// // 	potential_elt_3D (harmoniques, finesse, R, "potential_elt_3D");
// 	int verbose=100;
// 	if (rankWorld==0){
// ////=============================================================////
//     ////=======================  Mesh building  =====================////
//     ////=============================================================////
//     if (verbose>0){
//         std::cout<<"Construction du maillage"<<std::endl;
//     }
// 
// 	
//     gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc,verbose);
//     gmsh_disc  (("disc_"+NbrToStr(lc)).c_str(),R*0.9,lc,verbose);
// 	}
//     ////=============================================================////
//     ////=======================  Mesh loading  ======================////
//     ////=============================================================////
//     if (verbose>0){
//         std::cout<<"Chargement du maillage"<<std::endl;
//     }
// 	MPI_Barrier(MPI_COMM_WORLD);
//     geometry geom,vol;
//     load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());
//     load_node_gmsh(vol ,("disc_"+NbrToStr(lc)).c_str());
// 
//     mesh_1D Omega(geom);
//     mesh_2D Vol(vol);
// 
//     load_elt_gmsh(Omega,0);
//     load_elt_gmsh(Vol,0);
// 	
// 	if (rankWorld==0){
// 	
// //     gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
// //     gmsh_clean(("disc_"+NbrToStr(lc)).c_str());
// 	}
//     ////=============================================================////
//     ////================== Calcul de la normale =====================////
//     ////=============================================================////
// 
//     nrml_1D n_(Omega);
// 
//     ////=============================================================////
//     ////============ Matrice pour le champs rayonne =================////
//     ////=============================================================////
//     const vect<R3>& node = get_node(vol);
//     int nbpt  = size(node);
// 
// 	const vect<R3>& node_g = get_node(geom);
// //     int nbpt  = size(node);
// 
// 	
//     int nbelt = nb_elt(Omega);
//     P1_1D dof(Omega);
//     int nbdof = nb_dof(dof);
// 
//     Real kappa=10;
// 
//     gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
//     potential<P1_1D,SLP_2D> SLPop(kappa,n_);
//     potential<P1_1D,DLP_2D> DLPop(kappa,n_);
// 	MyMatrix vmat_SLP(SLPop,vol,nbpt,nbdof);
// // 	MyMatrix vmat_DLP(DLPop,nbpt,nbdof);
// 	
// 	
// 	std::vector<htool::Real> rt(nbpt);
// 	std::vector<htool::Real> rs(nbdof);
// 	std::vector<int> tabt(nbpt);
// 	std::vector<int> tabs(nbdof);
// 	
// 	std::vector<htool::R3> xt(nbpt);
// 	std::vector<htool::R3> xs(nbdof);
// 	
// 	for (int i=0;i<nbpt;i++){
// 		tabt[i]=i;
// 		htool::R3 pt;pt[0]=node[i][0];pt[1]=node[i][1];pt[2]=node[i][2];
// 		xt[i]=pt;
// 		rt[i]=0;
// 	}
// 	
// 	for (int i=0;i<nbdof;i++){
// 		tabs[i]=i;
// 		htool::R3 pt;pt[0]=node_g[i][0];pt[1]=node_g[i][1];pt[2]=node_g[i][2];
// 		xs[i]=pt;
// 		rs[i]=0;
// 	}
// 	
// 	
// 	htool::HMatrix hmat_SLP(vmat_SLP,xt,rt,tabt,xs,rs,tabs);
// // 	htool::HMatrix hmat_DLP(vmat_DLP,xt,rt,tabt,xs,rs,tabs);
// 	
// //     progress bar("assembly", nbpt*nbelt,verbose);
//     for (int j=0; j<nbpt ;j++){
//         const N1& jj = j;
//         for (int k=0;k<nbelt;k++){
//             const elt_1D& tk = Omega[k];
//             const N2&     kk = dof[k];
// 
//             SLP (jj,kk) += SLPop(node[j],tk) ;
// 			DLP (jj,kk) += DLPop(node[j],tk) ;
// 		}
// 		
//     }
// //     bar.end();
// 	
// // 	for (int l=0;l<harmonics.size();l++){
// 		Real p =1;
// 		////=============================================================////
// 		////======================= Trace solution ======================////
// 		////=============================================================////
// 		vect<Cplx> TraceDirichlet;resize(TraceDirichlet,nbdof);fill(TraceDirichlet,0.);
// 		vect<Cplx> TraceNeumann;  resize(TraceNeumann,nbdof);  fill(TraceNeumann,0.);
// 
// 		for (int j=0 ; j<nbelt ; j++){
// 			const elt_1D& seg = Omega[j];
// 			const N2&     I   = dof[j];
// 
// 			R3 X0; X0[0] =  seg[0][0];X0[1]=seg[0][1]; X0[2]=0;
// 			R3 X1; X1[0] =  seg[1][0];X1[1]=seg[1][1]; X1[2]=0;
// 
// 			Real theta0 = atan (X0[1]/X0[0]);
// 			Real theta1 = atan (X1[1]/X1[0]);
// 
// 			if (X0[0]<0 & X0[1]>=0){
// 				theta0 += M_PI;
// 			}
// 			if (X0[0]<0 & X0[1]<0){
// 				theta0 -= M_PI;
// 			}
// 			if (X1[0]<0 & X1[1]>=0){
// 				theta1 += M_PI;
// 			}
// 			if (X1[0]<0 & X1[1]<0){
// 				theta1 -= M_PI;
// 			}
// 
// 			C2 Vinc;
// 
// 			Vinc[0] = exp( iu*p*theta0 );
// 			Vinc[1] = exp( iu*p*theta1 );
// 
// 			TraceDirichlet [I] += 0.5*Vinc;
// 
// 
// 			Vinc[0] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta0 );
// 			Vinc[1] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta1 );
// 
// 			TraceNeumann[I] += 0.5*Vinc;
// 
// 		}
// 
// 		////=============================================================////
// 		////=================== Solution de référence ===================////
// 		////=============================================================////
// 		vect<Cplx> Ref;resize(Ref,nbpt);fill(Ref,0.);
// 		for (int j=0;j<nbpt;j++){
// 			R3 X; X[0] = node[j][0];X[1]= node[j][1]; X[2]= node[j][2];;
// 
// 			Real r     = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
// 			Real theta = atan (X[1]/X[0]);
// 
// 			if (X[0]<0 & X[1]>=0){
// 				theta += M_PI;
// 			}
// 			if (X[0]<0 & X[1]<0){
// 				theta -= M_PI;
// 			}
// 
// 			Ref[j]= boost::math::cyl_bessel_j(p,kappa*r)/boost::math::cyl_bessel_j(p,kappa*R)*exp( iu*p*theta );
// 		}
// 
// 		////=============================================================////
// 		////================ Construction de la solution ================////
// 		////=============================================================////
// 		vect<Cplx> S;resize(S,nbpt);fill(S,0.);
// 		vect<Cplx> Sh;resize(Sh,nbpt);fill(Sh,0.);
// 		vect<Cplx> D;resize(D,nbpt);fill(D,0.);
// 		vect<Real> Out;resize(Out,nbpt);fill(Out,0.);
// 		vect<Real> Out_ref;resize(Out_ref,nbpt);fill(Out_ref,0.);
// 		std::string output_name="test";
// 
// 		Real time=MPI_Wtime();
// 		mv_prod(S,SLP,TraceNeumann);
// 		time=MPI_Wtime()-time;
// 		
// 		mv_prod(D,DLP,TraceDirichlet);
// 		for (int i=0;i<nbpt;i++){
// 			Out[i]= std::abs(S[i]);
// 			Out_ref[i]= std::abs(Ref[i]);
// 		}
// 
// 		if (rankWorld==0){
// 		std::cout << "MvProd dense "<<time << std::endl;
// 		
// 		write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc));
// 		write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc)+"_ref");
// 		}
// 		
// 		MvProdMPI(Sh,hmat_SLP,TraceNeumann);
// // 		MvProdMPI(D,hmat_DLP,TraceDirichlet);
// 		for (int i=0;i<nbpt;i++){
// 			Out[i]= std::abs(Sh[i]);
// // 			Out_ref[i]= std::abs(Ref[i]);
// 		}
// 		
// 		Real err=0;
// 		Real norm=0;
// 		for (int i=0;i<nbpt;i++){
// 			err+=std::pow(std::abs(S[i]-Sh[i]),2);
// 			norm+=std::pow(std::abs(S[i]),2);
// 		}
// 		err=std::sqrt(err/norm);
// 		
// 		Real compression=CompressionRate(hmat_SLP);
// 		int nb_lr = nb_lrmats(hmat_SLP);
// 		int nb_dense = nb_densemats(hmat_SLP);
// 		Real err_frob = squared_absolute_error(hmat_SLP, vmat_SLP);
// 		if (rankWorld==0){
// 			std::cout << "nb_lrmats : "<<nb_lr<<std::endl;
// 			std::cout << "nb_densemats : "<<nb_dense<<std::endl;
// 			std::cout << "err_frob : "<<err_frob<<std::endl;
// 			std::cout << "err_l2 : "<<err<<std::endl;
// 			
// 			
// 		write_gmsh(Vol,Out,output_name+"h_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc));
// // 		write_gmsh(Vol,Out,output_name+"h_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc)+"_ref");
// 		
// 			std::cout << compression << std::endl;
// 		}
// 
// // 	}
// 		
// 		MPI_Finalize();
// }