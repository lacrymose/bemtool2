#include "tests_sphere_3D.hpp"
#include <math.h>
//#include <boost/math/special_functions/spheric_harmonic.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

///==================================================================================////
///==========================First kind Dirichlet====================================////
///==================================================================================////

// void first_kind_dirichlet_3D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose){
//     ////=============================================================////
//     ////=======================  Mesh building  =====================////
//     ////=============================================================////
//     if (verbose>0){
//         std::cout<<"Construction du maillage"<<std::endl;
//     }
//     gmsh_sphere(("sphere_"+NbrToStr(lc)).c_str(),R,lc);

   
    
//     ////=============================================================////
//     ////=======================  Mesh loading  ======================////
//     ////=============================================================////
//     if (verbose>0){
//         std::cout<<"Chargement du maillage"<<std::endl;
//     }
//     Real kappa=1.;
    
//     geometry geom;
//     load_node_gmsh(geom,("sphere_"+NbrToStr(lc)).c_str());
    
//     mesh_2D Omega(geom);
//     load_elt_gmsh(Omega,1);
//     gmsh_clean(("sphere_"+NbrToStr(lc)).c_str());




//     ////=============================================================////
//     ////================== Calcul de la normale =====================////
//     ////=============================================================////
//     if (verbose>0){
//         std::cout<<"Calcul de la normale"<<std::endl;
//     }

//     nrml_2D n_(Omega);

//     // 	for (int j=0;j < nb_elt(Omega);j++){
//     // 		elt_1D seg = Omega[j];
//     // 		R3 G;
//     // 		G[0]= 0.5 * ( seg[0][0] + seg[1][0] );
//     // 		G[1]= 0.5 * ( seg[0][1] + seg[1][1] );
//     // 		G[2]= 0;
//     //
//     // 		std::cout<<(G,n_[j])<<std::endl;
//     // 	}
    
//     // swap(n_);
//     ////=============================================================////
//     ////================ Assemblage de la matrice ===================////
//     ////=============================================================////
//     if (verbose>0){
//         std::cout<<"Assemblage operateurs integraux"<<std::endl;
//     }

//     int nbelt = nb_elt(Omega);
//     P1_2D dof(Omega);
//     int nbdof = nb_dof(dof);

//     if (verbose>0){
//         std::cout<<"Matrice de taille: ";
//         std::cout<<nbdof<< std::endl;
//         std::cout<<"Nel: ";
//         std::cout<<nbelt<< std::endl;
//     }
    
    
//     gmm_dense V(nbdof,nbdof),K(nbdof,nbdof),M(nbdof,nbdof),CS(nbdof,nbdof);
//     bem<P1_2D,P1_2D, SLP_3D>   Vop(kappa,n_,n_);
//     bem<P1_2D,P1_2D, DLP_3D>   Kop(kappa,n_,n_);

//     bem<P1_2D,P1_2D, CST_3D>   Cop(kappa,n_,n_);
    
//     progress bar("assembly", nbelt*nbelt, verbose);
//     for(int j=0; j<nbelt; j++){
//         const elt_2D& tj = Omega[j];
//         const N3&     jj = dof[j];
        
//         for(int k=0; k<nbelt; k++,bar++){
//             const elt_2D& tk = Omega[k];
//             const N3&     kk = dof[k];

//             V(jj,kk) += Vop (tj,tk);
//             //CS(jj,kk) += Cop (tj,tk);
//             //std::cout<<Cop (tj,tk)<<std::endl;
//             K(jj,kk) += Kop (tj,tk);
//             //std::cout<< tj;
//             //std::cout<< tk;
//             //std::cout<< Vop (tj,tk) << std::endl << std::endl;

//         }

        
//         M(jj,jj) += MassP1(tj);
//         //std::cout<< MassP1(tj) << std::endl << std::endl;
//     }
//     bar.end();




    
    
//     for (int l=0;l<harmonics.size();l++){
//         Real p = harmonics[l];
        
//         ////=============================================================////
//         ////================== Harmonique de Fourier ====================////
//         ////=============================================================////
        
//         vect<Cplx> Ep ;resize(Ep ,nbdof);fill(Ep ,0.);
// 		vect<Cplx> Ref;resize(Ref,nbdof);fill(Ref,0.);
//         for (int j=0 ; j<nbelt ; j++){
//             const elt_2D& seg = Omega[j];
//             const N3&     I   = dof[j];
            
//             R3 X0; X0[0] =  seg[0][0];X0[1]=seg[0][1]; X0[2]=seg[0][2];
//             R3 X1; X1[0] =  seg[1][0];X1[1]=seg[1][1]; X1[2]=seg[1][2];
//             R3 X2; X2[0] =  seg[2][0];X2[1]=seg[2][1]; X2[2]=seg[2][2];
            
//             //aa

//             //Real theta0 = atan (X0[1]/X0[0]);
//             //Real theta1 = atan (X1[1]/X1[0]);


//             // OK POUR UN RAYON DE 1 ATTENTION

//             Real theta0 = acos(X0[2]);
//             Real theta1 = acos(X1[2]);
//             Real theta2 = acos(X2[2]);

//             Real phi0 = atan(X0[1]/X0[0]);
//             Real phi1 = atan(X1[1]/X1[0]);
//             Real phi2 = atan(X2[1]/X2[0]);


//             //Real phi0 = atan2(X0[1],X0[0]);
//             //Real phi1 = atan2(X1[1],X1[0]);
//             //Real phi2 = atan2(X2[1],X2[0]);

//             //std::cout << "---------------"<< std::endl;
//             //std::cout << X0 << std::endl;
//             ////std::cout << X1 << std::endl;
//             ////std::cout << X2 << std::endl;
//             //std::cout << theta0 << "  " << phi0 << std::endl;

//             //std::cout << "---------------"<< std::endl;
                

//             std::cout << phi0 << "phi0 avant ch" <<std::endl;
            
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

//             std::cout << phi0 << "phi0 apres ch" <<std::endl;
            
            
//             C3 Vinc;
            
//             //Vinc[0] = exp( iu*p*theta0 );
//             //Vinc[1] = exp( iu*p*theta1 );

//             //std::cout << p <<"theta0:"<< theta0 <<"phi0" << phi0 << std::endl;
//             Vinc[0] = boost::math::spherical_harmonic(p,p,theta0,phi0);
//             Vinc[1] = boost::math::spherical_harmonic(p,p,theta1,phi1);
//             Vinc[2] = boost::math::spherical_harmonic(p,p,theta2,phi2);
            
            
//             //Vinc[0] = exp( iu*kappa*X0[0]);
//             //Vinc[1] = exp( iu*kappa*X1[0]);
//             //Vinc[2] = exp( iu*kappa*X2[0]);


//             Ep [I] += 1./3 * Vinc;

            
			
// 			//Vinc[0] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta0 );
//             //Vinc[1] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta1 );
//             Cplx snn;
//             snn = (iu * kappa * boost::math::sph_bessel(p,kappa*R) * boost::math::sph_hankel_1(p,kappa*R));


//             Vinc[0] *= snn;
//             Vinc[1] *= snn;
//             Vinc[2] *= snn;

// 			Ref[I] += 1./3 * Vinc;
//         }

        
//         ////=============================================================////
//         ////================ Assemblage du second membre ================////
//         ////=============================================================////
    
//         vect<Cplx> F; resize(F,nbdof); fill(F,0);
//         vect<Cplx> Ftemp; resize(Ftemp,nbdof); fill(Ftemp,0);
//         vect<Cplx> gD; resize(gD,nbdof); fill(gD,0);
    
//         // Boundary condition
//         gD=Ep;
//         //fill(gD,Cplx(0.));
    
//         mv_prod(F,K,gD);
//         mv_prod(Ftemp,M,gD);
    
//         for(int j=0; j<nbelt; j++){
//             F[j] = 0.5*Ftemp[j] - F[j];
            
//         }

//         ////=============================================================////
//         ////================ Résolution système linéaire ================////
//         ////=============================================================////
//         if (verbose>0){
//         std::cout<<"Appel du solveur"<<std::endl;
//         }
//         vect<Cplx> U;
//         resize(U,nbdof);
    
//         // 	gmm_dense LU(nbddl,nbddl);
//         // 		lu_factor(J,LU);
//         // 		lu_solve(LU,U,F);
//         gmres_solve(V,U,F,40,verbose);

	
//         ////=============================================================////
//         ////===================== Calcul de l'erreur ====================////
//         ////=============================================================////
//         vect<Cplx> Err, Err2, Norme;
//         resize(Err, nbdof);
// 		resize(Err2,nbdof);
// 		resize(Norme,nbdof);
// 		for(int j=0; j<nbdof; j++){
// 			Err[j] =  U[j]-Ref[j];
// 		}
// 		mv_prod(Err2,M,Err);
// 		mv_prod(Norme,M,Ref);
		
// 		Cplx val=0;
// 		Cplx norme =0.;
// 		Real erreur=0;
// 		for(int j=0; j<nbdof; j++){
// 			val += Err2[j]*conj(Err[j]);
// 			norme  += Norme[j]*conj(Ref[j]);
// 		}
// 		erreur=abs(val/norme);
// 		if (verbose>0){
// 			std::cout << "erreur:\t" << erreur << std::endl;
// 		}
// 		////=============================================================////
// 		////======================== Sauvegardes ========================////
// 		////=============================================================////
// 		std::ofstream output((output_name+"_"+NbrToStr<Real>(p)+".txt").c_str(),std::ios::app);
// 		if (verbose>0){
// 			std::cout<<"Output in "<<output_name<<std::endl;
// 		}
// 		if (!output){
// 			std::cerr<<"Output file cannot be created"<<std::endl;
// 			exit(1);
// 		}
// 		else{
// 			output<<lc<<" "<<erreur<<std::endl;
// 		}
// 		output.close();
// 	}

  
// }



///==================================================================================////
///============================Fourier harmonics=====================================////
///==================================================================================////

void fourier_harmonic_3D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name,int verbose){
    
    ////=============================================================////
    ////=======================  Mesh building  =====================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_sphere(("sphere_"+NbrToStr(lc)).c_str(),R,lc,verbose);
    
    ////=============================================================////
    ////=======================  Mesh loading  ======================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Chargement du maillage"<<std::endl;
    }
    Real kappa=1.;
    geometry geom;
    load_node_gmsh(geom,("sphere_"+NbrToStr(lc)).c_str());
    
    mesh_2D Omega(geom);
    load_elt_gmsh(Omega,0);
    //gmsh_clean(("sphere_"+NbrToStr(lc)).c_str());
    ////=============================================================////
    ////================== Calcul de la normale =====================////
    ////=============================================================////
    
    nrml_2D n_(Omega);
    
    ////=============================================================////
    ////================ Assemblage de la matrice ===================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Assemblage operateurs integraux"<<std::endl;
    }
    int nbelt = nb_elt(Omega);
    P1_2D dof(Omega);
    int nbdof = nb_dof(dof);

    if (verbose>0){
        std::cout<<"Matrice de taille: ";
        std::cout<<nbdof<< std::endl;
        std::cout<<"Nel: ";
        std::cout<<nbelt<< std::endl;
    }

    gmm_dense V(nbdof,nbdof),W(nbdof,nbdof),K(nbdof,nbdof),M(nbdof,nbdof);
    bem<P1_2D,P1_2D, SLP_3D>   Vop (kappa,n_,n_);
    bem<P1_2D,P1_2D, HSP_3D>   Wop (kappa,n_,n_);
    bem<P1_2D,P1_2D, DLP_3D>   Kop (kappa,n_,n_);
    
    progress bar("assembly", nbelt*nbelt,verbose);
    for(int j=0; j<nbelt; j++){
        const elt_2D& tj = Omega[j];
        const N3&     jj = dof[j];
        
        for(int k=0; k<nbelt; k++,bar++){
            const elt_2D& tk = Omega[k];
            const N3&     kk = dof[k];
            V (jj,kk) += Vop  (tj,tk);
            W (jj,kk) += Wop  (tj,tk);
            K (jj,kk) += Kop  (tj,tk);
    
        }
        
        M(jj,jj) += MassP1(tj);
    }
    bar.end();
    

    // for (int l=0;l<harmonics.size();l++){
    //     Real p = harmonics[l];
        
    //     ////=============================================================////
    //     ////================== Harmonique de Fourier ====================////
    //     ////=============================================================////
        
    //     vect<Cplx> Ep;resize(Ep,nbdof);fill(Ep,0.);
        
    //     for (int j=0 ; j<nbelt ; j++){
    //         const elt_2D& seg = Omega[j];
    //         const N3&     I   = dof[j];
                
    //         R3 X0; X0[0] =  seg[0][0];X0[1]=seg[0][1]; X0[2]=seg[0][2];
    //         R3 X1; X1[0] =  seg[1][0];X1[1]=seg[1][1]; X1[2]=seg[1][2];
    //         R3 X2; X2[0] =  seg[2][0];X2[1]=seg[2][1]; X2[2]=seg[2][2];

    //         Real phi0 = atan2(X0[1],X0[0]);
    //         Real phi1 = atan2(X1[1],X1[0]);
    //         Real phi2 = atan2(X2[1],X2[0]);

    //         Real theta0 = acos(X0[2]);
    //         Real theta1 = acos(X1[2]);
    //         Real theta2 = acos(X2[2]);

    //         if (phi0<0){
    //             phi0 += 2 * M_PI;
    //         }
    //         if (phi1<0){
    //             phi1 += 2 * M_PI;
    //         }
    //         if (phi2<0){
    //             phi2 += 2 * M_PI;
    //         }

    //         C3 Vinc;

    //         Vinc[0] =  boost::math::spherical_harmonic(p,p,theta0,phi0);
    //         Vinc[1] =  boost::math::spherical_harmonic(p,p,theta1,phi1);
    //         Vinc[2] =  boost::math::spherical_harmonic(p,p,theta2,phi2);

    //         Ep[I] = Vinc;   
    //     }

    //     vect<Cplx> Eph;resize(Eph,nbdof);fill(Eph,0.);
    //     mv_prod(Eph,M,Ep);


	
	
	
	
    // 	////=============================================================////
    //     ////================== Calcul valeurs propres ===================////
    //     ////=============================================================////
        
        
    //     ////===================================////
    //     ////===== Test avec orthogonalité =====////
        
    //     Cplx val =0;
    //     Cplx err =0;
    //     for(int j=0; j<nbdof; j++){
    //     	val += Ep[j]*conj(Eph[j]);
    //     }
    
    //     err = abs(val)-1;
    //     std::cout<<"PS : "<<err<<std::endl;

        
    //     ////===================================////
    //     ////=========== Test avec V ===========////
        
        
    //     val =0;
    //     Real err_V =0;
        
    //     Cplx j_p = boost::math::sph_bessel(p, kappa*R);
    //     Cplx y_p = boost::math::sph_neumann(p, kappa*R);
    //     Cplx h_p = boost::math::sph_hankel_1(p, kappa*R);
        
    //     Cplx j_p_prime = boost::math::sph_bessel_prime(p, kappa*R);
    //     Cplx y_p_prime = boost::math::sph_neumann_prime(p, kappa*R);
        
    //     Cplx h_p_prime = j_p_prime + iu * y_p_prime;
        
    //     Cplx ref = iu * kappa * j_p * h_p;


    //    vect<Cplx> Temp;resize(Temp,nbdof);fill(Temp,0.);
    //    mv_prod(Temp,V,Ep);
       
    //    for(int j=0; j<nbdof; j++){
    //        val    += conj(Ep[j])*Temp[j];
    //    }
    //    //
    //    err_V = abs(val-ref) /abs(ref);
    //    if (verbose>0){
    //        std::cout<<"V : "<<err_V<<std::endl;
    //    }

        
        
    //     ////===================================////
    //     ////=========== Test avec W ===========////
        
    //     val =0;
    //     Real err_W =0;
    //     ref = - iu * pow(kappa,3) * j_p_prime * h_p_prime;
        
    //     fill(Temp,0.);
    //     mv_prod(Temp,W,Ep);
    //     for(int j=0; j<nbdof; j++){
    //         val    += conj(Ep[j])*Temp[j];
    //     }
        
    //     err_W = abs(val-ref) /abs(ref);
    //     if (verbose>0){
    //         std::cout<<"W : "<<err_W<<std::endl;
    //     }
    //     ////===================================////
    //     ////=========== Test avec K ===========////
        
    //     val =0;
    //     Real err_K =0;
    //     ref = iu /2. * kappa*kappa * (j_p_prime * h_p + j_p * h_p_prime);


        
    //     fill(Temp,0.);
    //     mv_prod(Temp,K,Ep);

    //     for(int j=0; j<nbdof; j++){
    //         val    += conj(Ep[j])*Temp[j];
    //     }
    //     err_K = abs(val-ref) /abs(ref);
    //     if (verbose>0){
    //         std::cout<<"K : "<<err_K<<std::endl;
    //     }
        
    //     ////=============================================================////
    //     ////======================== Sauvegardes ========================////
    //     ////=============================================================////
    //     std::ofstream output((output_name+"_"+NbrToStr<Real>(p)+".txt").c_str(),std::ios::app);
    //     if (!output){
    //         std::cerr<<"Output file cannot be created"<<std::endl;
    //         exit(1);
    //     }
    //     else{
    //         output<<lc<<" "<<err_V<<" "<<err_W<<" "<<err_K<<std::endl;
    //     }

    //     output.close();
        
    
    // if (l == 0)
    // {
    //     std::cout << p << "P" << std::endl;
    // }
    // }




    
}
