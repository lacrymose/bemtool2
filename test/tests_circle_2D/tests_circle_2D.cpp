#include "tests_circle_2D.hpp"
#include <math.h>



///==================================================================================////
///==========================First kind Dirichlet====================================////
///==================================================================================////

void first_kind_dirichlet_2D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose){
    ////=============================================================////
    ////=======================  Mesh building  =====================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc);
    
    
    ////=============================================================////
    ////=======================  Mesh loading  ======================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Chargement du maillage"<<std::endl;
    }
    Real kappa=1.;
    
    geometry geom;
    load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());
    
    mesh_1D Omega(geom);
    load_elt_gmsh(Omega,0);
    gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
    ////=============================================================////
    ////================== Calcul de la normale =====================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Calcul de la normale"<<std::endl;
    }
    nrml_1D n_(Omega);
    
    // 	for (int j=0;j < nb_elt(Omega);j++){
    // 		elt_1D seg = Omega[j];
    // 		R3 G;
    // 		G[0]= 0.5 * ( seg[0][0] + seg[1][0] );
    // 		G[1]= 0.5 * ( seg[0][1] + seg[1][1] );
    // 		G[2]= 0;
    //
    // 		std::cout<<(G,n_[j])<<std::endl;
    // 	}
    
    // swap(n_);
    ////=============================================================////
    ////================ Assemblage de la matrice ===================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Assemblage operateurs integraux"<<std::endl;
    }
    int nbelt = nb_elt(Omega);
    P1_1D dof(Omega);
    int nbdof = nb_dof(dof);
    
    
    gmm_dense V(nbdof,nbdof),K(nbdof,nbdof),M(nbdof,nbdof);
    bem<P1_1D,P1_1D, SLP_2D>   Vop(kappa,n_,n_);
    bem<P1_1D,P1_1D, DLP_2D>   Kop(kappa,n_,n_);
    
    progress bar("assembly", nbelt*nbelt, verbose);
    for(int j=0; j<nbelt; j++){
        const elt_1D& tj = Omega[j];
        const N2&     jj = dof[j];
        
        for(int k=0; k<nbelt; k++,bar++){
            const elt_1D& tk = Omega[k];
            const N2&     kk = dof[k];
            
            V(jj,kk) += Vop (tj,tk);
            K(jj,kk) += Kop (tj,tk);
            
        }
        
        M(jj,jj) += MassP1(tj);
    }
    bar.end();
    
    for (int l=0;l<harmonics.size();l++){
        Real p = harmonics[l];
        
        ////=============================================================////
        ////================== Harmonique de Fourier ====================////
        ////=============================================================////
        
        vect<Cplx> Ep ;resize(Ep ,nbdof);fill(Ep ,0.);
		vect<Cplx> Ref;resize(Ref,nbdof);fill(Ref,0.);
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
            
            Ep [I] += 0.5*Vinc;
			
			Vinc[0] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta0 );
            Vinc[1] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta1 );
            
			Ref[I] += 0.5*Vinc;
            
        }
        
        ////=============================================================////
        ////================ Assemblage du second membre ================////
        ////=============================================================////
    
        vect<Cplx> F; resize(F,nbdof); fill(F,0);
        vect<Cplx> Ftemp; resize(Ftemp,nbdof); fill(Ftemp,0);
        vect<Cplx> gD; resize(gD,nbdof); fill(gD,0);
    
        // Boundary condition
        gD=Ep;
        //fill(gD,Cplx(0.));
    
        mv_prod(F,K,gD);
        mv_prod(Ftemp,M,gD);
    
        for(int j=0; j<nbelt; j++){
            F[j] = 0.5*Ftemp[j] - F[j];
        }

        ////=============================================================////
        ////================ Résolution système linéaire ================////
        ////=============================================================////
        if (verbose>0){
        std::cout<<"Appel du solveur"<<std::endl;
        }
        vect<Cplx> U;
        resize(U,nbdof);
    
        // 	gmm_dense LU(nbddl,nbddl);
        // 		lu_factor(J,LU);
        // 		lu_solve(LU,U,F);
        gmres_solve(V,U,F,40,verbose);
    
    
        ////=============================================================////
        ////===================== Calcul de l'erreur ====================////
        ////=============================================================////
        vect<Cplx> Err, Err2, Norme;
        resize(Err, nbdof);
		resize(Err2,nbdof);
		resize(Norme,nbdof);
		for(int j=0; j<nbdof; j++){
			Err[j] =  U[j]-Ref[j];
		}
		mv_prod(Err2,M,Err);
		mv_prod(Norme,M,Ref);
		
		Cplx val=0;
		Cplx norme =0.;
		Real erreur=0;
		for(int j=0; j<nbdof; j++){
			val += Err2[j]*conj(Err[j]);
			norme  += Norme[j]*conj(Ref[j]);
		}
		erreur=abs(val/norme);
		if (verbose>0){
			std::cout << "erreur:\t" << erreur << std::endl;
		}
		////=============================================================////
		////======================== Sauvegardes ========================////
		////=============================================================////
		std::ofstream output((output_name+"_"+NbrToStr<Real>(p)+".txt").c_str(),std::ios::app);
		if (verbose>0){
			std::cout<<"Output in "<<output_name<<std::endl;
		}
		if (!output){
			std::cerr<<"Output file cannot be created"<<std::endl;
			exit(1);
		}
		else{
			output<<lc<<" "<<erreur<<std::endl;
		}
		output.close();
	}
}


///==================================================================================////
///===========================First kind Neumann=====================================////
///==================================================================================////

void first_kind_neumann_2D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose){
    
    ////=============================================================////
    ////=======================  Mesh building  =====================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc,verbose);
    
    ////=============================================================////
    ////=======================  Mesh loading  ======================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Chargement du maillage"<<std::endl;
    }
    Real kappa=1.;
    
    geometry geom;
    load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());
    
    mesh_1D Omega(geom);
    load_elt_gmsh(Omega,0);
    gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
    
    ////=============================================================////
    ////================== Calcul de la normale =====================////
    ////=============================================================////
    
    nrml_1D n_(Omega);
    
    // 	for (int j=0;j < nb_elt(Omega);j++){
    // 		elt_1D seg = Omega[j];
    // 		R3 G;
    // 		G[0]= 0.5 * ( seg[0][0] + seg[1][0] );
    // 		G[1]= 0.5 * ( seg[0][1] + seg[1][1] );
    // 		G[2]= 0;
    
    // 		std::cout<<(G,n_[j])<<std::endl;
    // 	}
    
//     swap(n_);
    ////=============================================================////
    ////================ Assemblage de la matrice ===================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Assemblage operateurs integraux"<<std::endl;
    }
    int nbelt = nb_elt(Omega);
    P1_1D dof(Omega);
    int nbdof = nb_dof(dof);
    
    
    gmm_dense W(nbdof,nbdof),TK(nbdof,nbdof),M(nbdof,nbdof);
    bem<P1_1D,P1_1D, HSP_2D>    Wop (kappa,n_,n_);
    bem<P1_1D,P1_1D, TDLP_2D>   TKop(kappa,n_,n_);
    
    progress bar("assembly", nbelt*nbelt,verbose);
    for(int j=0; j<nbelt; j++){
        const elt_1D& tj = Omega[j];
        const N2&     jj = dof[j];
        
        for(int k=0; k<nbelt; k++,bar++){
            const elt_1D& tk = Omega[k];
            const N2&     kk = dof[k];
            
            W (jj,kk) += Wop  (tj,tk);
            TK(jj,kk) += TKop (tj,tk);
            
        }
        
        M(jj,jj) += MassP1(tj);
    }
    bar.end();
    
	for (int l=0;l<harmonics.size();l++){
		Real p = harmonics[l];
	
	    ////=============================================================////
        ////================== Harmonique de Fourier ====================////
        ////=============================================================////
        
        vect<Cplx> Ep ;resize(Ep ,nbdof);fill(Ep ,0.);
		vect<Cplx> Ref;resize(Ref,nbdof);fill(Ref,0.);
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
            
            Ep [I] += 0.5*Vinc;
			
			Vinc[0] = (1./kappa)*boost::math::cyl_bessel_j(p,kappa*R)/((p/(kappa*R))*boost::math::cyl_bessel_j(p,kappa*R)-boost::math::cyl_bessel_j(p+1,kappa*R))*exp( iu*p*theta0 );
            Vinc[1] = (1./kappa)*boost::math::cyl_bessel_j(p,kappa*R)/((p/(kappa*R))*boost::math::cyl_bessel_j(p,kappa*R)-boost::math::cyl_bessel_j(p+1,kappa*R))*exp( iu*p*theta1 );
            
			Ref[I] += 0.5*Vinc;
            
        }
        
		////=============================================================////
		////================ Assemblage du second membre ================////
		////=============================================================////
		
		vect<Cplx> F; resize(F,nbdof); fill(F,0);
		vect<Cplx> Ftemp; resize(Ftemp,nbdof); fill(Ftemp,0);
		vect<Cplx> gN; resize(gN,nbdof); fill(gN,0);
		
		// Boundary condition
		gN=Ep;
		// 	fill(gD,Cplx(0.));
		
		mv_prod(F,TK,gN);
		mv_prod(Ftemp,M,gN);
		
		for(int j=0; j<nbelt; j++){
			F[j] = 0.5*Ftemp[j] - F[j];
		}
		

		////=============================================================////
		////================ Résolution système linéaire ================////
		////=============================================================////
		if (verbose>0){
			std::cout<<"Appel du solveur"<<std::endl;
		}
		vect<Cplx> U;
		resize(U,nbdof);
		
		// 	gmm_dense LU(nbddl,nbddl);
		// 		lu_factor(J,LU);
		// 		lu_solve(LU,U,F);
		gmres_solve(W,U,F,40,verbose);
		
		
		////=============================================================////
		////===================== Calcul de l'erreur ====================////
		////=============================================================////
		vect<Cplx> Err, Err2, Norme;
		resize(Err, nbdof);
		resize(Err2,nbdof);
		resize(Norme,nbdof);
		for(int j=0; j<nbdof; j++){
			Err[j] =  U[j]-Ref[j];
		}
		mv_prod(Err2,M,Err);
		mv_prod(Norme,M,Ref);
		
		Cplx val=0;
		Cplx norme =0.;
		Real erreur =0;
		for(int j=0; j<nbdof; j++){
			val += Err2[j]*conj(Err[j]);
			norme  += Norme[j]*conj(Ref[j]);
		}
		erreur=abs(val/norme);
		if (verbose>0){
			std::cout << "erreur:\t" << erreur << std::endl;
		}
		////=============================================================////
		////======================== Sauvegardes ========================////
		////=============================================================////
		std::ofstream output((output_name+"_"+NbrToStr<Real>(p)+".txt").c_str(),std::ios::app);
		if (verbose>0){
			std::cout<<"Output in "<<output_name<<std::endl;
		}
		if (!output){
			std::cerr<<"Output file cannot be created"<<std::endl;
			exit(1);
		}
		else{
			output<<lc<<" "<<erreur<<std::endl;
		}
		output.close();
	}
}

///==================================================================================////
///==========================Second kind Dirichlet===================================////
///==================================================================================////

void second_kind_dirichlet_2D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose){
    ////=============================================================////
    ////=======================  Mesh building  =====================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc,verbose);
    
    
    ////=============================================================////
    ////=======================  Mesh loading  ======================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Chargement du maillage"<<std::endl;
    }
    Real kappa=1.;
    
    geometry geom;
    load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());
    
    mesh_1D Omega(geom);
    load_elt_gmsh(Omega,0);
    gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
    
    ////=============================================================////
    ////================== Calcul de la normale =====================////
    ////=============================================================////
    
    nrml_1D n_(Omega);
    
    // 	for (int j=0;j < nb_elt(Omega);j++){
    // 		elt_1D seg = Omega[j];
    // 		R3 G;
    // 		G[0]= 0.5 * ( seg[0][0] + seg[1][0] );
    // 		G[1]= 0.5 * ( seg[0][1] + seg[1][1] );
    // 		G[2]= 0;
    
    // 		cout<<(G,n_[j])<<endl;
    // 	}
    
//     swap(n_);
    ////=============================================================////
    ////================ Assemblage de la matrice ===================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Assemblage operateurs integraux"<<std::endl;
    }
    int nbelt = nb_elt(Omega);
    P1_1D dof(Omega);
    int nbdof = nb_dof(dof);
    
    
    gmm_dense W(nbdof,nbdof),TT(nbdof,nbdof),M(nbdof,nbdof);
    bem<P1_1D,P1_1D, HSP_2D>    Wop (kappa,n_,n_);
    bem<P1_1D,P1_1D, TDLP_2D>   TKop(kappa,n_,n_);
    
    progress bar("assembly", nbelt*nbelt,verbose);
    for(int j=0; j<nbelt; j++){
        const elt_1D& tj = Omega[j];
        const N2&     jj = dof[j];
        
        for(int k=0; k<nbelt; k++,bar++){
            const elt_1D& tk = Omega[k];
            const N2&     kk = dof[k];
            
            W (jj,kk) += Wop  (tj,tk);
            TT(jj,kk) += TKop (tj,tk);
            
        }
        
        M (jj,jj) +=      MassP1(tj);
        TT(jj,jj) += -0.5*MassP1(tj);
    }
    bar.end();
    for (int l=0;l<harmonics.size();l++){
		Real p = harmonics[l];
		
        ////=============================================================////
        ////================== Harmonique de Fourier ====================////
        ////=============================================================////
        
        vect<Cplx> Ep ;resize(Ep ,nbdof);fill(Ep ,0.);
		vect<Cplx> Ref;resize(Ref,nbdof);fill(Ref,0.);
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
            
            Ep [I] += 0.5*Vinc;
			
			Vinc[0] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta0 );
            Vinc[1] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta1 );
            
			Ref[I] += 0.5*Vinc;
            
        }
        
		////=============================================================////
		////================ Assemblage du second membre ================////
		////=============================================================////
		vect<Cplx> F; resize(F,nbdof); fill(F,0);
		vect<Cplx> Ftemp; resize(Ftemp,nbdof); fill(Ftemp,0);
		vect<Cplx> gD; resize(gD,nbdof); fill(gD,0);
		
		// Boundary condition
		gD= Ep;
		// 	fill(gD,Cplx(0.));
		
		
		mv_prod(Ftemp,W,gD);
		for(int j=0; j<nbelt; j++){
			F[j] += -Ftemp[j];
		}
		
		////=============================================================////
		////================ Résolution système linéaire ================////
		////=============================================================////
		if (verbose>0){
			std::cout<<"Appel du solveur"<<std::endl;
		}
		
		vect<Cplx> U;
		resize(U,nbdof);
		
		// 	gmm_dense LU(nbddl,nbddl);
		// 		lu_factor(J,LU);
		// 		lu_solve(LU,U,F);
		gmres_solve(TT,U,F,40,verbose);
		
		
		////=============================================================////
		////===================== Calcul de l'erreur ====================////
		////=============================================================////
		vect<Cplx> Err, Err2, Norme;
		resize(Err, nbdof);
		resize(Err2,nbdof);
		resize(Norme,nbdof);
		for(int j=0; j<nbdof; j++){
			Err[j] =  U[j]-Ref[j];
		}
		mv_prod(Err2,M,Err);
		mv_prod(Norme,M,Ref);
		
		Real erreur=0;
		Cplx norme =0.;
		Cplx val=0;
		for(int j=0; j<nbdof; j++){
			val += Err2[j]*conj(Err[j]);
			norme  += Norme[j]*conj(Ref[j]);
		}
		erreur=abs(val/norme);
		if (verbose>0){
			std::cout << "erreur:\t" << erreur << std::endl;
		}
		////=============================================================////
		////======================== Sauvegardes ========================////
		////=============================================================////
		std::ofstream output((output_name+"_"+NbrToStr<Real>(p)+".txt").c_str(),std::ios::app);
		if (verbose>0){
			std::cout<<"Output in "<<output_name<<std::endl;
		}
		if (!output){
			std::cerr<<"Output file cannot be created"<<std::endl;
			exit(1);
		}
		else{
			output<<lc<<" "<<erreur<<std::endl;
		}
		output.close();
	}
}

///==================================================================================////
///===========================Second kind Neumann====================================////
///==================================================================================////

void second_kind_neumann_2D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose){
    ////=============================================================////
    ////=======================  Mesh building  =====================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc,verbose);
    
    
    ////=============================================================////
    ////=======================  Mesh loading  ======================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Chargement du maillage"<<std::endl;
    }
    Real kappa=1.;
    
    geometry geom;
    load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());
    
    mesh_1D Omega(geom);
    load_elt_gmsh(Omega,0);
    gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
    
    ////=============================================================////
    ////================== Calcul de la normale =====================////
    ////=============================================================////
    
    nrml_1D n_(Omega);
    
    // 	for (int j=0;j < nb_elt(Omega);j++){
    // 		elt_1D seg = Omega[j];
    // 		R3 G;
    // 		G[0]= 0.5 * ( seg[0][0] + seg[1][0] );
    // 		G[1]= 0.5 * ( seg[0][1] + seg[1][1] );
    // 		G[2]= 0;
    
    // 		cout<<(G,n_[j])<<std::endl;
    // 	}
    
//     swap(n_);
    ////=============================================================////
    ////================ Assemblage de la matrice ===================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Assemblage operateurs integraux"<<std::endl;
    }
    int nbelt = nb_elt(Omega);
    P1_1D dof(Omega);
    int nbdof = nb_dof(dof);
    
    
    gmm_dense V(nbdof,nbdof),T(nbdof,nbdof),M(nbdof,nbdof);
    bem<P1_1D,P1_1D, SLP_2D>   Vop(kappa,n_,n_);
    bem<P1_1D,P1_1D, DLP_2D>   Kop(kappa,n_,n_);
    
    progress bar("assembly", nbelt*nbelt,verbose);
    for(int j=0; j<nbelt; j++){
        const elt_1D& tj = Omega[j];
        const N2&     jj = dof[j];
        
        for(int k=0; k<nbelt; k++,bar++){
            const elt_1D& tk = Omega[k];
            const N2&     kk = dof[k];
            
            V(jj,kk) += Vop (tj,tk);
            T(jj,kk) += Kop (tj,tk);
            
        }
        
        M(jj,jj) +=      MassP1(tj);
        T(jj,jj) += -0.5*MassP1(tj);
    }
    bar.end();
    
	for (int l=0;l<harmonics.size();l++){
		Real p = harmonics[l];
	    ////=============================================================////
        ////================== Harmonique de Fourier ====================////
        ////=============================================================////
        
        vect<Cplx> Ep ;resize(Ep ,nbdof);fill(Ep ,0.);
		vect<Cplx> Ref;resize(Ref,nbdof);fill(Ref,0.);
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
            
            Ep [I] += 0.5*Vinc;
			
			Vinc[0] = (1./kappa)*boost::math::cyl_bessel_j(p,kappa*R)/((p/(kappa*R))*boost::math::cyl_bessel_j(p,kappa*R)-boost::math::cyl_bessel_j(p+1,kappa*R))*exp( iu*p*theta0 );
            Vinc[1] = (1./kappa)*boost::math::cyl_bessel_j(p,kappa*R)/((p/(kappa*R))*boost::math::cyl_bessel_j(p,kappa*R)-boost::math::cyl_bessel_j(p+1,kappa*R))*exp( iu*p*theta1 );
            
			Ref[I] += 0.5*Vinc;
            
        }
        
        
		////=============================================================////
		////================ Assemblage du second membre ================////
		////=============================================================////
		
		vect<Cplx> F; resize(F,nbdof); fill(F,0);
		vect<Cplx> Ftemp; resize(Ftemp,nbdof); fill(Ftemp,0);
		vect<Cplx> gN; resize(gN,nbdof); fill(gN,0);
		
		// Boundary condition
		gN= Ep;
		// 	fill(gD,Cplx(0.));
		
		
		mv_prod(Ftemp,V,gN);
		for(int j=0; j<nbelt; j++){
			F[j] += -Ftemp[j];
		}

		////=============================================================////
		////================ Résolution système linéaire ================////
		////=============================================================////
		if (verbose>0){
			std::cout<<"Appel du solveur"<<std::endl;
		}
		vect<Cplx> U;
		resize(U,nbdof);
		
		// 	gmm_dense LU(nbddl,nbddl);
		// 		lu_factor(J,LU);
		// 		lu_solve(LU,U,F);
		gmres_solve(T,U,F,40,verbose);
		
		
		////=============================================================////
		////===================== Calcul de l'erreur ====================////
		////=============================================================////
		vect<Cplx> Err, Err2, Norme;
		resize(Err, nbdof);
		resize(Err2,nbdof);
		resize(Norme,nbdof);
		for(int j=0; j<nbdof; j++){
			Err[j] =  U[j]-Ref[j];
		}
		mv_prod(Err2,M,Err);
		mv_prod(Norme,M,Ref);
		
		Real erreur=0;
		Cplx norme =0.;
		Cplx val=0;
		for(int j=0; j<nbdof; j++){
			val += Err2[j]*conj(Err[j]);
			norme  += Norme[j]*conj(Ref[j]);
		}
		erreur=abs(val/norme);
		if (verbose>0){
			std::cout << "erreur:\t" << erreur << std::endl;
		}
		////=============================================================////
		////======================== Sauvegardes ========================////
		////=============================================================////
		std::ofstream output((output_name+"_"+NbrToStr<Real>(p)+".txt").c_str(),std::ios::app);
		if (verbose>0){
			std::cout<<"Output in "<<output_name<<std::endl;
		}
		if (!output){
			std::cerr<<"Output file cannot be created"<<std::endl;
			exit(1);
		}
		else{
			output<<lc<<" "<<erreur<<std::endl;
		}
		output.close();
	}
}


///==================================================================================////
///============================Fourier harmonics=====================================////
///==================================================================================////

void fourier_harmonic_2D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name,int verbose){
    
    ////=============================================================////
    ////=======================  Mesh building  =====================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc,verbose);
    
    ////=============================================================////
    ////=======================  Mesh loading  ======================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Chargement du maillage"<<std::endl;
    }
    Real kappa=1.;
    geometry geom;
    load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());
    
    mesh_1D Omega(geom);
    load_elt_gmsh(Omega,0);
    gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
    ////=============================================================////
    ////================== Calcul de la normale =====================////
    ////=============================================================////
    
    nrml_1D n_(Omega);
    
    // 	for (int j=0;j < nb_elt(Omega);j++){
    // 		elt_1D seg = Omega[j];
    // 		R3 G;
    // 		G[0]= 0.5 * ( seg[0][0] + seg[1][0] );
    // 		G[1]= 0.5 * ( seg[0][1] + seg[1][1] );
    // 		G[2]= 0;
    //
    // 		cout<<(G,n_[j])<<endl;
    // 	}
    
    // 	swap(n_);
    ////=============================================================////
    ////================ Assemblage de la matrice ===================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Assemblage operateurs integraux"<<std::endl;
    }
    int nbelt = nb_elt(Omega);
    P1_1D dof(Omega);
    int nbdof = nb_dof(dof);
    
    gmm_dense V(nbdof,nbdof),W(nbdof,nbdof),K(nbdof,nbdof),M(nbdof,nbdof);
    bem<P1_1D,P1_1D, SLP_2D>   Vop (kappa,n_,n_);
    bem<P1_1D,P1_1D, HSP_2D>   Wop (kappa,n_,n_);
    bem<P1_1D,P1_1D, DLP_2D>   Kop (kappa,n_,n_);
    
    progress bar("assembly", nbelt*nbelt,verbose);
    for(int j=0; j<nbelt; j++){
        const elt_1D& tj = Omega[j];
        const N2&     jj = dof[j];
        
        for(int k=0; k<nbelt; k++,bar++){
            const elt_1D& tk = Omega[k];
            const N2&     kk = dof[k];
            
            V (jj,kk) += Vop  (tj,tk);
            W (jj,kk) += Wop  (tj,tk);
            K (jj,kk) += Kop  (tj,tk);
            
        }
        
        M(jj,jj) += MassP1(tj);
    }
    bar.end();
    
    for (int l=0;l<harmonics.size();l++){
        Real p = harmonics[l];
        
        ////=============================================================////
        ////================== Harmonique de Fourier ====================////
        ////=============================================================////
        
        vect<Cplx> Ep;resize(Ep,nbdof);fill(Ep,0.);
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
            
            Ep[I] += 0.5*Vinc;
            
        }
        //        vect<Cplx> Eph;resize(Eph,nbdof);fill(Eph,0.);
        //        mv_prod(Eph,M,Ep);
        
        
        ////=============================================================////
        ////================== Calcul valeurs propres ===================////
        ////=============================================================////
        
        
        ////===================================////
        ////===== Test avec orthogonalité =====////
        
        // 	Cplx val =0;
        // 	Real err =0;
        // 	for(int j=0; j<nbdof; j++){
        // 		val += Eph[j]*conj(Ep[j]);
        // 	}
        //
        // 	err = abs(val-2*M_PI)/(2*M_PI);
        //
        //
        // 	cout<<"PS : "<<err<<endl;
        
        ////===================================////
        ////=========== Test avec V ===========////
        
        Cplx val =0;
        Real err_V =0;
        
        Cplx Bessel_1_p = boost::math::cyl_bessel_j(p,kappa*R);
        Cplx Bessel_2_p = boost::math::cyl_neumann (p,kappa*R);
        Cplx Hankel_1_p = Bessel_1_p+ iu*Bessel_2_p;
        Cplx ref      = iu * M_PI * M_PI * R * R * Hankel_1_p *Bessel_1_p;
        
        vect<Cplx> Temp;resize(Temp,nbdof);fill(Temp,0.);
        mv_prod(Temp,V,Ep);
        for(int j=0; j<nbdof; j++){
            val    += conj(Ep[j])*Temp[j];
        }
        
        err_V = abs(val-ref) /abs(ref);
        if (verbose>0){
            std::cout<<"V : "<<err_V<<std::endl;
        }
        ////===================================////
        ////=========== Test avec W ===========////
        
        val =0;
        Real err_W =0;
        
        Cplx d_Bessel_1_p = (p/(kappa*R))*Bessel_1_p-boost::math::cyl_bessel_j(p+1,kappa*R);
        Cplx d_Bessel_2_p = (p/(kappa*R))*Bessel_2_p-boost::math::cyl_neumann(p+1,kappa*R);
        Cplx d_Hankel_1_p = d_Bessel_1_p+ iu*d_Bessel_2_p;
        ref      = - kappa * kappa * R * R * iu * M_PI * M_PI * d_Hankel_1_p *d_Bessel_1_p;
        
        fill(Temp,0.);
        mv_prod(Temp,W,Ep);
        for(int j=0; j<nbdof; j++){
            val    += conj(Ep[j])*Temp[j];
        }
        
        err_W = abs(val-ref) /abs(ref);
        if (verbose>0){
            std::cout<<"W : "<<err_W<<std::endl;
        }
        ////===================================////
        ////=========== Test avec K ===========////
        
        val =0;
        Real err_K =0;
        ref      = - iu * R * R * kappa  * M_PI * M_PI * (d_Hankel_1_p *Bessel_1_p + Hankel_1_p * d_Bessel_1_p)/2.;
        
        fill(Temp,0.);
        mv_prod(Temp,K,Ep);
        for(int j=0; j<nbdof; j++){
            val    += conj(Ep[j])*Temp[j];
        }
        
        err_K = abs(val-ref) /abs(ref);
        if (verbose>0){
            std::cout<<"K : "<<err_K<<std::endl;
        }
        
        ////=============================================================////
        ////======================== Sauvegardes ========================////
        ////=============================================================////
        std::ofstream output((output_name+"_"+NbrToStr<Real>(p)+".txt").c_str(),std::ios::app);
        if (!output){
            std::cerr<<"Output file cannot be created"<<std::endl;
            exit(1);
        }
        else{
            output<<lc<<" "<<err_V<<" "<<err_W<<" "<<err_K<<std::endl;
        }
        output.close();
    }
    
}

///==================================================================================////
///===========================Plane wave harmonics===================================////
///==================================================================================////

void plane_wave_harmonics_2D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose){
    
    ////=============================================================////
    ////=======================  Mesh building  =====================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc,verbose);
    
    ////=============================================================////
    ////=======================  Mesh loading  ======================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Chargement du maillage"<<std::endl;
    }
    Real kappa=1.;
    geometry geom;
    load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());
    
    mesh_1D Omega(geom);
    load_elt_gmsh(Omega,0);
    gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
    
    ////=============================================================////
    ////================== Calcul de la normale =====================////
    ////=============================================================////
    
    nrml_1D n_(Omega);
    
    // 	for (int j=0;j < nb_elt(Omega);j++){
    // 		elt_1D seg = Omega[j];
    // 		R3 G;
    // 		G[0]= 0.5 * ( seg[0][0] + seg[1][0] );
    // 		G[1]= 0.5 * ( seg[0][1] + seg[1][1] );
    // 		G[2]= 0;
    //
    // 		std::cout<<(G,n_[j])<<std::endl;
    // 	}
    
    // 	swap(n_);
    ////=============================================================////
    ////================ Assemblage de la matrice ===================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Assemblage operateurs integraux"<<std::endl;
    }
    int nbelt = nb_elt(Omega);
    P1_1D dof(Omega);
    int nbdof = nb_dof(dof);
    
    gmm_dense V(nbdof,nbdof),W(nbdof,nbdof),K(nbdof,nbdof),M(nbdof,nbdof);
    bem<P1_1D,P1_1D, SLP_2D>   Vop (kappa,n_,n_);
    bem<P1_1D,P1_1D, HSP_2D>   Wop (kappa,n_,n_);
    bem<P1_1D,P1_1D, DLP_2D>   Kop (kappa,n_,n_);
    
    progress bar("assembly", nbelt*nbelt,verbose);
    for(int j=0; j<nbelt; j++){
        const elt_1D& tj = Omega[j];
        const N2&     jj = dof[j];
        
        for(int k=0; k<nbelt; k++,bar++){
            const elt_1D& tk = Omega[k];
            const N2&     kk = dof[k];
            
            V (jj,kk) += Vop  (tj,tk);
            W (jj,kk) += Wop  (tj,tk);
            K (jj,kk) += Kop  (tj,tk);
            
        }
        
        M(jj,jj) += MassP1(tj);
    }
    bar.end();
    
    for (int l=0;l<harmonics.size();l++){
        Real p=harmonics[l];
        
        ////=============================================================////
        ////========== Onde plane et Harmonique de Fourier ==============////
        ////=============================================================////
        
        R3 dir; dir[0]=sqrt(2)/2.;dir[1]=sqrt(2)/2.;dir[2]=0;
        vect<Cplx> Uinc; resize(Uinc,nbdof); fill(Uinc,0.);
        
        vect<Cplx> Ep;resize(Ep,nbdof);fill(Ep,0.);
        
        
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
            
            Ep[I] += 0.5*Vinc;
            
            
            Vinc[0] = exp( iu*R*cos(theta0) );
            Vinc[1] = exp( iu*R*cos(theta1) );
            
            Uinc[I] += 0.5*Vinc;
            
        }
        vect<Cplx> Eph;resize(Eph,nbdof);fill(Eph,0.);
        mv_prod(Eph,M,Ep);
        
        
        ////=============================================================////
        ////================== Calcul valeurs propres ===================////
        ////=============================================================////
        
        
        ////===================================////
        ////===== Test avec orthogonalité =====////
        
        // 	Cplx val =0;
        // 	Real err =0;
        // 	for(int j=0; j<nbdof; j++){
        // 		val += Eph[j]*conj(Ep[j]);
        // 	}
        //
        // 	err = abs(val-2*M_PI)/(2*M_PI);
        //
        //
        // 	std::cout<<"PS : "<<err<<std::endl;
        
        ////===================================////
        ////=========== Test avec V ===========////
        
        Cplx val =0;
        Real err_V =0;
        
        Cplx Bessel_1_p = boost::math::cyl_bessel_j(p,kappa*R);
        Cplx Bessel_2_p = boost::math::cyl_neumann (p,kappa*R);
        Cplx Hankel_1_p = Bessel_1_p + iu*Bessel_2_p;
        Cplx ref      = pow(iu, p+1) * pow(R,2) * M_PI * M_PI * Hankel_1_p *Bessel_1_p*boost::math::cyl_bessel_j(p,R);
        
        vect<Cplx> Temp;resize(Temp,nbdof);fill(Temp,0.);
        mv_prod(Temp,V,Uinc);
        for(int j=0; j<nbdof; j++){
            val    += conj(Ep[j])*Temp[j];
        }
        
        err_V = abs(val-ref) /abs(ref);
        if (verbose>0){
            std::cout<<"V : "<<err_V<<std::endl;
        }
        ////===================================////
        ////=========== Test avec W ===========////
        val =0;
        Real err_W =0;
        
        Cplx d_Bessel_1_p = (p/(kappa*R))*Bessel_1_p-boost::math::cyl_bessel_j(p+1,kappa*R);
        Cplx d_Bessel_2_p = (p/(kappa*R))*Bessel_2_p-boost::math::cyl_neumann (p+1,kappa*R);
        Cplx d_Hankel_1_p = d_Bessel_1_p+ iu*d_Bessel_2_p;
        ref      = - kappa * kappa * R * R * pow(iu, p+1) * M_PI * M_PI * d_Hankel_1_p *d_Bessel_1_p*boost::math::cyl_bessel_j(p,R);
        
        fill(Temp,0.);
        mv_prod(Temp,W,Uinc);
        for(int j=0; j<nbdof; j++){
            val    += conj(Ep[j])*Temp[j];
        }
        
        err_W = abs(val-ref) /abs(ref);
        if (verbose>0){
            std::cout<<"W : "<<err_W<<std::endl;
        }
        ////===================================////
        ////=========== Test avec K ===========////
        
        val =0;
        Real err_K =0;
        ref      = - pow(iu, p+1) * R * R * kappa  * M_PI * M_PI * (d_Hankel_1_p *Bessel_1_p + Hankel_1_p * d_Bessel_1_p)/2. * boost::math::cyl_bessel_j(p,R);
        
        fill(Temp,0.);
        mv_prod(Temp,K,Uinc);
        for(int j=0; j<nbdof; j++){
            val    += conj(Ep[j])*Temp[j];
        }
        
        err_K = abs(val-ref) /abs(ref);
        if (verbose>0){
            std::cout<<"K : "<<err_K<<std::endl;
        }
        ////=============================================================////
        ////======================== Sauvegardes ========================////
        ////=============================================================////
        std::ofstream output((output_name+"_"+NbrToStr<Real>(p)+".txt").c_str(),std::ios::app);
        if (!output){
            std::cerr<<"Output file cannot be created"<<std::endl;
            exit(1);
        }
        else{
            output<<lc<<" "<<err_V<<" "<<err_W<<" "<<err_K<<std::endl;
        }
        output.close();
    }
}

///==================================================================================////
///===============================Champs rayonné=====================================////
///==================================================================================////

void champs_rayonne_2D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose){
	////=============================================================////
    ////=======================  Mesh building  =====================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc,verbose);
	gmsh_disc  (("disc_"+NbrToStr(lc)).c_str(),R,lc,verbose);
    
    ////=============================================================////
    ////=======================  Mesh loading  ======================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Chargement du maillage"<<std::endl;
    }
    Real kappa=1.;
    
    geometry geom,vol;
    load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());
	std::cout<<meshfile(geom)<<std::endl;
    load_node_gmsh(vol ,("disc_"+NbrToStr(lc)).c_str());
	std::cout<<get_node(geom,10)<<std::endl;
	std::cout<<get_node(vol,50)<<std::endl;
	std::cout<<meshfile(geom)<<std::endl;
	std::cout<<meshfile(vol)<<std::endl;
	
    mesh_1D Omega(geom);
	
    load_elt_gmsh(Omega,0);
	std::cout<<"ok"<<std::endl;
	std::cout<<Omega[0]<<std::endl;
//     gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
// 	gmsh_clean(("disc_"+NbrToStr(lc)).c_str());
    
    ////=============================================================////
    ////================== Calcul de la normale =====================////
    ////=============================================================////
    std::cout<<"ok"<<std::endl;
    nrml_1D n_(Omega);
	std::cout<<"ok"<<std::endl;
	////=============================================================////
    ////================== Calcul de la normale =====================////
    ////=============================================================////
    
	
	
	
	
	
}
