#include "tests_potential.hpp"

///==================================================================================////
///===============================Champs rayonné=====================================////
///==================================================================================////

void potential_elt_2D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose){
    ////=============================================================////
    ////=======================  Mesh building  =====================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc,verbose);
    gmsh_disc  (("disc_"+NbrToStr(lc)).c_str(),R*0.9,lc,verbose);

    ////=============================================================////
    ////=======================  Mesh loading  ======================////
    ////=============================================================////
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

    gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
    gmsh_clean(("disc_"+NbrToStr(lc)).c_str());

    ////=============================================================////
    ////================== Calcul de la normale =====================////
    ////=============================================================////

    nrml_1D n_(Omega);

    ////=============================================================////
    ////============ Matrice pour le champs rayonne =================////
    ////=============================================================////
    const vect<R3>& node = get_node(vol);
    int nbpt  = size(node);

    int nbelt = nb_elt(Omega);
    P1_1D dof(Omega);
    int nbdof = nb_dof(dof);

    Real kappa=10;

    gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
    potential<P1_1D,SLP_2D> SLPop(kappa,n_);
    potential<P1_1D,CST_2D> DLPop(kappa,n_);

    progress bar("assembly", nbpt*nbelt,verbose);
    for (int j=0; j<nbpt ;j++){
        const N1& jj = j;
        for (int k=0;k<nbelt;k++, bar++){
            const elt_1D& tk = Omega[k];
            const N2&     kk = dof[k];

            SLP (jj,kk) += SLPop(node[j],tk) ;
			DLP (jj,kk) += DLPop(node[j],tk) ;
		}
		
    }
    bar.end();
	
	for (int l=0;l<harmonics.size();l++){
		Real p = harmonics[l];
		////=============================================================////
		////======================= Trace solution ======================////
		////=============================================================////
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

		////=============================================================////
		////=================== Solution de référence ===================////
		////=============================================================////
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

		////=============================================================////
		////================ Construction de la solution ================////
		////=============================================================////
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
		write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc)+"_ref");

	}
}


void potential_node_2D(std::vector<Real> harmonics, Real lc, Real R, std::string output_name, int verbose){
    ////=============================================================////
    ////=======================  Mesh building  =====================////
    ////=============================================================////
    if (verbose>0){
        std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc,verbose);
    gmsh_disc  (("disc_"+NbrToStr(lc)).c_str(),R*0.9,lc,verbose);

    ////=============================================================////
    ////=======================  Mesh loading  ======================////
    ////=============================================================////
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

    gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
    gmsh_clean(("disc_"+NbrToStr(lc)).c_str());

    ////=============================================================////
    ////================== Calcul de la normale =====================////
    ////=============================================================////

    nrml_1D n_(Omega);

    ////=============================================================////
    ////============ Matrice pour le champs rayonne =================////
    ////=============================================================////
    const vect<R3>& node = get_node(vol);
    int nbpt  = size(node);

    int nbelt = nb_elt(Omega);
    P1_1D dof(Omega);
    int nbdof = nb_dof(dof);

    Real kappa=10;

    gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
    potential<P1_1D,SLP_2D> SLPop(kappa,n_);
    potential<P1_1D,CST_2D> DLPop(kappa,n_);

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
		////=============================================================////
		////======================= Trace solution ======================////
		////=============================================================////
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

		////=============================================================////
		////=================== Solution de référence ===================////
		////=============================================================////
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

		////=============================================================////
		////================ Construction de la solution ================////
		////=============================================================////
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
		write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc)+"_ref");

	}
}