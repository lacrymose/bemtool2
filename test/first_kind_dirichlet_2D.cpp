#include "tests_circle_2D.hpp"

void first_kind_dirichlet(Real lc, string output_name){

////=============================================================////
////=======================  Mesh building  =====================////
////=============================================================////
	cout<<"Construction du maillage"<<endl;
	gmm_circle(("circle_"+NbrToStr(lc)+".geo").c_str(),1,lc);
	
	
////=============================================================////
////=======================  Mesh loading  ======================////
////=============================================================////
	cout<<"Chargement du maillage"<<endl;
	Real kappa=1.;
	
	load_node(("circle_"+NbrToStr(lc)+".msh").c_str());
	const vect<R3>& node = get_node();  
	int nb_node = size(node);    

	mesh_1D Omega;  
	load(Omega,0);
	gmm_clean(("circle_"+NbrToStr(lc)).c_str());
////=============================================================////
////================== Calcul de la normale =====================////
////=============================================================////
	cout<<"Calcul de la normale"<<endl;
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
	
	swap(n_);
////=============================================================////
////================ Assemblage de la matrice ===================////
////=============================================================////
	cout<<"Assemblage operateurs integraux"<<endl;
	int nbelt = nb_elt(Omega);
	P1_1D dof; dof.attach_to(Omega);
	int nbdof = nb_dof(dof);
	
	
	gmm_dense V(nbdof,nbdof),K(nbdof,nbdof),M(nbdof,nbdof);
	bem<P1_1D,P1_1D, SLP_2D>   Vop(kappa,n_,n_);
	bem<P1_1D,P1_1D, DLP_2D>   Kop(kappa,n_,n_);
	
	progress bar("assembly", nbelt*nbelt);
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

////=============================================================////
////================ Assemblage du second membre ================////
////=============================================================////
	cout<<"Assemblage du second membre"<<endl;
	
	R3 dir; dir[0]=sqrt(2)/2.;dir[1]=sqrt(2)/2.;dir[2]=0;
	
	vect<Cplx> Udinc; resize(Udinc,nbdof); fill(Udinc,0.);
	vect<Cplx> Uninc; resize(Uninc,nbdof); fill(Uninc,0.);
	
	for (int j=0 ; j<nbelt ; j++){
		const elt_1D& seg = Omega[j];
		const N2&     I   = dof[j];
		
		R3 X0; X0[0] =  seg[0][0];X0[1]=seg[0][1]; X0[2]=0;
		R3 X1; X1[0] =  seg[1][0];X1[1]=seg[1][1]; X1[2]=0;
		C2 Vinc;
		
		Vinc[0] = exp( iu*kappa*(X0,dir) );
		Vinc[1] = exp( iu*kappa*(X1,dir) );

		Udinc[I] += 0.5*Vinc;
		
		Vinc[0] = (dir,n_[j])*iu*kappa*exp( iu*kappa*(X0,dir) );
		Vinc[1] = (dir,n_[j])*iu*kappa*exp( iu*kappa*(X1,dir) );
		
		Uninc[I] += 0.5*Vinc;
		
	}

	vect<Cplx> F; resize(F,nbdof); fill(F,0);
	vect<Cplx> Ftemp; resize(Ftemp,nbdof); fill(Ftemp,0);
	vect<Cplx> gD; resize(gD,nbdof); fill(gD,0);
	
	// Boundary condition
	gD=Udinc;
// 	fill(gD,Cplx(0.));
	
	mv_prod(F,K,gD);
	mv_prod(Ftemp,M,gD);
	
	for(int j=0; j<nbelt; j++){
		F[j] = 0.5*Ftemp[j] - F[j];
	}
	
	// Source term
// 	mv_prod(Ftemp,K,Udinc);
// 	for(int j=0; j<nbelt; j++){
// 		F[j] += Ftemp[j];
// 	}
// 	mv_prod(Ftemp,M,Udinc);
// 	for(int j=0; j<nbelt; j++){
// 		F[j] -= 0.5*Ftemp[j];
// 	}
// 	mv_prod(Ftemp,V,Uninc);
// 	for(int j=0; j<nbelt; j++){
// 		F[j] += Ftemp[j];
// 	}
	mv_prod(Ftemp,M,Udinc);
	for(int j=0; j<nbelt; j++){
		F[j] -= Ftemp[j];
	}
////=============================================================////
////================ Résolution système linéaire ================////
////=============================================================////
	cout<<"Appel du solveur"<<endl;
	
	vect<Cplx> U;
	resize(U,nbdof);

// 	gmm_dense LU(nbddl,nbddl);
// 		lu_factor(J,LU);
// 		lu_solve(LU,U,F);
	gmres_solve(V,U,F,40);

	
////=============================================================////
////===================== Calcul de l'erreur ====================////
////=============================================================////
	vect<Cplx> Err, Err2, Norme;
	resize(Err, nbdof);
	resize(Err2,nbdof);
	resize(Norme,nbdof);
	for(int j=0; j<nbdof; j++){
		Err[j] =  U[j]-Uninc[j];
	}
	mv_prod(Err2,M,Err);
	mv_prod(Norme,M,Uninc);

	Cplx val=0;
	Cplx norme =0.;
	Real erreur=0;
	for(int j=0; j<nbdof; j++){
		val += Err2[j]*conj(Err[j]);
		norme  += Norme[j]*conj(Uninc[j]);
	}
	erreur=abs(val/norme);
	cout << "erreur:\t" << erreur << endl;
	
////=============================================================////
////======================== Sauvegardes ========================////
////=============================================================////
	ofstream output(output_name.c_str(),ios::app);
	cout<<"Output in "<<output_name<<endl;
	if (!output){
		cerr<<"Output file cannot be created"<<endl;
		exit(1);
	}
	else{
		output<<lc<<" "<<erreur<<endl;
	}
	output.close();
	
}
