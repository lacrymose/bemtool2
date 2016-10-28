#include <bemtool2/tools.h>
#include <math.h>       /* atan */
#include <boost/math/special_functions/bessel.hpp>

using namespace std;
using namespace boost::math;

int main(int argc, char *argv[]){
	
////=============================================================////
////===========================  Input ==========================////
////=============================================================////
	
	// Check the number of parameters
	if (argc < 2) {
		// Tell the user how to run the program
		cerr << "Usage: " << argv[0] << " harmonique finesse rayon output_name" << endl;
		/* "Usage messages" are a conventional way of telling the user
			* how to run a program if they enter the command incorrectly.
			*/
		return 1;
	}
	
	Real p  = StrToReal(argv[1]);
	Real lc = StrToReal(argv[2]);
	Real R  = StrToReal(argv[3]);
	string output_name = argv[4];
	
////=============================================================////
////=======================  Mesh building  =====================////
////=============================================================////
	cout<<"Construction du maillage"<<endl;
	
	gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc);
	
////=============================================================////
////=======================  Mesh loading  ======================////
////=============================================================////
	cout<<"Chargement du maillage"<<endl;
	Real kappa=1.;
	geometry geom;
	load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());
	const vect<R3>& node = get_node(geom);  
	int nb_node = size(node);    

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
	cout<<"Assemblage operateurs integraux"<<endl;
	int nbelt = nb_elt(Omega);
	P1_1D dof(geom); cout<<"attaching"<<endl;dof.attach_to(Omega);
	int nbdof = nb_dof(dof);
	cout<<"attaching finished"<<endl;
	gmm_dense V(nbdof,nbdof),W(nbdof,nbdof),K(nbdof,nbdof),M(nbdof,nbdof);
	cout<<"V"<<endl;
	bem<P1_1D,P1_1D, SLP_2D>   Vop (kappa,n_,n_);
// 	bem<P1_1D,P1_1D, HSP_2D>   Wop (kappa,n_,n_);
// 	bem<P1_1D,P1_1D, DLP_2D>   Kop (kappa,n_,n_);
	cout<<"Debut boucle"<<endl;	
	progress bar("assembly", nbelt*nbelt);
	for(int j=0; j<nbelt; j++){
		const elt_1D& tj = Omega[j];
		const N2&     jj = dof[j];
		
		for(int k=0; k<nbelt; k++,bar++){
			const elt_1D& tk = Omega[k];
// 			const N2&     kk = dof[k];

// 			V (jj,kk) += Vop  (tj,tk);
// 			W (jj,kk) += Wop  (tj,tk);
// 			K (jj,kk) += Kop  (tj,tk);
			
		}
		
// 		M(jj,jj) += MassP1(tj);
	}
	bar.end();
// ////=============================================================////
// ////========== Onde plane et Harmonique de Fourier ==============////
// ////=============================================================////
// 	
// // 	R3 dir; dir[0]=sqrt(2)/2.;dir[1]=sqrt(2)/2.;dir[2]=0;
// 	vect<Cplx> Uinc; resize(Uinc,nbdof); fill(Uinc,0.);
// 	
// 	vect<Cplx> Ep;resize(Ep,nbdof);fill(Ep,0.);
// 	
// 	
// 	for (int j=0 ; j<nbelt ; j++){
// 		const elt_1D& seg = Omega[j];
// 		const N2&     I   = dof[j];
// 		
// 		R3 X0; X0[0] =  seg[0][0];X0[1]=seg[0][1]; X0[2]=0;
// 		R3 X1; X1[0] =  seg[1][0];X1[1]=seg[1][1]; X1[2]=0;
// 		
// 		Real theta0 = atan (X0[1]/X0[0]);
// 		Real theta1 = atan (X1[1]/X1[0]);
// 		
// 		if (X0[0]<0 & X0[1]>=0){
// 			theta0 += M_PI;
// 		}
// 		if (X0[0]<0 & X0[1]<0){
// 			theta0 -= M_PI;
// 		}
// 		if (X1[0]<0 & X1[1]>=0){
// 			theta1 += M_PI;
// 		}
// 		if (X1[0]<0 & X1[1]<0){
// 			theta1 -= M_PI;
// 		}
// 
// 		C2 Vinc;
// 		
// 		Vinc[0] = exp( iu*p*theta0 );
// 		Vinc[1] = exp( iu*p*theta1 );
// 
// 		Ep[I] += 0.5*Vinc;
// 		
// 		
// 		Vinc[0] = exp( iu*R*cos(theta0) );
// 		Vinc[1] = exp( iu*R*cos(theta1) );
// 			
// 		Uinc[I] += 0.5*Vinc;
// 		
// // 		cout<<theta0*180/M_PI<<" "<<theta1*180/M_PI<<endl;
// 	}
// 	vect<Cplx> Eph;resize(Eph,nbdof);fill(Eph,0.);
// 	mv_prod(Eph,M,Ep);
// 	
// 	
// ////=============================================================////
// ////================== Calcul valeurs propres ===================////
// ////=============================================================////
// 	
// 	
// 	////===================================////
// 	////===== Test avec orthogonalité =====////
// 		
// // 	Cplx val =0;
// // 	Real err =0;
// // 	for(int j=0; j<nbdof; j++){
// // 		val += Eph[j]*conj(Ep[j]);
// // 	}
// // 	
// // 	err = abs(val-2*M_PI)/(2*M_PI);
// // 	
// // 	
// // 	cout<<"PS : "<<err<<endl;
// 	
// 	////===================================////
// 	////=========== Test avec V ===========////
// 	
// 	Cplx val =0;
// 	Real err_V =0;
// 
// 	Cplx Bessel_1_p_kappaR = cyl_bessel_j(p,kappa*R);
// 	Cplx Bessel_2_p_kappaR = cyl_neumann (p,kappa*R);
// 	Cplx Hankel_1_p_kappaR = Bessel_1_p_kappaR + iu*Bessel_2_p_kappaR;
// 	Cplx ref      = pow(iu, p+1) * pow(R,2) * M_PI * M_PI * Hankel_1_p_kappaR *Bessel_1_p_kappaR*cyl_bessel_j(p,R);
// 	
// 	vect<Cplx> Temp;resize(Temp,nbdof);fill(Temp,0.);
// 	mv_prod(Temp,V,Uinc);
// 	for(int j=0; j<nbdof; j++){
// 		val    += conj(Ep[j])*Temp[j];
// 	}
// 	
// 	err_V = abs(val-ref) /abs(ref);
// 	cout<<"V : "<<err_V<<endl;
// // 	////===================================////
// // 	////=========== Test avec W ===========////
// // 	
// // 	val =0;
// // 	Real err_W =0;
// // 	
// // 	Cplx d_Bessel_1_p = (p/kappa)*Bessel_1_p-cyl_bessel_j(p+1,kappa);
// // 	Cplx d_Bessel_2_p = (p/kappa)*Bessel_2_p-cyl_neumann(p+1,kappa);;
// // 	Cplx d_Hankel_1_p = d_Bessel_1_p+ iu*d_Bessel_2_p;
// // 	ref      = - kappa * kappa * iu * M_PI * M_PI * d_Hankel_1_p *d_Bessel_1_p;
// // 
// // 	fill(Temp,0.);
// // 	mv_prod(Temp,W,Ep);
// // 	for(int j=0; j<nbdof; j++){
// // 		val    += conj(Ep[j])*Temp[j];
// // 	}
// // 	
// // 	err_W = abs(val-ref) /abs(ref);
// // // 	cout<<"W : "<<err_W<<endl;
// // 	////===================================////
// // 	////=========== Test avec K ===========////
// // 	
// // 	val =0;
// // 	Real err_K =0;
// // 	ref      = - iu * kappa  * M_PI * M_PI * (d_Hankel_1_p *Bessel_1_p + Hankel_1_p * d_Bessel_1_p)/2.;
// // 	
// // 	fill(Temp,0.);
// // 	mv_prod(Temp,K,Ep);
// // 	for(int j=0; j<nbdof; j++){
// // 		val    += conj(Ep[j])*Temp[j];
// // 	}
// // 	
// // 	err_K = abs(val-ref) /abs(ref);
// // // 	cout<<"K : "<<err_K<<endl;
	
}