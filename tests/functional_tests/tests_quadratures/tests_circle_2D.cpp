#include <bemtool2/tools.h>
#include <bemtool2/gmm_wrap.h>
#include <gmm/gmm.h>
#include <mpi.h>

using namespace bemtool;
using namespace std;

int main(int argc, char const *argv[]) {
  ////====================== Initialization MPI ===================////
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int sizeWorld;
  MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);

  // Get the rank of the process
  int rankWorld;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);

  ////======================= Some variables  =====================////
  bool test =0;
  Real kappa=1.;
  Real lc = 0.1;
  Real R=1;
  std::vector<Real> harmonics(3);harmonics[0]=1;harmonics[1]=2;harmonics[2]=3;
  // htool::SetNdofPerElt(1);
	// htool::SetEpsilon(1e-6);
	// htool::SetEta(0.1);
  // htool::SetMinClusterSize(10);


  ////=======================  Mesh building  =====================////
  if (rankWorld==0){
      std::cout<<"Construction du maillage"<<std::endl;
  }
  gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc);
  MPI_Barrier(MPI_COMM_WORLD);

  ////=======================  Mesh loading  ======================////
  if (rankWorld==0){
  std::cout<<"Chargement du maillage"<<std::endl;
  }

  geometry geom;
  load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());

  mesh_1D Omega(geom);
  load_elt_gmsh(Omega,0);
  gmsh_clean(("circle_"+NbrToStr(lc)).c_str());

  ////================== Calcul de la normale =====================////

  nrml_1D n_(Omega);


  ////================ Assemblage de la matrice ===================////
  if (rankWorld==0){
      std::cout<<"Assemblage operateurs integraux"<<std::endl;
  }
  int nbelt = nb_elt(Omega);
  P1_1D dof(Omega);
  int nbdof = nb_dof(dof);


  gmm_dense V(nbdof,nbdof),K(nbdof,nbdof),M(nbdof,nbdof);
  bem<P1_1D,P1_1D, SLP_PH_2D>   Vop(kappa,n_,n_);
  bem<P1_1D,P1_1D, DLP_PH_2D>   Kop(kappa,n_,n_);

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

  for (int l=0;l<1;l++){
    Real p = harmonics[l];

    ////================== Harmonique de Fourier ====================////

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

    ////================ Assemblage du second membre ================////

    vect<Cplx> F; resize(F,nbdof); fill(F,0);
    vect<Cplx> Ftemp; resize(Ftemp,nbdof); fill(Ftemp,0);
    vect<Cplx> gD; resize(gD,nbdof); fill(gD,0);

    // Boundary condition
    gD=Ep;
    //fill(gD,Cplx(0.));

    mv_prod(F,K,gD);
    mv_prod(Ftemp,M,gD);
    cout << "nbdof : "<<nbdof<< endl;
    for(int j=0; j<nbdof; j++){
        F[j] = 0.5*Ftemp[j] - F[j];
        cout << F[j] << endl;
    }

    ////================ Résolution système linéaire ================////

    if (rankWorld==0){
    std::cout<<"Appel du solveur"<<std::endl;
    }
    vect<Cplx> U;
    resize(U,nbdof);

    // 	gmm_dense LU(nbddl,nbddl);
    // 		lu_factor(J,LU);
    // 		lu_solve(LU,U,F);
    gmres_solve(V,U,F,40);



    ////===================== Calcul de l'erreur ====================////
    vect<Cplx> Err, Err2, Norme;
    resize(Err, nbdof);
		resize(Err2,nbdof);
		resize(Norme,nbdof);
		for(int j=0; j<nbdof; j++){
      cout << U[j] <<" "<<Ref[j]<< endl;
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
		if (rankWorld==0){
			std::cout << "erreur:\t" << erreur << std::endl;
		}
		////=============================================================////
		////======================== Sauvegardes ========================////
		////=============================================================////
		// std::ofstream output((output_name+"_"+NbrToStr<Real>(p)+".txt").c_str(),std::ios::app);
		// if (verbose>0){
		// 	std::cout<<"Output in "<<output_name<<std::endl;
		// }
		// if (!output){
		// 	std::cerr<<"Output file cannot be created"<<std::endl;
		// 	exit(1);
		// }
		// else{
		// 	output<<lc<<" "<<erreur<<std::endl;
		// }
		// output.close();
	}


  return test;
}
