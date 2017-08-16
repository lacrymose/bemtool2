#include <bemtool/tools.h>
#include <htool/htool.hpp>


int main(int argc, char const *argv[]) {
  ////==== Initialization MPI
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int sizeWorld;
  MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);

  // Get the rank of the process
  int rankWorld;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);

  ////=== Some variables
  bool test =0;
  Real kappa=1.;
  Real lc = 0.1;
  Real R=1;
  std::vector<Real> harmonics(3);harmonics[0]=1;harmonics[1]=2;harmonics[2]=3;

  // Htool variables
  htool::SetNdofPerElt(1);
	htool::SetEpsilon(1e-6);
	htool::SetEta(10);
  htool::SetMinClusterSize(30);
	htool::SetMaxBlockSize(10000000);

  // HPDDM setting
  HPDDM::Option& opt = *HPDDM::Option::get();
	opt.parse(argc, argv, rankWorld == 0);
	if(rankWorld != 0)
		opt.remove("verbosity");


  ////=== Mesh
  if (rankWorld==0){
    gmsh_circle(("circle",R,lc);
  }

  geometry geom;
  load_node_gmsh(geom,("circle");

  mesh_1D Omega(geom);
  load_elt_gmsh(Omega,0);
  gmsh_clean(("circle");

  nrml_1D n_(Omega);
  int nbelt = nb_elt(Omega);
  P1_1D dof(Omega);
  int nbdof = nb_dof(dof);

  ////=== Assembling
  htool::Matrix M(nbdof,nbdof);
  for(int j=0; j<nbelt; j++){
    const elt_1D& tj = Omega[j];
    const N2&     jj = dof[j];
    for (int k=0;k<;k++){
      for (int l=0;l<;l++){
        M(jj[k],jj[l])
      }
    }
    M(jj,jj) += MassP1(tj);
  }

  bem<P1_1D,P1_1D, SLP_2D>   Vop(kappa,n_,n_);
  modified_bem<P1_1D,P1_1D, DLP_2D>   M_Kop(kappa,n_,n_,-0.5);

  std::vector<int> tab(nbdof);
  std::vector<htool::R3> x(nbdof);


  for (int i =0 ; i<nbelt;i++){
  	for (int j =0;j<dof.dim_loc;j++){
  		x[dof[i][j]][0]=Omega[i][j][0];
  		x[dof[i][j]][1]=Omega[i][j][1];
  		x[dof[i][j]][2]=Omega[i][j][2];
  	}
  }

  MyMatrix V(Vop,nbdof);
  MyMatrix K(Kop,nbdof);
  htool::HMatrix<htool::partialACA,Cplx> HV(V,x,tab);
  htool::HMatrix<htool::partialACA,Cplx> HK(M_Kop,x,tab);

  ////=== Analytical solution
  for (int l=0;l<harmonics.size();l++){
    Real p = harmonics[0];
    std::vector<Cplx> Ep(nbdof,0), Ref(nbdof,0);
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

      for (int i=0;i<2;i++){
        Ep [I[i]] += 0.5*Vinc[i];
      }

      Vinc[0] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta0 );
      Vinc[1] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta1 );

      for (int i=0;i<2;i++){
        Ref [I[i]] += 0.5*Vinc[i];
      }

    }

    ////=== Right-hand side

    std::vector<Cplx> f_global(nbdof,0),Ftemp(nbdof,0),gD(nbdof,0);
    // Boundary condition
    gD=Ep;
    //fill(gD,Cplx(0.));
    HK.mvprod_global(gD,f_global);

    mv_prod(Ftemp,M,gD);

    for(int j=0; j<nbdof; j++){
        f_global[j] = - f_global[j];
    }

    ////=== Solve
    std::vector<Cplx> U(nbdof);
    int n_local= HA.get_local_size_cluster();

    std::vector<complex<double>> x_ref_local(n_local),x_local(n_local,0),f_local(n_local);
    for (int i=0;i<n_local;i++){
    	f_local[i]=f_global[MasterClusters[rankWorld][i]];
      x_ref_local[i]=Ref[MasterClusters[rankWorld][i]];
    }


    htool::HPDDMOperator<htool::partialACA,Cplx> A_HPDDM(HA,P);
    htool::Identity<Cplx> P(cluster_to_ovr_subdomain.size());
    complex<double>* const f = &(f_local[0]);
    complex<double>* sol = &(x_local[0]);
    HPDDM::IterativeMethod::solve(A_HPDDM, f, sol, 1,HA.get_comm());

    //=== Error
    double error2=0;
    for (int i=0;i<x_local.size();i++){
      error2+=pow(std::abs(x_local[i]-x_ref_local[i]),2);
    }
    MPI_Allreduce(MPI_IN_PLACE, &error2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (rankWorld==0){
      cout << "Test "<<sqrt(error2)<<endl;
    }
  }


  return test;
}
