#include <bemtool2/tools.h>
#include <htool/htool.hpp>

using namespace bemtool;
using namespace std;

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
  Real lc = 0.005;
  Real R=2;
  std::vector<Real> harmonics(3);harmonics[0]=1;harmonics[1]=2;harmonics[2]=3;

  // Htool variables
  htool::SetNdofPerElt(1);
	htool::SetEpsilon(1e-6);
	htool::SetEta(100);
  htool::SetMinClusterSize(10);
	htool::SetMaxBlockSize(10000000);

  // HPDDM setting
  HPDDM::Option& opt = *HPDDM::Option::get();
	opt.parse(argc, argv, rankWorld == 0);
	if(rankWorld != 0)
		opt.remove("verbosity");


  ////=== Mesh
  if (rankWorld==0){
    gmsh_circle("circle",R,lc);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  geometry geom;
  load_node_gmsh(geom,"circle");

  mesh_1D Omega(geom);
  load_elt_gmsh(Omega,0);
  MPI_Barrier(MPI_COMM_WORLD);

  if (rankWorld==0){
    gmsh_clean("circle");
  }


  nrml_1D n_(Omega);
  int nbelt = nb_elt(Omega);
  P1_1D dof(Omega);
  int nbdof = nb_dof(dof);

  ////=== Assembling
  bem<P1_1D,P1_1D, SLP_PH_2D>   Vop(kappa,n_,n_);
  modified_bem<P1_1D,P1_1D, DLP_PH_2D>   M_Kop(kappa,n_,n_,-0.5);

  std::vector<int> tab(nbdof);
  std::vector<htool::R3> x(nbdof);


  for (int i =0 ; i<nbelt;i++){
  	for (int j =0;j<dof.dim_loc;j++){
  		x[dof[i][j]][0]=Omega[i][j][0];
  		x[dof[i][j]][1]=Omega[i][j][1];
  		x[dof[i][j]][2]=Omega[i][j][2];
  	}
  }
  for (int i=0;i<nbdof;i++){
    tab[i]=i;
  }


  HMatrix<bem<P1_1D,P1_1D, SLP_PH_2D>> V(Vop,nbdof);
  HMatrix<modified_bem<P1_1D,P1_1D, DLP_PH_2D>> MK(M_Kop,nbdof);
  htool::HMatrix<htool::partialACA,Cplx> HV(V,x,tab);
  htool::HMatrix<htool::partialACA,Cplx> HMK(MK,x,tab);

  ////=== Analytical solution
  for (int l=0;l<1;l++){
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

      Vinc[0] = (1./kappa)*boost::math::cyl_bessel_j(p,kappa*R)/((p/(kappa*R))*boost::math::cyl_bessel_j(p,kappa*R)-boost::math::cyl_bessel_j(p+1,kappa*R))*exp( iu*p*theta0 );
      Vinc[1] = (1./kappa)*boost::math::cyl_bessel_j(p,kappa*R)/((p/(kappa*R))*boost::math::cyl_bessel_j(p,kappa*R)-boost::math::cyl_bessel_j(p+1,kappa*R))*exp( iu*p*theta1 );

      for (int i=0;i<2;i++){
        Ref [I[i]] += 0.5*Vinc[i];
      }

    }

    ////=== Right-hand side

    std::vector<Cplx> f_global(nbdof,0),gN(nbdof,0);
    // Boundary condition
    gN=Ep;
    //fill(gD,Cplx(0.));
    HV.mvprod_global(gN.data(),f_global.data());

    // cout << "nbdof : "<<nbdof<< endl;
    for(int j=0; j<nbdof; j++){
        f_global[j] = - f_global[j];
        // cout << f_global[j] << endl;

    }

    ////=== Solve
    const std::vector<std::vector<int>>& MasterClusters=HMK.get_MasterClusters();
    int n_local= HMK.get_local_size_cluster();

    std::vector<Cplx> x_ref_local(n_local),x_local(n_local,0),f_local(n_local);
    for (int i=0;i<n_local;i++){
    	f_local[i]=f_global[MasterClusters[rankWorld][i]];
      x_ref_local[i]=Ref[MasterClusters[rankWorld][i]];
    }

    htool::Identity<Cplx> P(n_local);
    htool::HPDDMOperator<htool::partialACA,Cplx> A_HPDDM(HMK,P);

    complex<double>* const f = &(f_local[0]);
    complex<double>* sol = &(x_local[0]);
    HPDDM::IterativeMethod::solve(A_HPDDM, f, sol, 1,HMK.get_comm());

		//=== Local to global
		std::vector<Cplx> x_global(nbdof);
		std::vector<Cplx> rcv(nbdof);


		std::vector<int> recvcounts(sizeWorld);
		std::vector<int>  displs(sizeWorld);

		displs[0] = 0;

		for (int i=0; i<sizeWorld; i++) {
			recvcounts[i] = MasterClusters[i].size();
			if (i > 0)
				displs[i] = displs[i-1] + recvcounts[i-1];
		}


		MPI_Allgatherv(&(x_local.front()), recvcounts[rankWorld], MPI_DOUBLE_COMPLEX, &(rcv.front()), &(recvcounts.front()), &(displs.front()), MPI_DOUBLE_COMPLEX, HMK.get_comm());

		for (int i=0; i<sizeWorld; i++)
		for (int j=0; j< MasterClusters[i].size(); j++)
			x_global[MasterClusters[i][j]] = rcv[displs[i]+j];

    //=== Error
    // Global mass matrix
    htool::Matrix<Cplx> Mass(nbdof,nbdof);

		for(int j=0; j<nbelt; j++){
			const elt_1D& tj = Omega[j];
			const N2&     jj = dof[j];
			mat<2,2,Real> Melt=MassP1(tj);
			for (int i =0;i<2;i++){
				for (int j =0;j<2;j++){
					Mass(jj[i],jj[j]) +=Melt(i,j);
				}
			}
		}

    // Computation
    std::vector<Cplx> Err(nbdof), Err2(nbdof), Norm(nbdof);
    for (int i=0;i<x_global.size();i++){
      // cout << x_local[i] <<" "<<x_ref_local[i]<<endl;
      Err[i] = x_global[i]-Ref[i];
    }
    Err2=Mass*Err;
    Norm=Mass*x_global;

    Cplx val=0;
    Cplx norm =0.;
    double error2=0;
    for(int j=0; j<nbdof; j++){
      val += Err2[j]*conj(Err[j]);
      norm  += Norm[j]*conj(Ref[j]);
    }
  	error2=abs(val/norm);
		double compression = HMK.compression();
		int nb_lrmat = HMK.get_nlrmat();
		int nb_dmat  = HMK.get_ndmat();
		// MPI_Barrier(MPI_COMM_WORLD);
    HMK.print_stats();
		if (rankWorld==0){
			cout << "Error: "<<sqrt(error2)<<endl;
			cout << "Compression rate : "<<compression<<endl;
			cout << "nbr lr : "<<nb_lrmat<<endl;
			cout << "nbr dense : "<<nb_dmat<<endl;

    }
  }


  return test;
}
