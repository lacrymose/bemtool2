#include <bemtool2/tools.h>
#include <htool/htool.hpp>
#include <boost/math/special_functions.hpp>

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
  Real lc = 0.01;
  Real R=0.5;
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
  // BEM matrices with hmatrix
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


  bem<P1_1D,P1_1D, SLP_PH_2D>   Vop (kappa,n_,n_);
  bem<P1_1D,P1_1D, HSP_PH_2D>   Wop (kappa,n_,n_);
  bem<P1_1D,P1_1D, DLP_PH_2D>   Kop (kappa,n_,n_);

  HMatrix<bem<P1_1D,P1_1D, SLP_PH_2D>> V(Vop,nbdof);
  HMatrix<bem<P1_1D,P1_1D, HSP_PH_2D>> W(Wop,nbdof);
  HMatrix<bem<P1_1D,P1_1D, DLP_PH_2D>> K(Kop,nbdof);

  htool::HMatrix<htool::partialACA,Cplx> HV(V,x,tab);
  htool::HMatrix<htool::partialACA,Cplx> HW(W,x,tab);
  htool::HMatrix<htool::partialACA,Cplx> HK(K,x,tab);


  // Mass matrix
  htool::Matrix<Cplx>  M(nbdof,nbdof);

  for(int j=0; j<nbelt; j++){
      const elt_1D& tj = Omega[j];
      const N2&     jj = dof[j];
      mat<2,2,Real> Melt=MassP1(tj);
      for (int i =0;i<2;i++){
        for (int j =0;j<2;j++){
          M(jj[i],jj[j]) +=Melt(i,j);
        }
      }
  }

  // For each harmonic
  for (int l=0;l<harmonics.size();l++){
    Real p = harmonics[l];

    // Eigenvector
    std::vector<Cplx> Ep(nbdof,0);

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

      Ep[I[0]] += 0.5*Vinc[0];
      Ep[I[1]] += 0.5*Vinc[1];


    }


    Cplx norm;
    Cplx val_test=0;
    std::vector<Cplx> temp(nbdof);
    temp = M*Ep;
    for(int j=0; j<nbdof; j++){
        val_test    += conj(Ep[j])*temp[j];
    }
    norm = val_test;
    cout << "TeST : "<<norm<<endl;

    //=== Tests with V

    Cplx val =0;
    Real err_V =0;

    Cplx Bessel_1_p = boost::math::cyl_bessel_j(p,kappa*R);
    // Cplx Bessel_2_p = boost::math::cyl_neumann (p,kappa*R);
    // Cplx Hankel_1_p = Bessel_1_p+ iu*Bessel_2_p;
    Cplx Hankel_1_p = boost::math::cyl_hankel_1(p, kappa*R);
    Cplx ref    = iu * M_PI * M_PI * R * R * Hankel_1_p *Bessel_1_p;

    std::vector<Cplx> V_eigenvector(nbdof);
    HV.mvprod_global(Ep.data(),V_eigenvector.data());
    V_eigenvector=M*V_eigenvector;
    for(int j=0; j<nbdof; j++){
        val    += conj(Ep[j])*V_eigenvector[j];
    }

    err_V = abs(sqrt(val/norm)-ref);
    std::cout<<"V : "<<err_V<<std::endl;


        //
        // ////===================================////
        // ////=========== Test avec W ===========////
        //
        // val =0;
        // Real err_W =0;
        //
        // Cplx d_Bessel_1_p = (p/(kappa*R))*Bessel_1_p-boost::math::cyl_bessel_j(p+1,kappa*R);
        // Cplx d_Bessel_2_p = (p/(kappa*R))*Bessel_2_p-boost::math::cyl_neumann(p+1,kappa*R);
        // Cplx d_Hankel_1_p = d_Bessel_1_p+ iu*d_Bessel_2_p;
        // ref      = - kappa * kappa * R * R * iu * M_PI * M_PI * d_Hankel_1_p *d_Bessel_1_p;
        //
        // fill(Temp,0.);
        // mv_prod(Temp,W,Ep);
        // for(int j=0; j<nbdof; j++){
        //     val    += conj(Ep[j])*Temp[j];
        // }
        //
        // err_W = abs(val-ref) /abs(ref);
        // if (verbose>0){
        //     std::cout<<"W : "<<err_W<<std::endl;
        // }
        // ////===================================////
        // ////=========== Test avec K ===========////
        //
        // val =0;
        // Real err_K =0;
        // ref      = - iu * R * R * kappa  * M_PI * M_PI * (d_Hankel_1_p *Bessel_1_p + Hankel_1_p * d_Bessel_1_p)/2.;
        //
        // fill(Temp,0.);
        // mv_prod(Temp,K,Ep);
        // for(int j=0; j<nbdof; j++){
        //     val    += conj(Ep[j])*Temp[j];
        // }
        //
        // err_K = abs(val-ref) /abs(ref);
        // if (verbose>0){
        //     std::cout<<"K : "<<err_K<<std::endl;
        // }
        //
        // ////=============================================================////
        // ////======================== Sauvegardes ========================////
        // ////=============================================================////
        // std::ofstream output((output_name+"_"+NbrToStr<Real>(p)+".txt").c_str(),std::ios::app);
        // if (!output){
        //     std::cerr<<"Output file cannot be created"<<std::endl;
        //     exit(1);
        // }
        // else{
        //     output<<lc<<" "<<err_V<<" "<<err_W<<" "<<err_K<<std::endl;
        // }
        // output.close();
  }
}
