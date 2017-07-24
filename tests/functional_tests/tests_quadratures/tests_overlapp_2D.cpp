#include <bemtool2/tools.h>
#include <bemtool2/overlap.hpp>
#include <set>
#include <mpi.h>
#include <htool/htool.hpp>


using namespace bemtool;
using namespace std;



class MyMatrix: public htool::IMatrix<complex<double>>{
	bem<P1_1D,P1_1D, SLP_2D>& Vop;

public:
	MyMatrix(bem<P1_1D,P1_1D, SLP_2D>& Vop0, int nbdof0 ):IMatrix(nbdof0,nbdof0),Vop(Vop0){}

	Cplx get_coef(const int& i, const int& j)const {return Vop(i,j);}

 	htool::SubMatrix<Cplx> get_submatrix(const std::vector<int>& J, const std::vector<int>& K) const{
 		htool::SubMatrix<Cplx> submat(J,K);
 		Vop(J,K,submat);
 		return submat;
 	}

};

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

	// HPDDM verbosity
	HPDDM::Option& opt = *HPDDM::Option::get();
	opt.parse(argc, argv, rankWorld == 0);
	if(rankWorld != 0)
		opt.remove("verbosity");

  ////======================= Some variables  =====================////
  bool test =0;
  Real kappa=1.;
  Real lc = 0.01;
  Real R=0.5;
  std::vector<Real> harmonics(3);harmonics[0]=1;harmonics[1]=2;harmonics[2]=3;
  // htool::SetNdofPerElt(1);
	// htool::SetEpsilon(1e-6);
	htool::SetEta(0.01);
  // htool::SetMinClusterSize(10);


  ////=======================  Mesh building  =====================////
  if (rankWorld==0){
      std::cout<<"Construction du maillage"<<std::endl;
      gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  ////=======================  Mesh loading  ======================////
  if (rankWorld==0){
  std::cout<<"Chargement du maillage"<<std::endl;
  }

  geometry geom;
  load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());

  mesh_1D Omega(geom);
  load_elt_gmsh(Omega,0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rankWorld==0){
    // gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
  }
  ////================== Calcul de la normale =====================////

  nrml_1D n_(Omega);


  ////================ Assemblage de la matrice ===================////
  if (rankWorld==0){
      std::cout<<"Assemblage operateurs integraux"<<std::endl;
  }
  int nbelt = nb_elt(Omega);
  P1_1D dof(Omega);
  int nbdof = nb_dof(dof);
  // const vect<R3>& node = get_node(geom);
  // int nbpt = size(node);

	// htool::Matrix<Cplx> V(nbdof,nbdof),K(nbdof,nbdof);
	bem<P1_1D,P1_1D, SLP_2D>   Vop(kappa,n_,n_);
	// bem<P1_1D,P1_1D, DLP_2D>   Kop(kappa,n_,n_);
	//
	// progress bar("assembly", nbelt*nbelt);
	// for(int j=0; j<nbelt; j++){
	// 		const elt_1D& tj = Omega[j];
	// 		const N2&     jj = dof[j];
	//
	// 		for(int k=0; k<nbelt; k++,bar++){
	// 				const elt_1D& tk = Omega[k];
	// 				const N2&     kk = dof[k];
	// 				mat<2,2, Cplx> vmat;
	// 				mat<2,2, Cplx>	kmat;
	// 				vmat = Vop (tj,tk);
	// 				kmat = Kop (tj,tk);
	//
	// 				for (int j=0;j<2;j++){
	// 					for (int k=0;k<2;k++){
	// 						V(jj[j],kk[k])+= vmat(j,k);
	// 						K(jj[j],kk[k])+= kmat(j,k);
	// 					}
	// 				}
	//
	//
	// 		}
	//
	// 		// M(jj,jj) += MassP1(tj);
	// }
	// bar.end();

  std::vector<int> tab(nbdof);
  std::vector<htool::R3> x(nbdof);

// if (rankWorld==0){
// 	for(int j=0; j<nbpt; j++){
//     cout << j<<" "<<dof[j][0]<<" "<<dof[j][1]<< endl;
// 	}
// }
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
	MyMatrix V(Vop,nbdof);
  htool::HMatrix<htool::fullACA,complex<double>> HA(V,x,tab);
MPI_Barrier(MPI_COMM_WORLD);
HA.print_stats();
double compression = HA.compression();
int nb_lrmat = HA.get_nlrmat();
int nb_dmat  = HA.get_ndmat();
if (rankWorld==0){
	cout << "Compression rate : "<<compression<<endl;
	cout << "nbr lr : "<<nb_lrmat<<endl;
	cout << "nbr dense : "<<nb_dmat<<endl;
}
  ////================ Partition ===================////

  const std::vector<std::vector<int>>& MasterClusters=HA.get_MasterClusters();


	std::vector<int> cluster_to_ovr_subdomain;
	std::vector<int> ovr_subdomain_to_global;
	std::vector<int> neighbors;
	std::vector<std::vector<int> > intersections;

	Partition(MasterClusters,dof,cluster_to_ovr_subdomain,ovr_subdomain_to_global,neighbors,intersections);
	// int count=0;
	// cout << "proc : "<<rankWorld<<endl;
	// for(auto iter=neighbors.begin(); iter!=neighbors.end();++iter) {
	// 	// if (rankWorld==0){
	//
	// 		cout << "voisin : "<<*iter<<" | ";
	// 		for (int i=0;i<intersections[count].size();i++){
	// 			cout << intersections[count][i]<<" ";
	// 		}
	// 		cout<<endl;
	// 	// }
	// 	count+=1;
  // }
	// cout << "=================="<<endl;


	std::vector<int> part_overlap(nbdof,0);
	for (int i=0;i<ovr_subdomain_to_global.size();i++){
		part_overlap[ovr_subdomain_to_global[i]]=1;
	}
	for (int i =0;i<cluster_to_ovr_subdomain.size();i++){
		part_overlap[ovr_subdomain_to_global[cluster_to_ovr_subdomain[i]]]+=1;
	}
	write_gmsh_2(Omega,dof,part_overlap,"part_ovlerap_"+NbrToStr(rankWorld));

	htool::Preconditioner<Cplx> P(
		V,ovr_subdomain_to_global,cluster_to_ovr_subdomain,neighbors,intersections);

	P.num_fact();

	// Global vectors
	std::vector<complex<double>> x_ref(nbdof,1),f_global(nbdof,0);
	HA.mvprod_global(x_ref.data(),f_global.data());

	// Local vectors
	int n_local= HA.get_local_size_cluster();
	std::vector<complex<double>> x_ref_local(n_local,1),x_local(n_local,0),f_local(n_local,1);
	for (int i=0;i<n_local;i++){
		f_local[i]=f_global[MasterClusters[rankWorld][i]];
	}

	// Solve
  htool::HPDDMOperator<htool::fullACA,complex<double>> A_HPDDM(HA,P);
  complex<double>* const f = &(f_local[0]);
  complex<double>* sol = &(x_local[0]);
  HPDDM::IterativeMethod::solve(A_HPDDM, f, sol, 1,HA.get_comm());

	double error2=0;
	for (int i=0;i<x_local.size();i++){
		error2+=pow(std::abs(x_local[i]-x_ref_local[i]),2);
	}
	MPI_Allreduce(MPI_IN_PLACE, &error2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if (rankWorld==0){
		cout << "Test "<<sqrt(error2)<<endl;
	}

  // gmm_dense V(nbdof,nbdof),K(nbdof,nbdof),M(nbdof,nbdof);
  // bem<P1_1D,P1_1D, SLP_2D>   Vop(kappa,n_,n_);
  // bem<P1_1D,P1_1D, DLP_2D>   Kop(kappa,n_,n_);
  //
  // progress bar("assembly", nbelt*nbelt);
  // for(int j=0; j<nbelt; j++){
  //   const elt_1D& tj = Omega[j];
  //   const N2&     jj = dof[j];
  //
  //   for(int k=0; k<nbelt; k++,bar++){
  //     const elt_1D& tk = Omega[k];
  //     const N2&     kk = dof[k];
  //
  //     V(jj,kk) += Vop (tj,tk);
  //     K(jj,kk) += Kop (tj,tk);
  //
  //   }
  //
  //   M(jj,jj) += MassP1(tj);
  // }
  // bar.end();
  //
  // for (int l=0;l<harmonics.size();l++){
  //   Real p = harmonics[l];
  //
  //   ////================== Harmonique de Fourier ====================////
  //
  //   vect<Cplx> Ep ;resize(Ep ,nbdof);fill(Ep ,0.);
  //   vect<Cplx> Ref;resize(Ref,nbdof);fill(Ref,0.);
  //   for (int j=0 ; j<nbelt ; j++){
  //     const elt_1D& seg = Omega[j];
  //     const N2&     I   = dof[j];
  //
  //     R3 X0; X0[0] =  seg[0][0];X0[1]=seg[0][1]; X0[2]=0;
  //     R3 X1; X1[0] =  seg[1][0];X1[1]=seg[1][1]; X1[2]=0;
  //
  //     Real theta0 = atan (X0[1]/X0[0]);
  //     Real theta1 = atan (X1[1]/X1[0]);
  //
  //     if (X0[0]<0 & X0[1]>=0){
  //         theta0 += M_PI;
  //     }
  //     if (X0[0]<0 & X0[1]<0){
  //         theta0 -= M_PI;
  //     }
  //     if (X1[0]<0 & X1[1]>=0){
  //         theta1 += M_PI;
  //     }
  //     if (X1[0]<0 & X1[1]<0){
  //           theta1 -= M_PI;
  //     }
  //
  //     C2 Vinc;
  //
  //     Vinc[0] = exp( iu*p*theta0 );
  //     Vinc[1] = exp( iu*p*theta1 );
  //
  //     Ep [I] += 0.5*Vinc;
  //
  //     Vinc[0] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta0 );
  //     Vinc[1] = kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R))*exp( iu*p*theta1 );
  //
  //     Ref[I] += 0.5*Vinc;
  //
  //   }
  //
  //   ////================ Assemblage du second membre ================////
  //
  //   vect<Cplx> F; resize(F,nbdof); fill(F,0);
  //   vect<Cplx> Ftemp; resize(Ftemp,nbdof); fill(Ftemp,0);
  //   vect<Cplx> gD; resize(gD,nbdof); fill(gD,0);
  //
  //   // Boundary condition
  //   gD=Ep;
  //   //fill(gD,Cplx(0.));
  //
  //   mv_prod(F,K,gD);
  //   mv_prod(Ftemp,M,gD);
  //
  //   for(int j=0; j<nbelt; j++){
  //       F[j] = 0.5*Ftemp[j] - F[j];
  //   }
  //
  //   ////================ Résolution système linéaire ================////
  //
  //   if (rankWorld==0){
  //   std::cout<<"Appel du solveur"<<std::endl;
  //   }
  //   vect<Cplx> U;
  //   resize(U,nbdof);
  //
  //   // 	gmm_dense LU(nbddl,nbddl);
  //   // 		lu_factor(J,LU);
  //   // 		lu_solve(LU,U,F);
  //   gmres_solve(V,U,F,40);
  //
  //
  //
  //   ////===================== Calcul de l'erreur ====================////
  //   vect<Cplx> Err, Err2, Norme;
  //   resize(Err, nbdof);
	// 	resize(Err2,nbdof);
	// 	resize(Norme,nbdof);
	// 	for(int j=0; j<nbdof; j++){
	// 		Err[j] =  U[j]-Ref[j];
	// 	}
	// 	mv_prod(Err2,M,Err);
	// 	mv_prod(Norme,M,Ref);
  //
	// 	Cplx val=0;
	// 	Cplx norme =0.;
	// 	Real erreur=0;
	// 	for(int j=0; j<nbdof; j++){
	// 		val += Err2[j]*conj(Err[j]);
	// 		norme  += Norme[j]*conj(Ref[j]);
	// 	}
	// 	erreur=abs(val/norme);
	// 	if (rankWorld==0){
	// 		std::cout << "erreur:\t" << erreur << std::endl;
	// 	}
	// 	////=============================================================////
	// 	////======================== Sauvegardes ========================////
	// 	////=============================================================////
	// 	// std::ofstream output((output_name+"_"+NbrToStr<Real>(p)+".txt").c_str(),std::ios::app);
	// 	// if (verbose>0){
	// 	// 	std::cout<<"Output in "<<output_name<<std::endl;
	// 	// }
	// 	// if (!output){
	// 	// 	std::cerr<<"Output file cannot be created"<<std::endl;
	// 	// 	exit(1);
	// 	// }
	// 	// else{
	// 	// 	output<<lc<<" "<<erreur<<std::endl;
	// 	// }
	// 	// output.close();
	// }

MPI_Barrier(MPI_COMM_WORLD);
  return test;
}
