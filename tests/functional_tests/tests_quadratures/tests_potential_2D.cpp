#include <bemtool2/tools.h>
#include <bemtool2/htool_wrap.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

using namespace std;
using namespace bemtool;
// using namespace htool;

int main(int argc, char* argv[]){

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
	std::vector<Real> harmonics(3);harmonics[0]=1;harmonics[1]=2;harmonics[2]=3;
  Real lc_s = 0.005;
	Real lc_v = 0.005;
	Real R=0.5;
  htool::SetNdofPerElt(1);
	htool::SetEpsilon(1e-6);
	htool::SetEta(0.1);
  htool::SetMinClusterSize(10);
  ////=======================  Mesh building  =====================////
  if (rankWorld==0){
    std::cout<<"Construction du maillage"<<std::endl;
    gmsh_circle(("circle_"+NbrToStr(lc_s)).c_str(),R,lc_s);
    gmsh_disc  (("disc_"+NbrToStr(lc_v)).c_str(),R*0.5,lc_v);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  ////=======================  Mesh loading  ======================////
  if (rankWorld==0){
  std::cout<<"Chargement du maillage"<<std::endl;
  }

  geometry geom,vol;
  load_node_gmsh(geom,("circle_"+NbrToStr(lc_s)).c_str());
  load_node_gmsh(vol ,("disc_"+NbrToStr(lc_v)).c_str());

  mesh_1D Omega(geom);
  mesh_2D Vol(vol);

  load_elt_gmsh(Omega,0);
  load_elt_gmsh(Vol,0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rankWorld==0){
    gmsh_clean(("circle_"+NbrToStr(lc_s)).c_str());
    gmsh_clean(("disc_"+NbrToStr(lc_v)).c_str());
  }

  ////================== Calcul de la normale =====================////


  nrml_1D n_(Omega);


  ////============ Matrice pour le champs rayonne =================////

  const vect<R3>& node = get_node(vol);
  int nbpt  = size(node);

  const vect<R3>& node_g = get_node(geom);

  int nbelt = nb_elt(Omega);
  P1_1D dof(Omega);
  int nbdof = nb_dof(dof);

  Real kappa=10;

  potential<P1_1D,SLP_2D> SLPop(kappa,n_);
  potential<P1_1D,DLP_2D> DLPop(kappa,n_);
  MyPotential<P1_1D,SLP_2D> vmat_SLP(SLPop,vol,nbpt,nbdof);
  MyPotential<P1_1D,DLP_2D> vmat_DLP(DLPop,vol,nbpt,nbdof);

  std::vector<int> tabt(nbpt);
  std::vector<int> tabs(nbdof);
  std::vector<htool::R3> xt(nbpt);
  std::vector<htool::R3> xs(nbdof);

  for (int i=0;i<nbpt;i++){
    tabt[i]=i;
    htool::R3 pt;pt[0]=node[i][0];pt[1]=node[i][1];pt[2]=node[i][2];
    xt[i]=pt;
  }

  for (int i=0;i<nbdof;i++){
    tabs[i]=i;
    htool::R3 pt;pt[0]=node_g[i][0];pt[1]=node_g[i][1];pt[2]=node_g[i][2];
    xs[i]=pt;
  }

  htool::HMatrix<htool::fullACA,Cplx> hmat_SLP(vmat_SLP,xt,tabt,xs,tabs,1);
  htool::HMatrix<htool::fullACA,Cplx> hmat_DLP(vmat_DLP,xt,tabt,xs,tabs,1);

  for (int l=0;l<1;l++){
    Real p = harmonics[l];

    ////======================= Trace solution ======================////

    vect<Cplx> TraceDirichlet;resize(TraceDirichlet,nbdof);fill(TraceDirichlet,0.);
    vect<Cplx> TraceNeumann;  resize(TraceNeumann,nbdof);  fill(TraceNeumann,0.);

    for (int j=0 ; j<nbelt ; j++){
      const elt_1D& seg = Omega[j];
      const N2&     I   = dof[j];

      R3 X0; X0[0] =  seg[0][0];X0[1]=seg[0][1]; X0[2]=0;
      R3 X1; X1[0] =  seg[1][0];X1[1]=seg[1][1]; X1[2]=0;

      Real theta0 = std::atan2 (X0[1],X0[0]);
      Real theta1 = std::atan2 (X1[1],X1[0]);

      C2 Vinc;

      Vinc[0] = exp( iu*p*theta0 );
      Vinc[1] = exp( iu*p*theta1 );

      TraceDirichlet [I] += 0.5*Vinc;


      Vinc[0] *= kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R));
      Vinc[1] *= kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R));

      TraceNeumann[I] += 0.5*Vinc;

    }


    ////=================== Solution de référence ===================////

    vect<Cplx> Ref;resize(Ref,nbpt);fill(Ref,0.);
    for (int j=0;j<nbpt;j++){
      R3 X; X[0] = node[j][0];X[1]= node[j][1]; X[2]= node[j][2];;

      Real r     = std::sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
      Real theta = std::atan2 (X[1],X[0]);

      Ref[j]= boost::math::cyl_bessel_j(p,kappa*r)/boost::math::cyl_bessel_j(p,kappa*R)*exp( iu*p*theta );
    }


    ////================ Construction de la solution ================////

    vect<Cplx> S;resize(S,nbpt);fill(S,0.);
    vect<Cplx> D;resize(D,nbpt);fill(D,0.);
    vect<Real> Out;resize(Out,nbpt);fill(Out,0.);
    vect<Real> Out_ref;resize(Out_ref,nbpt);fill(Out_ref,0.);
    vect<Cplx> error;resize(error,nbpt);
    Real error2=0;
    Real norm = 0;

    hmat_SLP.mvprod_global(&TraceNeumann[0],&S[0]);
    hmat_DLP.mvprod_global(&TraceDirichlet[0],&D[0]);
    for (int i=0;i<nbpt;i++){
      Out[i]= std::abs(S[i]+D[i]);
      Out_ref[i]= std::abs(Ref[i]);
      error[i]=S[i]+D[i]-Ref[i];
      error2+=pow(std::abs(error[i]),2);
      norm += pow(std::abs(Ref[i]),2);
    }
    error2=std::sqrt(error2)/std::sqrt(norm);
    if (rankWorld==0){
      cout<<error2<<endl;
    }

    cout << "rank : "<<rankWorld<<" "<<error2 <<" "<<hmat_SLP.compression()<<" "<<hmat_DLP.compression()<<" "<<error2<<endl;
    hmat_SLP.print_stats();
    hmat_DLP.print_stats();

    if (rankWorld==0){
    write_gmsh(Vol,Out,"test");
    write_gmsh(Vol,Out_ref,"ref");
    }
  }

  // Finalize the MPI environment.
  MPI_Finalize();
  return test;
}
