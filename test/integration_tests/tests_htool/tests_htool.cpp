#include <bemtool2/tools.h>
#include <mpi.h>
#include <bemtool2/htool_wrap.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

using namespace bemtool;

///==================================================================================////
///============================ potential_hmat_2D ===================================////
///==================================================================================////

void potential_hmat_2D(int argc, char* argv[], std::vector<Real> harmonics, Real lc_s, Real lc_v, Real R, std::string output_name, Real Epsilon, Real Eta, int MinClusterSize, int MaxBlockSize, int verbose=0){

  ////=====================  Setting for htool  ===================////

  MPI_Init(&argc, &argv);
  int rankWorld, sizeWorld;
  MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);
  htool::SetNdofPerElt(1);
  htool::SetEpsilon(Epsilon);
  htool::SetEta(Eta);
  htool::SetMinClusterSize(MinClusterSize);
  htool::SetMaxBlockSize(MaxBlockSize);
  
  ////=======================  Mesh building  =====================////

  if (rankWorld==0){
    if (verbose>0){
      std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc_s)).c_str(),R,lc_s,verbose);
    gmsh_disc  (("disc_"+NbrToStr(lc_v)).c_str(),R*0.9,lc_v,verbose);
  }

  ////=======================  Mesh loading  ======================////

  MPI_Barrier(MPI_COMM_WORLD);
  if (verbose>0){
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
std::cout << nbdof <<" "<< nbpt << " " <<nbdof * nbpt << std::endl;
  Real kappa=10;

  potential<P1_1D,SLP_2D> SLPop(kappa,n_);
  potential<P1_1D,DLP_2D> DLPop(kappa,n_);
  MyMatrix<P1_1D,SLP_2D> vmat_SLP(SLPop,vol,nbpt,nbdof);
  MyMatrix<P1_1D,DLP_2D> vmat_DLP(DLPop,vol,nbpt,nbdof);


  std::vector<htool::Real> rt(nbpt);
  std::vector<htool::Real> rs(nbdof);
  std::vector<int> tabt(nbpt);
  std::vector<int> tabs(nbdof);

  std::vector<htool::R3> xt(nbpt);
  std::vector<htool::R3> xs(nbdof);

  for (int i=0;i<nbpt;i++){
    tabt[i]=i;
    htool::R3 pt;pt[0]=node[i][0];pt[1]=node[i][1];pt[2]=node[i][2];
    xt[i]=pt;
    rt[i]=0;
  }

  for (int i=0;i<nbdof;i++){
    tabs[i]=i;
    htool::R3 pt;pt[0]=node_g[i][0];pt[1]=node_g[i][1];pt[2]=node_g[i][2];
    xs[i]=pt;
    rs[i]=0;
  }

  htool::HMatrix hmat_SLP(vmat_SLP,xt,rt,tabt,xs,rs,tabs,1);
  htool::HMatrix hmat_DLP(vmat_DLP,xt,rt,tabt,xs,rs,tabs,1);


  for (int l=0;l<harmonics.size();l++){
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

    vect<Cplx> Sh;resize(Sh,nbpt);fill(Sh,0.);
    vect<Cplx> Dh;resize(Dh,nbpt);fill(Dh,0.);
    vect<Real> Out;resize(Out,nbpt);fill(Out,0.);
    vect<Real> Out_ref;resize(Out_ref,nbpt);fill(Out_ref,0.);

    std::pair <double,double > mvp_stats_SLP = MvProdMPI(Sh,hmat_SLP,TraceNeumann);
    std::pair <double,double > mvp_stats_DLP =MvProdMPI(Dh,hmat_DLP,TraceDirichlet);

    for (int i=0;i<nbpt;i++){
      Out[i]= std::abs(Sh[i]+Dh[i]);
      Out_ref[i]= std::abs(Ref[i]);
    }

	add_stats(hmat_SLP,"MvProd_SLP (mean)",std::get<0>(mvp_stats_SLP));
	add_stats(hmat_SLP,"MvProd_SLP (max)",std::get<1>(mvp_stats_SLP));
// 	add_stats(hmat_SLP,"MvProd err",norm(ua-ub)/norm(ua));
	add_stats(hmat_SLP,"Compression_SLP",CompressionRate(hmat_SLP));
	add_stats(hmat_SLP,"Nb dense mats SLP",nb_densemats(hmat_SLP));
	add_stats(hmat_SLP,"Nb lr mats SLP",nb_lrmats(hmat_SLP));
	add_stats(hmat_SLP,"Relative Frob error SLP",sqrt(squared_absolute_error(hmat_SLP, vmat_SLP))/NormFrob(vmat_SLP));

	add_stats(hmat_DLP,"MvProd_SLP (mean)",std::get<0>(mvp_stats_DLP));
	add_stats(hmat_DLP,"MvProd_SLP (max)",std::get<1>(mvp_stats_DLP));
// 	add_stats(hmat_DLP,"MvProd err",norm(ua-ub)/norm(ua));
	add_stats(hmat_DLP,"Compression_SLP",CompressionRate(hmat_DLP));
	add_stats(hmat_DLP,"Nb dense mats SLP",nb_densemats(hmat_DLP));
	add_stats(hmat_DLP,"Nb lr mats SLP",nb_lrmats(hmat_DLP));
	add_stats(hmat_DLP,"Relative Frob error SLP",sqrt(squared_absolute_error(hmat_DLP, vmat_DLP))/NormFrob(vmat_DLP));

	
	std::cout<< "################# Harmonique n°"<<p<<" ###################"<<std::endl;
	std::cout<< "## SLP"<<std::endl; 
	print_stats(hmat_SLP);
	std::cout<< "## DLP"<<std::endl;
	print_stats(hmat_DLP);
	std::cout << std::endl;
	
	
    if (rankWorld==0){
    write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_hmat");
    write_gmsh(Vol,Out_ref,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_ref");
    }
  }
  MPI_Finalize();
}

///==================================================================================////
///=========================== potential_hmat_3D ====================================////
///==================================================================================////

void potential_hmat_3D(int argc, char* argv[], std::vector<Real> harmonics, Real lc_s, Real lc_v, Real R, std::string output_name, Real Epsilon, Real Eta, int MinClusterSize, int MaxBlockSize, int verbose=0){
  ////=====================  Setting for htool  ===================////

  MPI_Init(&argc, &argv);
  int rankWorld, sizeWorld;
  MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);
  htool::SetNdofPerElt(1);
  htool::SetEpsilon(Epsilon);
  htool::SetEta(Eta);
  htool::SetMinClusterSize(MinClusterSize);
  htool::SetMaxBlockSize(MaxBlockSize);
  
  ////=======================  Mesh building  =====================////
  if (rankWorld==0){
    if (verbose>0){
      std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_sphere(("sphere_"+NbrToStr(lc_s)).c_str(),R,lc_s,verbose);
    gmsh_ball  (("ball_"+NbrToStr(lc_v)).c_str(),R*0.5,lc_v,verbose);
  }

  ////=======================  Mesh loading  ======================////

  if (verbose>0){
    std::cout<<"Chargement du maillage"<<std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  geometry geom,vol;
  load_node_gmsh(geom,("sphere_"+NbrToStr(lc_s)).c_str());
  load_node_gmsh(vol ,("ball_"+NbrToStr(lc_v)).c_str());

  mesh_2D Omega(geom);
  mesh_3D Vol(vol);

  load_elt_gmsh(Omega,0);
  load_elt_gmsh(Vol,0);
  MPI_Barrier(MPI_COMM_WORLD);

  if (rankWorld==0){
    gmsh_clean(("sphere_"+NbrToStr(lc_s)).c_str());
    gmsh_clean(("ball_"+NbrToStr(lc_v)).c_str());
  }

  ////================== Calcul de la normale =====================////


  nrml_2D n_(Omega);


  ////============ Matrice pour le champs rayonne =================////

  const vect<R3>& node = get_node(vol);
  const vect<R3>& node_g = get_node(geom);
  int nbpt  = size(node);

  int nbelt = nb_elt(Omega);
  P1_2D dof(Omega);
  int nbdof = nb_dof(dof);
std::cout << nbdof <<" "<< nbpt << " " <<nbdof * nbpt << std::endl;
  Real kappa=10;

  potential<P1_2D,SLP_3D> SLPop(kappa,n_);
  potential<P1_2D,DLP_3D> DLPop(kappa,n_);
  MyMatrix<P1_2D,SLP_3D> vmat_SLP(SLPop,vol,nbpt,nbdof);
  MyMatrix<P1_2D,DLP_3D> vmat_DLP(DLPop,vol,nbpt,nbdof);


  std::vector<htool::Real> rt(nbpt);
  std::vector<htool::Real> rs(nbdof);
  std::vector<int> tabt(nbpt);
  std::vector<int> tabs(nbdof);

  std::vector<htool::R3> xt(nbpt);
  std::vector<htool::R3> xs(nbdof);

  for (int i=0;i<nbpt;i++){
    tabt[i]=i;
    htool::R3 pt;pt[0]=node[i][0];pt[1]=node[i][1];pt[2]=node[i][2];
    xt[i]=pt;
    rt[i]=0;
  }

  for (int i=0;i<nbdof;i++){
    tabs[i]=i;
    htool::R3 pt;pt[0]=node_g[i][0];pt[1]=node_g[i][1];pt[2]=node_g[i][2];
    xs[i]=pt;
    rs[i]=0;
  }


  htool::HMatrix hmat_SLP(vmat_SLP,xt,rt,tabt,xs,rs,tabs);
  htool::HMatrix hmat_DLP(vmat_DLP,xt,rt,tabt,xs,rs,tabs);

  for (int l=0;l<harmonics.size();l++){
    Real p = harmonics[l];

    ////======================= Trace solution ======================////

    vect<Cplx> TraceDirichlet;resize(TraceDirichlet,nbdof);fill(TraceDirichlet,0.);
    vect<Cplx> TraceNeumann;  resize(TraceNeumann,nbdof);  fill(TraceNeumann,0.);


    for (int j=0 ; j<nbelt ; j++){
      const elt_2D& tri = Omega[j];
      const N3&     I   = dof[j];

      R3 X0; X0[0] =  tri[0][0];X0[1]=tri[0][1]; X0[2]=tri[0][2];
      R3 X1; X1[0] =  tri[0][0];X1[1]=tri[0][1]; X1[2]=tri[0][2];
      R3 X2; X2[0] =  tri[0][0];X2[1]=tri[0][1]; X2[2]=tri[0][2];


      Real radius0 = sqrt(X0[0] * X0[0] + X0[1] * X0[1] + X0[2] * X0[2]);
      Real radius1 = sqrt(X1[0] * X1[0] + X1[1] * X1[1] + X1[2] * X1[2]);
      Real radius2 = sqrt(X2[0] * X2[0] + X2[1] * X2[1] + X2[2] * X2[2]);

      Real theta0 = std::acos(X0[2]/radius0);
      Real theta1 = std::acos(X1[2]/radius1);
      Real theta2 = std::acos(X2[2]/radius2);

      Real phi0 = std::atan2(X0[1],X0[0]);
      Real phi1 = std::atan2(X1[1],X1[0]);
      Real phi2 = std::atan2(X2[1],X2[0]);

      if (phi0<0){
        phi0 += 2 * M_PI;
      }
      if (phi1<0){
        phi1 += 2 * M_PI;
      }
      if (phi2<0){
        phi2 += 2 * M_PI;
      }

      C3 Vinc;
      Vinc[0] = boost::math::spherical_harmonic(p,p,theta0,phi0);
      Vinc[1] = boost::math::spherical_harmonic(p,p,theta1,phi1);
      Vinc[2] = boost::math::spherical_harmonic(p,p,theta2,phi2);

      TraceDirichlet [I] =  Vinc;

      Cplx derivee;
      derivee= kappa*((p/(kappa*R)) - boost::math::sph_bessel(p+1,kappa*R)/boost::math::sph_bessel(p,kappa*R));

      TraceNeumann[I] =  TraceDirichlet [I]*derivee;


    }

    ////=================== Solution de référence ===================////

    vect<Cplx> Ref;resize(Ref,nbpt);fill(Ref,0.);
    for (int j=0;j<nbpt;j++){
      R3 X; X[0] = node[j][0];X[1]= node[j][1]; X[2]= node[j][2];;

      Real r     = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
      Real theta;
      if (r<1e-10)
      theta=0; // to avoid nan value for points near the origin
      else
      theta = std::acos(X[2]/r);
      Real phi = std::atan2(X[1],X[0]);
      if (phi<0){
        phi += 2 * M_PI;
      }
      Ref[j]= boost::math::sph_bessel(p,kappa*r)/boost::math::sph_bessel(p,kappa*R)*boost::math::spherical_harmonic(p,p,theta,phi);

    }

    ////================ Construction de la solution ================////

    vect<Cplx> Sh;resize(Sh,nbpt);fill(Sh,0.);
    vect<Cplx> Dh;resize(Dh,nbpt);fill(Dh,0.);
    vect<Real> Out;resize(Out,nbpt);fill(Out,0.);
    vect<Real> Out_ref;resize(Out_ref,nbpt);fill(Out_ref,0.);

    std::pair <double,double > mvp_stats_SLP = MvProdMPI(Sh,hmat_SLP,TraceNeumann);
    std::pair <double,double > mvp_stats_DLP = MvProdMPI(Dh,hmat_DLP,TraceDirichlet);

    for (int i=0;i<nbpt;i++){
      Out[i]= std::abs(Sh[i]+Dh[i]);
      Out_ref[i]= std::abs(Ref[i]);
    }


	add_stats(hmat_SLP,"MvProd_SLP (mean)",std::get<0>(mvp_stats_SLP));
	add_stats(hmat_SLP,"MvProd_SLP (max)",std::get<1>(mvp_stats_SLP));
// 	add_stats(hmat_SLP,"MvProd err",norm(ua-ub)/norm(ua));
	add_stats(hmat_SLP,"Compression_SLP",CompressionRate(hmat_SLP));
	add_stats(hmat_SLP,"Nb dense mats SLP",nb_densemats(hmat_SLP));
	add_stats(hmat_SLP,"Nb lr mats SLP",nb_lrmats(hmat_SLP));
	add_stats(hmat_SLP,"Relative Frob error SLP",sqrt(squared_absolute_error(hmat_SLP, vmat_SLP))/NormFrob(vmat_SLP));

	add_stats(hmat_DLP,"MvProd_DLP (mean)",std::get<0>(mvp_stats_DLP));
	add_stats(hmat_DLP,"MvProd_DLP (max)",std::get<1>(mvp_stats_DLP));
// 	add_stats(hmat_DLP,"MvProd err",norm(ua-ub)/norm(ua));
	add_stats(hmat_DLP,"Compression_DLP",CompressionRate(hmat_DLP));
	add_stats(hmat_DLP,"Nb dense mats DLP",nb_densemats(hmat_DLP));
	add_stats(hmat_DLP,"Nb lr mats DLP",nb_lrmats(hmat_DLP));
	add_stats(hmat_DLP,"Relative Frob error DLP",sqrt(squared_absolute_error(hmat_DLP, vmat_DLP))/NormFrob(vmat_DLP));

	
	std::cout<< "################# Harmonique n°"<<p<<" ###################"<<std::endl;
	std::cout<< "## SLP"<<std::endl; 
	print_stats(hmat_SLP);
	std::cout<< "## DLP"<<std::endl;
	print_stats(hmat_DLP);
	std::cout << std::endl;

    if (rankWorld==0){



    write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_hmat");
    write_gmsh(Vol,Out_ref,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_ref");
    }

  }
  MPI_Finalize();
}

///==================================================================================////
///============================== compare_hmat_2D ===================================////
///==================================================================================////

void compare_hmat_2D  (int argc, char* argv[], std::vector<Real> harmonics, Real lc_s, Real lc_v, Real R, std::string output_name, Real Epsilon, Real Eta, int MinClusterSize, int MaxBlockSize, int verbose=0){

  ////=====================  Setting for htool  ===================////

  MPI_Init(&argc, &argv);
  int rankWorld, sizeWorld;
  MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);
  htool::SetNdofPerElt(1);
  htool::SetEpsilon(Epsilon);
  htool::SetEta(Eta);
  htool::SetMinClusterSize(MinClusterSize);
  htool::SetMaxBlockSize(MaxBlockSize);

  ////=======================  Mesh building  =====================////

  if (rankWorld==0){
    if (verbose>0){
      std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_circle(("circle_"+NbrToStr(lc_s)).c_str(),R,lc_s,verbose);
    gmsh_disc  (("disc_"+NbrToStr(lc_v)).c_str(),R*0.9,lc_v,verbose);
  }

  ////=======================  Mesh loading  ======================////

  if (verbose>0){
    std::cout<<"Chargement du maillage"<<std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
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

  gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
  potential<P1_1D,SLP_2D> SLPop(kappa,n_);
  potential<P1_1D,DLP_2D> DLPop(kappa,n_);
  MyMatrix<P1_1D,SLP_2D> vmat_SLP(SLPop,vol,nbpt,nbdof);
  MyMatrix<P1_1D,DLP_2D> vmat_DLP(DLPop,vol,nbpt,nbdof);


  std::vector<htool::Real> rt(nbpt);
  std::vector<htool::Real> rs(nbdof);
  std::vector<int> tabt(nbpt);
  std::vector<int> tabs(nbdof);

  std::vector<htool::R3> xt(nbpt);
  std::vector<htool::R3> xs(nbdof);

  for (int i=0;i<nbpt;i++){
    tabt[i]=i;
    htool::R3 pt;pt[0]=node[i][0];pt[1]=node[i][1];pt[2]=node[i][2];
    xt[i]=pt;
    rt[i]=0;
  }

  for (int i=0;i<nbdof;i++){
    tabs[i]=i;
    htool::R3 pt;pt[0]=node_g[i][0];pt[1]=node_g[i][1];pt[2]=node_g[i][2];
    xs[i]=pt;
    rs[i]=0;
  }

  Real time=MPI_Wtime();
  htool::HMatrix hmat_SLP(vmat_SLP,xt,rt,tabt,xs,rs,tabs);
  htool::HMatrix hmat_DLP(vmat_DLP,xt,rt,tabt,xs,rs,tabs);
  time=MPI_Wtime()-time;
  
  if (rankWorld==0){
    std::cout << "Assembling time for h matrices : "<<time<<std::endl;
  }

  time=MPI_Wtime();
  for (int j=0; j<nbpt ;j++){
    const N1& jj = j;
    for (int k=0;k<nbelt;k++){
      const elt_1D& tk = Omega[k];
      const N2&     kk = dof[k];

      SLP (jj,kk) += SLPop(node[j],tk) ;
      DLP (jj,kk) += DLPop(node[j],tk) ;
    }

  }
  time=MPI_Wtime()-time;
  if (rankWorld==0){
    std::cout << "Assembling time for regular matrices : "<<time<<std::endl;
  }


  for (int l=0;l<harmonics.size();l++){
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
    vect<Cplx> Sh;resize(Sh,nbpt);fill(Sh,0.);
    vect<Cplx> D;resize(D,nbpt);fill(D,0.);
    vect<Cplx> Dh;resize(Dh,nbpt);fill(Dh,0.);
    vect<Real> Out;resize(Out,nbpt);fill(Out,0.);
    vect<Real> Out_ref;resize(Out_ref,nbpt);fill(Out_ref,0.);


    time=MPI_Wtime();
    mv_prod(S,SLP,TraceNeumann);
    mv_prod(D,DLP,TraceDirichlet);
    time=MPI_Wtime()-time;

    if (rankWorld==0){
      std::cout << "Dense matrice vector product: "<<time<<std::endl;
    }

    for (int i=0;i<nbpt;i++){
      Out[i]= std::abs(S[i]+D[i]);
      Out_ref[i]= std::abs(Ref[i]);
    }


    if (rankWorld==0){
      write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_dense");
      write_gmsh(Vol,Out_ref,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_ref");
    }

    std::pair <double,double > mvp_stats_SLP = MvProdMPI(Sh,hmat_SLP,TraceNeumann);
    std::pair <double,double > mvp_stats_DLP = MvProdMPI(Dh,hmat_DLP,TraceDirichlet);


    for (int i=0;i<nbpt;i++){
      Out[i]= std::abs(Sh[i]+Dh[i]);
    }

    Real err_SLP=0;
    Real norm_SLP=0;
    Real err_DLP=0;
    Real norm_DLP=0;
    for (int i=0;i<nbpt;i++){
      err_SLP+=std::pow(std::abs(S[i]-Sh[i]),2);
      norm_SLP+=std::pow(std::abs(S[i]),2);
      err_DLP+=std::pow(std::abs(D[i]-Dh[i]),2);
      norm_DLP+=std::pow(std::abs(D[i]),2);
    }
    err_SLP=std::sqrt(err_SLP/norm_SLP);
    err_DLP=std::sqrt(err_DLP/norm_DLP);

	
	
	add_stats(hmat_SLP,"MvProd_SLP (mean)",std::get<0>(mvp_stats_SLP));
	add_stats(hmat_SLP,"MvProd_SLP (max)",std::get<1>(mvp_stats_SLP));
	add_stats(hmat_SLP,"MvProd err_SLP",err_SLP);
	add_stats(hmat_SLP,"Compression_SLP",CompressionRate(hmat_SLP));
	add_stats(hmat_SLP,"Nb dense mats SLP",nb_densemats(hmat_SLP));
	add_stats(hmat_SLP,"Nb lr mats SLP",nb_lrmats(hmat_SLP));
	add_stats(hmat_SLP,"Relative Frob error SLP",sqrt(squared_absolute_error(hmat_SLP, vmat_SLP))/NormFrob(vmat_SLP));

	add_stats(hmat_DLP,"MvProd_DLP (mean)",std::get<0>(mvp_stats_DLP));
	add_stats(hmat_DLP,"MvProd_DLP (max)",std::get<1>(mvp_stats_DLP));
	add_stats(hmat_DLP,"MvProd err_DLP",err_DLP);
	add_stats(hmat_DLP,"Compression_SLP",CompressionRate(hmat_DLP));
	add_stats(hmat_DLP,"Nb dense mats DLP",nb_densemats(hmat_DLP));
	add_stats(hmat_DLP,"Nb lr mats DLP",nb_lrmats(hmat_DLP));
	add_stats(hmat_DLP,"Relative Frob error DLP",sqrt(squared_absolute_error(hmat_DLP, vmat_DLP))/NormFrob(vmat_DLP));
	
    if (rankWorld==0){
		write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_hmat");
    
	}

	if (rankWorld==0){
	std::cout<< "################# Harmonique n°"<<p<<" ###################"<<std::endl;
	std::cout<< "## SLP"<<std::endl;
	}
	print_stats(hmat_SLP);
	if (rankWorld==0){
		std::cout<< "## DLP"<<std::endl;
	}
	print_stats(hmat_DLP);
	if (rankWorld==0){
		std::cout << std::endl;
	}
  }
  MPI_Finalize();
}

///==================================================================================////
///=========================== compare_hmat_3D ====================================////
///==================================================================================////

void compare_hmat_3D  (int argc, char* argv[], std::vector<Real> harmonics, Real lc_s, Real lc_v, Real R, std::string output_name, Real Epsilon, Real Eta, int MinClusterSize, int MaxBlockSize, int verbose=0){
  ////=====================  Setting for htool  ===================////

  MPI_Init(&argc, &argv);
  int rankWorld, sizeWorld;
  MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);
  htool::SetNdofPerElt(1);
  htool::SetEpsilon(Epsilon);
  htool::SetEta(Eta);
  htool::SetMinClusterSize(MinClusterSize);
  htool::SetMaxBlockSize(MaxBlockSize);
  
  ////=======================  Mesh building  =====================////
  if (rankWorld==0){
    if (verbose>0){
      std::cout<<"Construction du maillage"<<std::endl;
    }
    gmsh_sphere(("sphere_"+NbrToStr(lc_s)).c_str(),R,lc_s,verbose);
    gmsh_ball  (("ball_"+NbrToStr(lc_v)).c_str(),R*0.9,lc_v,verbose);
  }

  ////=======================  Mesh loading  ======================////

  if (verbose>0){
    std::cout<<"Chargement du maillage"<<std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  geometry geom,vol;
  load_node_gmsh(geom,("sphere_"+NbrToStr(lc_s)).c_str());
  load_node_gmsh(vol ,("ball_"+NbrToStr(lc_v)).c_str());

  mesh_2D Omega(geom);
  mesh_3D Vol(vol);

  load_elt_gmsh(Omega,0);
  load_elt_gmsh(Vol,0);
  MPI_Barrier(MPI_COMM_WORLD);

  if (rankWorld==0){
//     gmsh_clean(("sphere_"+NbrToStr(lc_s)).c_str());
//     gmsh_clean(("ball_"+NbrToStr(lc_v)).c_str());
  }

  ////================== Calcul de la normale =====================////


  nrml_2D n_(Omega);


  ////============ Matrice pour le champs rayonne =================////

  const vect<R3>& node = get_node(vol);
  const vect<R3>& node_g = get_node(geom);
  int nbpt  = size(node);

  int nbelt = nb_elt(Omega);
  P1_2D dof(Omega);
  int nbdof = nb_dof(dof);

  Real kappa=10;

  gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
  potential<P1_2D,SLP_3D> SLPop(kappa,n_);
  potential<P1_2D,DLP_3D> DLPop(kappa,n_);
  MyMatrix<P1_2D,SLP_3D> vmat_SLP(SLPop,vol,nbpt,nbdof);
  MyMatrix<P1_2D,DLP_3D> vmat_DLP(DLPop,vol,nbpt,nbdof);


  std::vector<htool::Real> rt(nbpt);
  std::vector<htool::Real> rs(nbdof);
  std::vector<int> tabt(nbpt);
  std::vector<int> tabs(nbdof);

  std::vector<htool::R3> xt(nbpt);
  std::vector<htool::R3> xs(nbdof);

  for (int i=0;i<nbpt;i++){
    tabt[i]=i;
    htool::R3 pt;pt[0]=node[i][0];pt[1]=node[i][1];pt[2]=node[i][2];
    xt[i]=pt;
    rt[i]=0;
  }

  for (int i=0;i<nbdof;i++){
    tabs[i]=i;
    htool::R3 pt;pt[0]=node_g[i][0];pt[1]=node_g[i][1];pt[2]=node_g[i][2];
    xs[i]=pt;
    rs[i]=0;
  }


  Real time=MPI_Wtime();
  htool::HMatrix hmat_SLP(vmat_SLP,xt,rt,tabt,xs,rs,tabs);
  htool::HMatrix hmat_DLP(vmat_DLP,xt,rt,tabt,xs,rs,tabs);
  time=MPI_Wtime()-time;
  
  if (rankWorld==0){
    std::cout << "Assembling time for h matrices : "<<time<<std::endl;
  }

  time=MPI_Wtime();
  for (int j=0; j<nbpt ;j++){
    const N1& jj = j;
    for (int k=0;k<nbelt;k++){
      const elt_2D& tk = Omega[k];
      const N3&     kk = dof[k];

      SLP (jj,kk) += SLPop(node[j],tk) ;
      DLP (jj,kk) += DLPop(node[j],tk) ;
    }

  }
  time=MPI_Wtime()-time;
  
  if (rankWorld==0){
    std::cout << "Assembling time for regular matrices : "<<time<<std::endl;
  }
  
  for (int l=0;l<harmonics.size();l++){
    Real p = harmonics[l];

    ////======================= Trace solution ======================////

    vect<Cplx> TraceDirichlet;resize(TraceDirichlet,nbdof);fill(TraceDirichlet,0.);
    vect<Cplx> TraceNeumann;  resize(TraceNeumann,nbdof);  fill(TraceNeumann,0.);


    for (int j=0 ; j<nbelt ; j++){
      const elt_2D& tri = Omega[j];
      const N3&     I   = dof[j];

      R3 X0; X0[0] =  tri[0][0];X0[1]=tri[0][1]; X0[2]=tri[0][2];
      R3 X1; X1[0] =  tri[0][0];X1[1]=tri[0][1]; X1[2]=tri[0][2];
      R3 X2; X2[0] =  tri[0][0];X2[1]=tri[0][1]; X2[2]=tri[0][2];


      Real radius0 = sqrt(X0[0] * X0[0] + X0[1] * X0[1] + X0[2] * X0[2]);
      Real radius1 = sqrt(X1[0] * X1[0] + X1[1] * X1[1] + X1[2] * X1[2]);
      Real radius2 = sqrt(X2[0] * X2[0] + X2[1] * X2[1] + X2[2] * X2[2]);

      Real theta0 = std::acos(X0[2]/radius0);
      Real theta1 = std::acos(X1[2]/radius1);
      Real theta2 = std::acos(X2[2]/radius2);

      Real phi0 = std::atan2(X0[1],X0[0]);
      Real phi1 = std::atan2(X1[1],X1[0]);
      Real phi2 = std::atan2(X2[1],X2[0]);

      if (phi0<0){
        phi0 += 2 * M_PI;
      }
      if (phi1<0){
        phi1 += 2 * M_PI;
      }
      if (phi2<0){
        phi2 += 2 * M_PI;
      }

      C3 Vinc;
      Vinc[0] = boost::math::spherical_harmonic(p,p,theta0,phi0);
      Vinc[1] = boost::math::spherical_harmonic(p,p,theta1,phi1);
      Vinc[2] = boost::math::spherical_harmonic(p,p,theta2,phi2);

      TraceDirichlet [I] =  Vinc;

      Cplx derivee;
      derivee= kappa*((p/(kappa*R)) - boost::math::sph_bessel(p+1,kappa*R)/boost::math::sph_bessel(p,kappa*R));

      TraceNeumann[I] =  TraceDirichlet [I]*derivee;


    }

    ////=================== Solution de référence ===================////

    vect<Cplx> Ref;resize(Ref,nbpt);fill(Ref,0.);
    for (int j=0;j<nbpt;j++){
      R3 X; X[0] = node[j][0];X[1]= node[j][1]; X[2]= node[j][2];;

      Real r     = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
      Real theta;
      if (r<1e-10)
      theta=0; // to avoid nan value for points near the origin
      else
      theta = std::acos(X[2]/r);
      Real phi = std::atan2(X[1],X[0]);
      if (phi<0){
        phi += 2 * M_PI;
      }
      Ref[j]= boost::math::sph_bessel(p,kappa*r)/boost::math::sph_bessel(p,kappa*R)*boost::math::spherical_harmonic(p,p,theta,phi);

    }

    ////================ Construction de la solution ================////
    vect<Cplx> S;resize(S,nbpt);fill(S,0.);
    vect<Cplx> Sh;resize(Sh,nbpt);fill(Sh,0.);
    vect<Cplx> D;resize(D,nbpt);fill(D,0.);
    vect<Cplx> Dh;resize(Dh,nbpt);fill(Dh,0.);
    vect<Real> Out;resize(Out,nbpt);fill(Out,0.);
    vect<Real> Out_ref;resize(Out_ref,nbpt);fill(Out_ref,0.);


    time=MPI_Wtime();
    mv_prod(S,SLP,TraceNeumann);
    mv_prod(D,DLP,TraceDirichlet);
    time=MPI_Wtime()-time;

    if (rankWorld==0){
      std::cout << "Dense matrice vector product: "<<time<<std::endl;
    }

    for (int i=0;i<nbpt;i++){
      Out[i]= std::abs(S[i]+D[i]);
      Out_ref[i]= std::abs(Ref[i]);
    }


    if (rankWorld==0){
      write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_dense");
      write_gmsh(Vol,Out_ref,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_ref");
    }

    std::pair <double,double > mvp_stats_SLP = MvProdMPI(Sh,hmat_SLP,TraceNeumann);
    std::pair <double,double > mvp_stats_DLP = MvProdMPI(Dh,hmat_DLP,TraceDirichlet);


    for (int i=0;i<nbpt;i++){
      Out[i]= std::abs(Sh[i]+Dh[i]);
    }

    Real err_SLP=0;
    Real norm_SLP=0;
    Real err_DLP=0;
    Real norm_DLP=0;
    for (int i=0;i<nbpt;i++){
      err_SLP+=std::pow(std::abs(S[i]-Sh[i]),2);
      norm_SLP+=std::pow(std::abs(S[i]),2);
      err_DLP+=std::pow(std::abs(D[i]-Dh[i]),2);
      norm_DLP+=std::pow(std::abs(D[i]),2);
    }
    err_SLP=std::sqrt(err_SLP/norm_SLP);
    err_DLP=std::sqrt(err_DLP/norm_DLP);
	
	add_stats(hmat_SLP,"MvProd_SLP (mean)",std::get<0>(mvp_stats_SLP));
	add_stats(hmat_SLP,"MvProd_SLP (max)",std::get<1>(mvp_stats_SLP));
	add_stats(hmat_SLP,"MvProd err_SLP",err_SLP);
	add_stats(hmat_SLP,"Compression_SLP",CompressionRate(hmat_SLP));
	add_stats(hmat_SLP,"Nb dense mats SLP",nb_densemats(hmat_SLP));
	add_stats(hmat_SLP,"Nb lr mats SLP",nb_lrmats(hmat_SLP));
	add_stats(hmat_SLP,"Relative Frob error SLP",sqrt(squared_absolute_error(hmat_SLP, vmat_SLP))/NormFrob(vmat_SLP));

	add_stats(hmat_DLP,"MvProd_DLP (mean)",std::get<0>(mvp_stats_DLP));
	add_stats(hmat_DLP,"MvProd_DLP (max)",std::get<1>(mvp_stats_DLP));
	add_stats(hmat_DLP,"MvProd err_DLP",err_DLP);
	add_stats(hmat_DLP,"Compression_SLP",CompressionRate(hmat_DLP));
	add_stats(hmat_DLP,"Nb dense mats DLP",nb_densemats(hmat_DLP));
	add_stats(hmat_DLP,"Nb lr mats DLP",nb_lrmats(hmat_DLP));
	add_stats(hmat_DLP,"Relative Frob error DLP",sqrt(squared_absolute_error(hmat_DLP, vmat_DLP))/NormFrob(vmat_DLP));
	
    if (rankWorld==0){
		write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_hmat");
    
	}

	if (rankWorld==0){
	std::cout<< "################# Harmonique n°"<<p<<" ###################"<<std::endl;
	std::cout<< "## SLP"<<std::endl; 
	}
	print_stats(hmat_SLP);
	if (rankWorld==0){
	std::cout<< "## DLP"<<std::endl;
	}
	print_stats(hmat_DLP);
	if (rankWorld==0){
		std::cout << std::endl;
	}
  }
  MPI_Finalize();
}



int main (int argc, char *argv[]){
	// Check the number of parameters
	if (argc < 5) {
		// Tell the user how to run the program
		std::cerr << "Usage: " << argv[1] << " eps" << std::endl;
		std::cerr << "Usage: " << argv[2] << " eta" << std::endl;
		std::cerr << "Usage: " << argv[3] << " eta" << std::endl;
		std::cerr << "Usage: " << argv[4] << " eta" << std::endl;
		
		/* "Usage messages" are a conventional way of telling the user
		 * how to run a program if they enter the command incorrectly.
		 */
		return 1;
	}

	// Load the inputs
	Real Epsilon = StrToNbr<Real>(argv[1]);
	Real Eta = StrToNbr<Real>(argv[2]);
	int MinClusterSize = StrToNbr<int>(argv[3]);
	int MaxBlockSize = StrToNbr<int>(argv[4]);
	
	std::cout << Epsilon << " " << Eta <<" "<<MinClusterSize<<" "<<MaxBlockSize<< std::endl;

  std::vector<Real> harmonics(1);harmonics[0]=1;//harmonics[1]=2;harmonics[2]=3;
  Real lc_s = 0.01;
  Real lc_v = 0.05;
  Real R=1;

  potential_hmat_2D(argc,argv, harmonics,lc_s,lc_v,  R, "potential_2D",Epsilon,Eta,MinClusterSize,MaxBlockSize,0);
//   compare_hmat_2D(argc,argv, harmonics,lc_s,lc_v,  R, "compare_2D",Epsilon,Eta,MinClusterSize,MaxBlockSize,0);

  
  lc_s = 0.07;
  lc_v = 0.1;
  R=0.5;
//   potential_hmat_3D(argc,argv, harmonics,lc_s,lc_v,  R, "potential_3D",Epsilon,Eta,MinClusterSize,MaxBlockSize,0);
//   compare_hmat_3D(argc,argv, harmonics,lc_s,lc_v,  R, "compare_3D",Epsilon,Eta,MinClusterSize,MaxBlockSize,1);


}
