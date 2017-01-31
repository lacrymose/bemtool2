#include <bemtool2/tools.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

using namespace bemtool;
///==================================================================================////
///============================= potential_elt_2D ===================================////
///==================================================================================////

void potential_elt_2D(std::vector<Real> harmonics, Real lc_s, Real lc_v, Real R, std::string output_name, int verbose=0){

  ////=======================  Mesh building  =====================////

  if (verbose>0){
    std::cout<<"Construction du maillage"<<std::endl;
  }
  gmsh_circle(("circle_"+NbrToStr(lc_s)).c_str(),R,lc_s,verbose);
  gmsh_disc  (("disc_"+NbrToStr(lc_v)).c_str(),R*0.9,lc_v,verbose);


  ////=======================  Mesh loading  ======================////

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

  gmsh_clean(("circle_"+NbrToStr(lc_s)).c_str());
  gmsh_clean(("disc_"+NbrToStr(lc_v)).c_str());


  ////================== Calcul de la normale =====================////


  nrml_1D n_(Omega);


  ////============ Matrice pour le champs rayonne =================////

  const vect<R3>& node = get_node(vol);
  int nbpt  = size(node);

  int nbelt = nb_elt(Omega);
  P1_1D dof(Omega);
  int nbdof = nb_dof(dof);

  Real kappa=10;

  gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
  potential<P1_1D,SLP_2D> SLPop(kappa,n_);
  potential<P1_1D,DLP_2D> DLPop(kappa,n_);

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

    ////======================= Trace solution ======================////

    vect<Cplx> TraceDirichlet;resize(TraceDirichlet,nbdof);fill(TraceDirichlet,0.);
    vect<Cplx> TraceNeumann;  resize(TraceNeumann,nbdof);  fill(TraceNeumann,0.);

    for (int j=0 ; j<nbelt ; j++){
      const elt_1D& seg = Omega[j];
      const N2&     I   = dof[j];

      R3 X0; X0[0] =  seg[0][0];X0[1]=seg[0][1]; X0[2]=0;
      R3 X1; X1[0] =  seg[0][0];X1[1]=seg[0][1]; X1[2]=0;
      Real theta0 = std::atan2 (X0[1],X0[0]);
      Real theta1 = std::atan2 (X1[1],X1[0]);

      C2 Vinc;
      Vinc[0]=exp( iu*p*theta0 );
      Vinc[1]=exp( iu*p*theta1 );

      TraceDirichlet [I] += 0.5*Vinc;

      Cplx derivee=kappa*((p/(kappa*R))-boost::math::cyl_bessel_j(p+1,kappa*R)/boost::math::cyl_bessel_j(p,kappa*R));
      Vinc[0]*=derivee;
      Vinc[1]*=derivee;


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


    mv_prod(S,SLP,TraceNeumann);
    mv_prod(D,DLP,TraceDirichlet);
    for (int i=0;i<nbpt;i++){
      Out[i]= std::abs(S[i]+D[i]);
      Out_ref[i]= std::abs(Ref[i]);
    }

    write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_elt");
    write_gmsh(Vol,Out_ref,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_ref");

  }
}

///==================================================================================////
///============================= potential_node_2D ==================================////
///==================================================================================////

void potential_node_2D(std::vector<Real> harmonics, Real lc_s, Real lc_v, Real R, std::string output_name, int verbose=0){

  ////=======================  Mesh building  =====================////

  if (verbose>0){
    std::cout<<"Construction du maillage"<<std::endl;
  }
  gmsh_circle(("circle_"+NbrToStr(lc_s)).c_str(),R,lc_s,verbose);
  gmsh_disc  (("disc_"+NbrToStr(lc_v)).c_str(),R*0.9,lc_v,verbose);


  ////=======================  Mesh loading  ======================////

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

  gmsh_clean(("circle_"+NbrToStr(lc_s)).c_str());
  gmsh_clean(("disc_"+NbrToStr(lc_v)).c_str());


  ////================== Calcul de la normale =====================////


  nrml_1D n_(Omega);


  ////============ Matrice pour le champs rayonne =================////

  const vect<R3>& node = get_node(vol);
  int nbpt  = size(node);

  int nbelt = nb_elt(Omega);
  P1_1D dof(Omega);
  int nbdof = nb_dof(dof);

  Real kappa=10;

  gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
  potential<P1_1D,SLP_2D> SLPop(kappa,n_);
  potential<P1_1D,DLP_2D> DLPop(kappa,n_);

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


    mv_prod(S,SLP,TraceNeumann);
    mv_prod(D,DLP,TraceDirichlet);
    for (int i=0;i<nbpt;i++){
      Out[i]= std::abs(S[i]+D[i]);
      Out_ref[i]= std::abs(Ref[i]);
    }

    write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_node");
    write_gmsh(Vol,Out_ref,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_ref");

  }
}

///==================================================================================////
///============================= potential_elt_3D ===================================////
///==================================================================================////

void potential_elt_3D(std::vector<Real> harmonics, Real lc_s, Real lc_v, Real R, std::string output_name, int verbose=0){

  ////=======================  Mesh building  =====================////

  if (verbose>0){
    std::cout<<"Construction du maillage"<<std::endl;
  }
  gmsh_sphere(("sphere_"+NbrToStr(lc_s)).c_str(),R,lc_s,verbose);
  gmsh_ball  (("ball_"+NbrToStr(lc_v)).c_str(),R*0.9,lc_v,verbose);


  ////=======================  Mesh loading  ======================////

  if (verbose>0){
    std::cout<<"Chargement du maillage"<<std::endl;
  }

  geometry geom,vol;
  load_node_gmsh(geom,("sphere_"+NbrToStr(lc_s)).c_str());
  load_node_gmsh(vol ,("ball_"+NbrToStr(lc_v)).c_str());

  mesh_2D Omega(geom);
  mesh_3D Vol(vol);

  load_elt_gmsh(Omega,0);
  load_elt_gmsh(Vol,0);

  gmsh_clean(("sphere_"+NbrToStr(lc_s)).c_str());
  gmsh_clean(("ball_"+NbrToStr(lc_v)).c_str());

  ////================== Calcul de la normale =====================////


  nrml_2D n_(Omega);


  ////============ Matrice pour le champs rayonne =================////

  const vect<R3>& node = get_node(vol);
  int nbpt  = size(node);

  int nbelt = nb_elt(Omega);
  P1_2D dof(Omega);
  int nbdof = nb_dof(dof);

  Real kappa=1;

  gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
  potential<P1_2D,SLP_3D> SLPop(kappa,n_);
  potential<P1_2D,DLP_3D> DLPop(kappa,n_);

  progress bar("assembly", nbpt*nbelt,verbose);
  for (int j=0; j<nbpt ;j++){
    const N1& jj = j;
    for (int k=0;k<nbelt;k++, bar++){
      const elt_2D& tk = Omega[k];
      const N3&     kk = dof[k];

      SLP (jj,kk) += SLPop(node[j],tk) ;
      DLP (jj,kk) += DLPop(node[j],tk) ;
    }

  }
  bar.end();

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
    vect<Cplx> D;resize(D,nbpt);fill(D,0.);
    vect<Real> Out;resize(Out,nbpt);fill(Out,0.);
    vect<Real> Out_ref;resize(Out_ref,nbpt);fill(Out_ref,0.);


    mv_prod(S,SLP,TraceNeumann);
    mv_prod(D,DLP,TraceDirichlet);
    for (int i=0;i<nbpt;i++){
      Out[i]= std::abs(S[i]+D[i]);
      Out_ref[i]= std::abs(Ref[i]);
    }

    write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_elt");
    write_gmsh(Vol,Out_ref,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_ref");

  }
}

///==================================================================================////
///============================= potential_node_3D ==================================////
///==================================================================================////

void potential_node_3D(std::vector<Real> harmonics, Real lc_s, Real lc_v, Real R, std::string output_name, int verbose=0){

  ////=======================  Mesh building  =====================////

  if (verbose>0){
    std::cout<<"Construction du maillage"<<std::endl;
  }
  gmsh_sphere(("sphere_"+NbrToStr(lc_s)).c_str(),R,lc_s,verbose);
  gmsh_ball  (("ball_"+NbrToStr(lc_v)).c_str(),R*0.9,lc_v,verbose);


  ////=======================  Mesh loading  ======================////

  if (verbose>0){
    std::cout<<"Chargement du maillage"<<std::endl;
  }

  geometry geom,vol;
  load_node_gmsh(geom,("sphere_"+NbrToStr(lc_s)).c_str());
  load_node_gmsh(vol ,("ball_"+NbrToStr(lc_v)).c_str());

  mesh_2D Omega(geom);
  mesh_3D Vol(vol);

  load_elt_gmsh(Omega,0);
  load_elt_gmsh(Vol,0);

  gmsh_clean(("sphere_"+NbrToStr(lc_s)).c_str());
  gmsh_clean(("ball_"+NbrToStr(lc_v)).c_str());

  ////================== Calcul de la normale =====================////


  nrml_2D n_(Omega);


  ////============ Matrice pour le champs rayonne =================////

  const vect<R3>& node = get_node(vol);
  int nbpt  = size(node);

  int nbelt = nb_elt(Omega);
  P1_2D dof(Omega);
  int nbdof = nb_dof(dof);

  Real kappa=1;

  gmm_dense SLP(nbpt,nbdof),DLP(nbpt,nbdof);
  potential<P1_2D,SLP_3D> SLPop(kappa,n_);
  potential<P1_2D,DLP_3D> DLPop(kappa,n_);

  progress bar("assembly", nbpt*nbelt,verbose);
  for (int j=0; j<nbpt ;j++){
    const N1& jj = j;
    for (int k=0;k<nbdof;k++, bar++){
      const N1&     kk = k;

      SLP (jj,kk) += SLPop(node[j],k) ;
      DLP (jj,kk) += DLPop(node[j],k) ;
    }

  }
  bar.end();

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
    vect<Cplx> D;resize(D,nbpt);fill(D,0.);
    vect<Real> Out;resize(Out,nbpt);fill(Out,0.);
    vect<Real> Out_ref;resize(Out_ref,nbpt);fill(Out_ref,0.);


    mv_prod(S,SLP,TraceNeumann);
    mv_prod(D,DLP,TraceDirichlet);
    for (int i=0;i<nbpt;i++){
      Out[i]= std::abs(S[i]+D[i]);
      Out_ref[i]= std::abs(Ref[i]);
    }

    write_gmsh(Vol,Out,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_node");
    write_gmsh(Vol,Out_ref,output_name+"_"+NbrToStr<Real>(p)+"_"+NbrToStr(lc_s)+"_"+NbrToStr(lc_v)+"_ref");

  }
}



int main(int argc, char* argv[]){
	std::vector<Real> harmonics(1);harmonics[0]=1;harmonics[1]=2;harmonics[2]=3;
	Real lc_s = 0.01;
	Real lc_v = 0.05;
	Real R=1;


  	potential_elt_2D(harmonics, lc_s,lc_v,  R, "potential_2D",0);
  	potential_node_2D(harmonics,lc_s,lc_v,  R, "potential_2D",0);
  
	lc_s = 0.07;
	lc_v = 0.1;
	R=0.5;
  	potential_elt_3D(harmonics,lc_s,lc_v,  R, "potential_3D",0);
  	potential_node_3D(harmonics,lc_s,lc_v,  R, "potential_3D",0);

}
