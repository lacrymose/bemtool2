#include <bemtool2/tools.h>
#include <bemtool2/hpddm_calls.h>

void hpddm_gmm_test(std::vector<Real> harmonics, Real lc, Real R, int verbose){
  ////=============================================================////
  ////=======================  Mesh building  =====================////
  ////=============================================================////
  if (verbose>0){
    std::cout<<"Construction du maillage"<<std::endl;
  }
  gmsh_circle(("circle_"+NbrToStr(lc)).c_str(),R,lc);


  ////=============================================================////
  ////=======================  Mesh loading  ======================////
  ////=============================================================////
  if (verbose>0){
    std::cout<<"Chargement du maillage"<<std::endl;
  }
  Real kappa=1.;

  geometry geom;
  load_node_gmsh(geom,("circle_"+NbrToStr(lc)).c_str());

  mesh_1D Omega(geom);
  load_elt_gmsh(Omega,0);
  gmsh_clean(("circle_"+NbrToStr(lc)).c_str());
  ////=============================================================////
  ////================== Calcul de la normale =====================////
  ////=============================================================////
  if (verbose>0){
    std::cout<<"Calcul de la normale"<<std::endl;
  }
  nrml_1D n_(Omega);

  ////=============================================================////
  ////================ Assemblage de la matrice ===================////
  ////=============================================================////
  if (verbose>0){
    std::cout<<"Assemblage operateurs integraux"<<std::endl;
  }
  int nbelt = nb_elt(Omega);
  P1_1D dof(Omega);
  int nbdof = nb_dof(dof);


  gmm_dense V(nbdof,nbdof),K(nbdof,nbdof),M(nbdof,nbdof);
  bem<P1_1D,P1_1D, SLP_2D>   Vop(kappa,n_,n_);
  bem<P1_1D,P1_1D, DLP_2D>   Kop(kappa,n_,n_);

  progress bar("assembly", nbelt*nbelt, verbose);
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

  for (int l=0;l<harmonics.size();l++){
    Real p = harmonics[l];

    ////=============================================================////
    ////================== Harmonique de Fourier ====================////
    ////=============================================================////

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

    ////=============================================================////
    ////================ Assemblage du second membre ================////
    ////=============================================================////

    vect<Cplx> F; resize(F,nbdof); fill(F,0);
    vect<Cplx> Ftemp; resize(Ftemp,nbdof); fill(Ftemp,0);
    vect<Cplx> gD; resize(gD,nbdof); fill(gD,0);

    // Boundary condition
    gD=Ep;
    //fill(gD,Cplx(0.));

    mv_prod(F,K,gD);
    mv_prod(Ftemp,M,gD);

    for(int j=0; j<nbelt; j++){
      F[j] = 0.5*Ftemp[j] - F[j];
    }

    ////=============================================================////
    ////================ Résolution système linéaire ================////
    ////=============================================================////
    if (verbose>0){
      std::cout<<"Appel du solveur"<<std::endl;
    }
    vect<Cplx> U_gmm;
    resize(U_gmm,nbdof);
    gmres_solve(V,U_gmm,F,40,verbose);

    vect<Cplx> U_hpddm;
    resize(U_hpddm,nbdof);
    HPDDMOperator_GMM A(V);
    A.setPrefix("V_");
    std::complex<double>* const rhs = &(F[0]);
    std::complex<double>* x = &(U_hpddm[0]);
    HPDDM::IterativeMethod::solve(A, rhs, x, 1);

    ////=============================================================////
    ////===================== Calcul de l'erreur ====================////
    ////=============================================================////
    vect<Cplx> Err, Err2, Norme;
    resize(Err, nbdof);
    resize(Err2,nbdof);
    resize(Norme,nbdof);
    for(int j=0; j<nbdof; j++){
      Err[j] =  U_gmm[j]-Ref[j];
    }
    mv_prod(Err2,M,Err);
    mv_prod(Norme,M,Ref);

    Cplx val=0;
    Cplx norme =0.;
    Real erreur_gmm=0;
    for(int j=0; j<nbdof; j++){
      val += Err2[j]*conj(Err[j]);
      norme  += Norme[j]*conj(Ref[j]);
    }
    erreur_gmm=abs(val/norme);

    for(int j=0; j<nbdof; j++){
      Err[j] =  U_hpddm[j]-Ref[j];
    }
    mv_prod(Err2,M,Err);
    mv_prod(Norme,M,Ref);

    val=0;
    norme =0.;
    Real erreur_hpddm=0;
    for(int j=0; j<nbdof; j++){
      val += Err2[j]*conj(Err[j]);
      norme  += Norme[j]*conj(Ref[j]);
    }
    erreur_hpddm=abs(val/norme);

    if (verbose>0){
      std::cout << "erreur gmm:\t" << erreur_gmm << std::endl;
      std::cout << "erreur hpddm:\t" << erreur_hpddm << std::endl;

    }
  }
}

int main(int argc, char const *argv[]) {
  HPDDM::Option::get()->parse(argc, argv, 1);

  //// Harmoniques
  std::vector<Real> harmoniques;harmoniques.push_back(1);harmoniques.push_back(2);harmoniques.push_back(3);harmoniques.push_back(4);harmoniques.push_back(5);

  hpddm_gmm_test(harmoniques, 0.1, 1., 1);
  return 0;
}
