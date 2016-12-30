#include "tests_eigen.hpp"

using namespace Eigen;
using namespace std;
int main()
{
  // MatrixXd m = MatrixXd::Random(3,3);
  // m = (m + MatrixXd::Constant(3,3,1.2)) * 50;
  // cout << "m =" << endl << m << endl;
  // VectorXd v(3);
  // v << 1, 2, 3;
  // cout << "m * v =" << endl << m * v << endl;
  eigen_dense matrice(3,3);
  cout << "nbr de lignes : "<<  nb_rows(matrice)<< endl;
  cout << "nbr de colonnes : " << nb_cols(matrice) << endl;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrice(i,j)=i*j;
    }
  }
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      cout <<matrice(i,j);
    }
    cout <<endl;
  }
  N2 rows,cols;
  cout <<"Exctracted matrix : "<<endl;
  rows[0]=0;rows[1]=2;
  cols[0]=0;cols[1]=2;

  cout << matrice(rows,cols)<<endl;
  vect<Cplx> x;resize(x,10);
  for (int i = 0; i < size(x); i++) {
    x[i]=i;
  }
  // Eigen::Map<Eigen::Matrix<Cplx,Eigen::Dynamic,1> > x_map (x,size(x),1);
  // cout << x << endl <<x_map<<endl;



}
