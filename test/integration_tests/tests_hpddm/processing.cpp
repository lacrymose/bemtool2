#include "tests_hpddm.hpp"

int main(int argc, char const *argv[]) {
  HPDDM::Option::get()->parse(argc, argv, 1);

  //// Harmoniques
  std::vector<Real> harmoniques;harmoniques.push_back(1);harmoniques.push_back(2);harmoniques.push_back(3);harmoniques.push_back(4);harmoniques.push_back(5);

  hpddm_gmm_test(harmoniques, 0.1, 1., 1);
  return 0;
}
