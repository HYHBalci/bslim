#include <Rcpp.h>
#include <stochtree/tree_sampler.h>
#include <algorithm>
#include <map>
#include <sstream>
using namespace StochTree;

static std::string state_key(Tree& tree, int node, const std::vector<double>& x) {
  if (tree.IsLeaf(node)) return "L";
  int k = std::upper_bound(x.begin(), x.end(), tree.Threshold(node)) - x.begin();
  return "(" + std::to_string(k) + ":" + state_key(tree, tree.LeftChild(node), x) +
    ":" + state_key(tree, tree.RightChild(node), x) + ")";
}

// Exercise the actual header-defined tree transition, linked to the recorded
// package DLL for data/tracker/leaf-likelihood classes. Leaf values are integrated
// out; the partial residual, Gaussian variance and leaf prior variance stay fixed.
extern "C" SEXP native_kernel_probe(SEXP args_) {
 BEGIN_RCPP
  Rcpp::List args(args_);
  auto x = Rcpp::as<std::vector<double>>(args["x"]);
  auto y = Rcpp::as<std::vector<double>>(args["y"]);
  int n = x.size(), draws = Rcpp::as<int>(args["draws"]), burn = Rcpp::as<int>(args["burn"]);
  double alpha = Rcpp::as<double>(args["alpha"]), beta = Rcpp::as<double>(args["beta"]);
  int minleaf = Rcpp::as<int>(args["minleaf"]), maxdepth = Rcpp::as<int>(args["maxdepth"]);
  std::vector<FeatureType> types{FeatureType::kNumeric};
  ForestDataset data;
  data.AddCovariates(x.data(), n, 1, false);
  ForestTracker tracker(data.GetCovariates(), types, 1, n);
  Tree tree; tree.Init(1);
  ColumnVector residual(y.data(), n);
  TreePrior prior(alpha, beta, minleaf, maxdepth);
  GaussianConstantLeafModel leaf(1.0);
  ForestContainer forests(1, 1, true, false);
  std::mt19937 rng(Rcpp::as<int>(args["seed"]));
  std::vector<double> weights{1.0};
  double xmin, xmax; VarSplitRange(tracker, data, 0, 0, 0, xmin, xmax);
  bool nonconstant = NodeNonConstant(data, tracker, 0, 0);
  std::vector<double> sorted_x = x; std::sort(sorted_x.begin(), sorted_x.end());
  Rcpp::CharacterVector states(draws);
  for (int i = 0; i < burn + draws; ++i) {
    MCMCSampleTreeOneIter<GaussianConstantLeafModel, GaussianConstantSuffStat>(
      &tree, tracker, forests, leaf, data, residual, prior, rng, weights, 0, 1.0);
    if (i >= burn) states[i - burn] = state_key(tree, 0, sorted_x);
  }
  return Rcpp::List::create(Rcpp::_ ["states"] = states, Rcpp::_ ["range"] =
    Rcpp::NumericVector::create(xmin, xmax), Rcpp::_ ["nonconstant"] = nonconstant);
 END_RCPP
}
