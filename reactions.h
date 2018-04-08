#include <armadillo>
#include <string>
#include <unordered_map>

using namespace arma;
using namespace std;

mat F_rxn_A(vector<mat> &config, mat &E, vector<double> steady_states);
mat G_rxn_A(vector<mat> &config, mat &E, vector<double> steady_states);
mat F_rxn_C(vector<mat> &config, mat &E, unordered_map<string, double> &RDParams);
mat G_rxn_C(vector<mat> &config, mat &E, unordered_map<string, double> &RDParams);

mat testEndoR(vector<mat> &config, mat &E);
mat testEndoS(vector<mat> &config, mat &E);
