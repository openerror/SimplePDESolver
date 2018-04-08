#include "reactions.h"

mat F_rxn_A(std::vector<mat> &config, mat &E, std::vector<double> steady_states){
    // double b = 0.1;
    // //double Rs = steady_states[0], Ss = steady_states[1];
    // //return b * (Rs/Ss * E * config[1] - config[0]);
    return (E % config[1] - config[0]);
}

mat G_rxn_A(std::vector<mat> &config, mat &E, std::vector<double> steady_states){
    double mu = 0.7, beta = 7;
    double Ss = steady_states[1];

    return (-beta + (beta-mu)*E) % config[1] + mu/Ss * (E % arma::pow(config[1],2));
}

mat F_rxn_C(std::vector<mat> &config, mat &E, unordered_map<string, double> &RDParams){
    double m1, m2, beta, mu, Rs, Ss;
    m1 = RDParams["m1"]; m2 = RDParams["m2"];
    beta = RDParams["beta"]; mu = RDParams["mu"];
    Rs = RDParams["Rs"]; Ss = RDParams["Ss"];

    // Made certain simplifcations, by assuming Rs == Ss, and b == 1.0
    return (-config[0] + m1*Rs*E  - (m1+m2)*E%config[0] + E%config[1] + m2*E%config[0]%config[1]/Rs);
}

mat G_rxn_C(std::vector<mat> &config, mat &E, unordered_map<string, double> &RDParams){
    double m1, m2, beta, mu, Rs, Ss;
    m1 = RDParams["m1"]; m2 = RDParams["m2"];
    beta = RDParams["beta"]; mu = RDParams["mu"];
    Rs = RDParams["Rs"]; Ss = RDParams["Ss"];

    return beta*(E*Ss-config[1]) + mu*(arma::pow(config[1], 2)/Ss - config[1]) % E;
}

mat testEndoR(std::vector<mat> &config, mat &E){
    double rate = 0.5;
    return (-rate * config[0]);
}

mat testEndoS(std::vector<mat> &config, mat &E){
    double rate = 0.1;
    return (-rate * config[1]);
}
