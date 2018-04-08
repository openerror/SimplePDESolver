#include "main.h"
#include "reactions.h"

using namespace arma;
using namespace std;

int main(int argc, char** argv){
    arma_rng::set_seed(0);                  // For reproducable results
    InitParams();                           // Read simulation parameters
    unordered_map<string, double> RDParams = InitRDParams(); // Read reaction-diffusion rates

    cout << "Armadillo version: " << arma_version::as_string() << endl;
    printf("DeltaX, DeltaT: %le, %le\n", DeltaX, DeltaT);

    /* Set up initial conditions */
	vector<mat> config(2, mat(m,n)); 	// R = config[0]; S = config[1]
	mat E(m,n);                         // Steric hindrance term

	vector<double> steady_states = {Rs, Ss};
	double D = 1 - Rs - Ss; // Recurring constant, given dummy name D

	for (int i = 0; i < config.size(); ++i){
		double ss = steady_states[i];
		double spread = ss / 10.0;
        config[i] = init_lattice(m, n, ss, spread);
	}

	/* TEST CASE: initial cap distribution */
    // config[0] = zeros<mat>(m,n);
    // config[1] = zeros<mat>(m,n);
	// for (int y = 0; y < m; y++){
	//     for (int x = 0; x < n; x++){
	//         if (y > m/4 and y < 3*m/4 and x > n/4 and x < 3*n/4){
    //             config[0](y,x) = 1.0;
    //             config[1](y,x) = 1.0;
    //         }
	//     }
	// }
	/* Initial conditions done */

    /* Begin main loop */
    // Define the differential changes in R and S
    mat dR(m,n), dS(m,n), laplacianR(m,n), laplacianS(m,n);
    double Vs = RDParams["Vs"]; // Ratio of diffusion coefficients, scaffold to receptor

    for (int i = 0; i < StepLimit; ++i){
		// Compute Laplacians here, for repeated uses
		// Also compute the steric hindrance term E
		E = (1 - config[0] - config[1])/D;
		laplacian(config[0], laplacianR, DeltaX);
		laplacian(config[1], laplacianS, DeltaX);

        // Compute reactions
        dR = F_rxn_C(config, E, RDParams) * DeltaT;
        dS = G_rxn_C(config, E, RDParams) * DeltaT;

        // Compute diffusive terms
        // dR += laplacianR * DeltaT;
        // dS += 0.05 * laplacianS * DeltaT;
        dR += ((1.0-config[1]) % laplacianR + config[0] % laplacianS) * DeltaT;
        dS += Vs * ((1.0-config[0]) % laplacianS + config[1] % laplacianR) * DeltaT;

        // Update occupancies
        config[0] += dR;
        config[1] += dS;

        if (i % StepAvg == 0){
            string step = to_string(i);
            for (int k = 0; k < config.size(); ++k){
                string configType = to_string(k);
                string outFile = ("output/" + step + "-" + configType + ".txt");
                config[k].save(outFile, raw_ascii);
            }

			cout << "Step " << i << endl;
        }
    }
    /* End main loop */
}

mat init_lattice(int &m, int &n, double &steady_state, double &spread){
    /* Generate a random matrix of size m by n, centered at
    steady_state and has range "spread" */
    mat config = randu<mat>(m,n);
    config *= spread;
    config += steady_state;
    return config;
}

void laplacian(mat &A, mat &Z, double &dx){
    /* Return the Laplacian of A; assume periodic boundary conditions */
    unsigned m = A.n_rows-1;
    unsigned n = A.n_cols-1;

    double top, bottom, left, right, center;

    for (unsigned int y = 0; y <= m; ++y){
        for (unsigned int x = 0; x <= n; ++x){
            if (y == 0){
                top = A(m,x); bottom = A(1,x);
            }
            else if (y == m){
                top = A(m-1,x); bottom = A(0,x);
            }
            else {
                top = A(y+1,x); bottom = A(y-1,x);
            }

            if (x == 0){
                left = A(y,n); right = A(y,1);
            }
            else if (x == n){
                left = A(y,n-1); right = A(y,0);
            }
            else {
                left = A(y,x-1); right = A(y,x+1);
            }

            center = A(y,x);
            Z(y,x) = (top+left+bottom+right-4*center)/(dx*dx);
        }
    }
}
