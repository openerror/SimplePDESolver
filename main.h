#include <armadillo>
#include <iostream>
#include <cmath>
#include <string>
#include <unordered_map>

using namespace arma;
using namespace std;

/* Input parameters (read from an input file in this order) *******************/
int m, n;                    /* m-by-n lattice on which PDEs are solved */
double DeltaX;               /* Spatial discretization (in reduced unit) */
double Rs, Ss;               /* Steady states of R and S */
unsigned int StepLimit;      /* Number of time steps to be simulated */
unsigned int StepAvg;        /* Reporting interval for statistical data */

/* Additional variables ******************************************************************/
double DeltaT;               /* Spatial discretization of lattice */

/* Functions & function prototypes ********************************************/
mat init_lattice(int &m, int &n, double &steady_state, double &spread);
void laplacian(mat &A, mat &Z, double &dx);

void InitParams() {
	// Initialize simulation parameters
	/* Reads control parameters */
	scanf("%d %d",&m,&n);
    scanf("%lf %lf",&Rs,&Ss);
	scanf("%lf",&DeltaX);
	scanf("%d",&StepLimit);
	scanf("%d",&StepAvg);

	/* Computes basic parameters */
	DeltaT = 0.2*DeltaX*DeltaX;
}

unordered_map<string, double> InitRDParams(){
	// Read reaction-diffusion rates from file
    FILE *RD = fopen("RDParams.in", "r");

    unordered_map<string, double> RDParams;
    char paramName_raw[15]; string paramName_str;
    double paramValue;

    if (RD == NULL){
        cout << "Could not open RDParams.in" << endl;
        exit(8);
    } else {
        while (!feof(RD)) {
            fscanf(RD, "%s %lf", &paramName_raw, &paramValue);
            paramName_str = string(paramName_raw);
            RDParams[paramName_str] = paramValue;
        }
    }

    fclose(RD);
    return RDParams;
}
