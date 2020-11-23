/* 22-11-20        first revision
Benoit Duchet, University of Oxford */

#include "mex.h"
#include "math.h"

/* sampling of random number according to a normal distribution */
double rand_realn(double mu, double sigma)
{
    double U1, U2, W, mult, output;
    static double X1;
          
        do
            
        {
            
            U1 = -1 + ((double) rand () / RAND_MAX) * 2;
            U2 = -1 + ((double) rand () / RAND_MAX) * 2;
            W = pow (U1, 2) + pow (U2, 2);
            
        }
        
        while (W >= 1 || W == 0);
        
        mult = sqrt ((-2 * log (W)) / W);
        X1 = U1 * mult;
        output =(mu + sigma * (double) X1);
        
	return output;
}

/*  Euler-Maruyama method */
void doEulerMaruyama(double wIE, double wEI, double wEE, double beta, double Tau,
                  double thetaE, double thetaI, double nMax, double sigma,
                  double dt, double E0, double I0, double *E, double *I)
{
  E[0] = E0;
  I[0] = I0;
  
  mwSize i;
  double Edot;
  double Idot;

  for (i = 0; i < (int)(nMax - 1.0); i++) {
    double randS = rand_realn(0,1);
    double randG = rand_realn(0,1);
    Edot = (-E[i] + 1.0 / (1.0 + exp(-beta * (((thetaE - wIE *
      I[i]) + wEE * E[i]) - 1.0)))) / Tau + sigma * randS / sqrt(dt);
    Idot = (-I[i] + 1.0 / (1.0 + exp(-beta * ((thetaI + wEI *
      E[i]) - 1.0)))) / Tau + sigma * randG / sqrt(dt);
    E[i+1] = E[i] + Edot *
      dt;
    I[i+1] = I[i] + Idot *
      dt;
  }
}

/* mex interface with Matlab */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
				double wIE;
				double wEI;
				double wEE;
                double beta;
                double Tau;
                double thetaE;
                double thetaI;
                double nMax;
                double sigma;
                double dt;
				double E0;
				double I0;
                double *E;
                double *I;

				wIE = mxGetScalar(prhs[0]);
				wEI = mxGetScalar(prhs[1]);
				wEE = mxGetScalar(prhs[2]);
				beta = mxGetScalar(prhs[3]);
				Tau = mxGetScalar(prhs[4]);
				thetaE = mxGetScalar(prhs[5]);
				thetaI = mxGetScalar(prhs[6]);
				nMax = mxGetScalar(prhs[7]);
				sigma = mxGetScalar(prhs[8]);
				dt = mxGetScalar(prhs[9]);
				E0 = mxGetScalar(prhs[10]);
				I0 = mxGetScalar(prhs[11]);

				plhs[0] = mxCreateDoubleMatrix(1,(mwSize)nMax,mxREAL);
				E = mxGetPr(plhs[0]);
				plhs[1] = mxCreateDoubleMatrix(1,(mwSize)nMax,mxREAL);
				I = mxGetPr(plhs[1]);
				
				doEulerMaruyama(wIE,wEI,wEE,beta,Tau,thetaE,thetaI,nMax,sigma,dt,E0,I0,E,I);
}

