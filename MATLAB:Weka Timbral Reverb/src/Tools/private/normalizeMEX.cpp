//  normalizeMEX   Perform data normalization.
//
//  Author  :  Tobias May, © 2008-2009 
//             TUe Eindhoven and Philips Research  
//             t.may@tue.nl      tobias.may@philips.com
//
//  History :  
//  v.0.1   2008/12/11
//  v.0.2   2009/10/12 cleaned up
//  ***********************************************************************

#include <math.h>
#include <cmath>
#include "mex.h"

/* Input Arguments */
#define   INPUT                prhs[0]  // input matrix
#define   NORMFLAG             prhs[1]  // normalization flag

/* Output Arguments */
#define   OUTPUT        	   plhs[0]  // output matrix
#define   NORM                 plhs[1]  // normalization statistics

void usage()
{
	mexPrintf(" normalizeMEX   Perform channel-dependent normalization of input data.\n"); 
	mexPrintf("\n"); 
	mexPrintf(" USAGE\n"); 
	mexPrintf("\t[out,scale] = normalizeMEX(input,normFlag)\n"); 
	mexPrintf("\n"); 
	mexPrintf(" INPUT ARGUMENTS\n");
	mexPrintf("\t   input    : data matrix arranged as [nSamples x nChannels] \n");
	mexPrintf("\t   normFlag : flag specifying normalization method \n");
	mexPrintf("\t              0 - no normalization \n");
	mexPrintf("\t              1 - mean normalization \n");
	mexPrintf("\t              2 - variance normalization \n");
	mexPrintf("\t              3 - mean and variance normalization (default) \n");
	mexPrintf("\n");
	mexPrintf(" OUTPUT ARGUMENTS\n");
	mexPrintf("\t        out : normalized data [nSamples x nChannels]\n");
	mexPrintf("\t       norm : mean and/or variance normalization factors per \n");
    mexPrintf("\t              channel [1|2 x nChannels]\n");
	mexPrintf("\n");
	mexPrintf("\t Author  :  Tobias May, © 2008-2009 \n");
	mexPrintf("\t            TUe Eindhoven and Philips Research   \n");
	mexPrintf("\t            t.may@tue.nl      tobias.may@philips.com \n");
	mexPrintf("\n"); 
}    


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *input, *output, *norm, mean, var;
	int    ii, jj, nSamples, nChannels, normFlag;

	// Check for proper number of arguments
	if (nrhs < 1){
		usage(); mexErrMsgTxt("Not enough input arguments.");
	}

	if (nrhs > 2){
		usage(); mexErrMsgTxt("Too many input arguments.");
	}

	if (nlhs > 2){
		usage(); mexErrMsgTxt("Too many output arguments.");
	}

	// Check dimensions of input data
	nSamples  = (int)mxGetM(INPUT);
	nChannels = (int)mxGetN(INPUT);

	// Assign pointer to the feature space 
	input = mxGetPr(INPUT);

	// Create a matrix for the return argument 
	OUTPUT = mxCreateDoubleMatrix(nSamples, nChannels, mxREAL);
	output = mxGetPr(OUTPUT);

	// Initialize normalization flag
	if (nrhs < 2){
		// Set default normalization flag (zero mean and unit variance)
		normFlag = 3;
	}
	else{
		// Get normalization flag
		normFlag = (int)mxGetScalar(NORMFLAG);

		// Check for proper normalization flag
		if (normFlag > 3 || normFlag < 0){
			mexErrMsgTxt("Invalid normalization flag. \"normFlag\" can be 0,1,2 or 3.");
		}
	}

	// Allocate memory for normalization factors
	if (nlhs > 1){
		if (normFlag == 3){
			// Return mean and variance per channel
			NORM = mxCreateDoubleMatrix(2, nChannels, mxREAL);
		}
		else{
			// Return mean or variance per channel
			NORM = mxCreateDoubleMatrix(1, nChannels, mxREAL);
		}

		// Assign pointer
		norm = mxGetPr(NORM);
	}

	// Bypass normalization
	if(normFlag == 0){
		// Loop over the number of channels
		for (ii = 0; ii < nChannels; ii++){
			// Loop over the number of samples
			for (jj = 0; jj < nSamples; jj++){
				// Copy data
				output[jj+ii*nSamples] = input[jj+ii*nSamples];
			}
		}
	}
	// Mean normalization
	else if(normFlag == 1){

		// Loop over the number of channels
		for (ii = 0; ii < nChannels; ii++){

			// Reset mean 
			mean = 0.0;

			// Loop over the number of samples
			for (jj = 0; jj < nSamples; jj++){
				mean += input[jj+ii*nSamples];
			}

			// Compute mean value
			mean /= nSamples;

			// Store normalization value
			if (nlhs > 1){
				norm[ii] = mean;
			}

			// Apply channel-dependent normalization 
			for (jj = 0; jj < nSamples; jj++){
				// Mean subtraction
				output[jj+ii*nSamples] = input[jj+ii*nSamples] - mean;
			}
		}
	}	
	// Variance and Mean & Variance normalization
	else{
		// Loop over the number of channels
		for (ii = 0; ii < nChannels; ii++){

			// Reset mean and variance estimates
			mean = 0.0;
			var  = 0.0;

			// Loop over the number of samples
			for (jj = 0; jj < nSamples; jj++){
				mean += input[jj+ii*nSamples];
				var  += (input[jj+ii*nSamples]*input[jj+ii*nSamples]);
			}

			// Compute mean and variance values
			mean /= nSamples;
			var   = var/nSamples - mean * mean;

			// Store normalization value
			if (nlhs > 1){            
				if (normFlag == 3){
					norm[ii+ii]   = mean;				
					norm[ii+1+ii] = var;
				}
				else{
					norm[ii] = var;
				}
			}

			// Apply channel-dependent normalization 
			for (jj = 0; jj < nSamples; jj++){

				// Mean subtraction
				if (normFlag == 3){
					output[jj+ii*nSamples] = input[jj+ii*nSamples] - mean;
				}
				else{
					output[jj+ii*nSamples] = input[jj+ii*nSamples];
				}

				// Variance normalization
				output[jj+ii*nSamples] /= sqrt(var);
			}
		}
	}
}
