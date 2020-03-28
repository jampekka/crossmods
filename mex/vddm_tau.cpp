#include "mex.h"

#include "../vddm.hpp"

double getParam(const mxArray* p, const char* key) {
	mxArray* fuckme = mxGetField(p, 0, key);
	if(!fuckme) {
		mexErrMsgIdAndTxt("nofield", "No value for %s", key);
		return 0.0;
	}

	return mxGetScalar(fuckme);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs < 2) {
		mexErrMsgIdAndTxt("noargs", "At least 2 arguments required");
		return;
	}
	
	const mxArray* tau = prhs[0];
	const mxArray* p = prhs[1];

	mexPrintf("%f\n", getParam(p, "std"));
	crossmods::Vddm model(
		getParam(p, "dt"),
		getParam(p, "std"),
		getParam(p, "damping"),
		getParam(p, "tau_threshold"),
		getParam(p, "pass_threshold"),
		getParam(p, "scale"),
		getParam(p, "act_threshold")
	);

	crossmods::Grid1d grid(-3.0, 3.0, 100); // TODO: Parameterize

	size_t len = mxGetN(tau);
	mxArray *pdf = mxCreateDoubleMatrix(1, len, mxREAL);
	double undecided = model.decisions(grid, mxGetPr(tau), mxGetPr(pdf), len);

	plhs[0] = pdf;
	if(nlhs >= 2) {
		plhs[1] = mxCreateDoubleScalar(undecided);
	}
	
}
