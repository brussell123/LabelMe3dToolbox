// partsGibbsSampler.cpp

#include <time.h>
#include <string.h>
#include <math.h>
#include "mex.h"

int MultiRand(double *p,int dim) {
  double r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)));
  double tot = p[0];
  for(int i = 0; i < dim-1; i++) {
    if(tot>r) return i;
    tot += p[i+1];
  }
  return dim-1;
}

// Inputs: evidence,imgNdx,objNdx,N,q,rp,alpha_r,alpha_a,eta_a,eta_o
// Outputs: N,q
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
  // Initialize random number generator:
  srand(time(0));

  if(nrhs != 10) {
    mexErrMsgTxt("Error: 10 args. needed.");
    return;
  }

  const mxArray *evidence = prhs[0];
  int *imgNdx = (int*)mxGetData(prhs[1]);
  int *objNdx = (int*)mxGetData(prhs[2]);
  int *Nin = (int*)mxGetData(prhs[3]);
  int *qin = (int*)mxGetData(prhs[4]);
  int *rp = (int*)mxGetData(prhs[5]);
  double alpha_r = (double)(*mxGetPr(prhs[6]));
  double alpha_a = (double)(*mxGetPr(prhs[7]));
  double eta_a = (double)(*mxGetPr(prhs[8]));
  double eta_o = (double)(*mxGetPr(prhs[9]));

  int Nobjects = mxGetM(prhs[3]);
  int Ndata = mxGetNumberOfElements(prhs[4]);

  // Allocate outputs:
  int dims[2];
  dims[0] = Nobjects; dims[1] = Nobjects+2;
  plhs[0] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
  int *N = (int*)mxGetData(plhs[0]);
  dims[0] = 1; dims[1] = Ndata;
  plhs[1] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
  int *q = (int*)mxGetData(plhs[1]);

  // Allocate temporary memory:
  double *mult = (double*)mxCalloc(Ndata,sizeof(double));
  int *Cobj = (int*)mxCalloc(Nobjects,sizeof(int));

  // Initialize outputs:
  for(int i = 0; i < Ndata; i++) q[i] = qin[i];
  for(int i = 0; i < Nobjects*(Nobjects+2); i++) N[i] = Nin[i];

  for(int i = 0; i < Ndata; i++) {
    // Get pointers to sample:
    int j = rp[i]-1;
    int ii = imgNdx[j]-1;
    int oo = objNdx[j]-1;

    // Get object class of sample:
    int Nobj = mxGetNumberOfElements(mxGetField(evidence,ii,"objNdx"));
    unsigned short *evObjNdx = (unsigned short*)mxGetData(mxGetField(evidence,ii,"objNdx"));
    float *relativeOverlap = (float*)mxGetData(mxGetField(evidence,ii,"relativeOverlap"));
    float *relativeArea = (float*)mxGetData(mxGetField(evidence,ii,"relativeArea"));
    int oj = (int)evObjNdx[oo]-1;

    // Decrement counts:
    if(q[j]>-1) {
      if(q[j]>0) {
	int op = (int)evObjNdx[q[j]-1]-1;
	N[oj+op*Nobjects]--;
        N[oj+(Nobjects+1)*Nobjects]--;
      }
      else N[oj+Nobjects*Nobjects]--;
    }

    // Histogram object classes in this image:
    int Uobj = 0;
    int totobj = 0;
    for(int k = 0; k < Nobj; k++) {
      int op = (int)evObjNdx[k]-1;
      if(k!=oo) {
	if(Cobj[op]==0) {
	  Uobj++;
	  totobj += N[oj+op*Nobjects];
	}
	Cobj[op]++;
      }
    }

    // Compute posterior:
    mult[0] = log(double(N[oj+Nobjects*Nobjects])+eta_a);
    for(int k = 0; k < Nobj; k++) {
      int op = (int)evObjNdx[k]-1;
      if(k==oo) mult[k+1] = log(0);
      else mult[k+1] = log(double(N[oj+(Nobjects+1)*Nobjects])+1) + log(double(N[oj+op*Nobjects])+eta_o) - log(double(totobj)+Uobj*eta_o) - log(double(Cobj[op])) - log(alpha_a)-(alpha_a-1)*log(relativeArea[k]) + log(alpha_r)+(alpha_r-1)*log(relativeOverlap[oo+Nobj*k]);
    }

    // Normalize multinomial:
    double maxval = mult[0];
    for(int k = 1; k < Nobj+1; k++) {
      if(mult[k]>maxval) maxval = mult[k];
    }
    double tot = 0;
    for(int k = 0; k < Nobj+1; k++) {
      mult[k] = exp(mult[k]-maxval);
      tot += mult[k];
    }
    for(int k = 0; k < Nobj+1; k++) mult[k] /= tot;

    // Sample from posterior:
    int s = MultiRand(mult,Nobj+1);
    
    // Update counts:
    q[j] = s; // s is in {0,...,Nobj_j}
    if(q[j]>0) {
      int op = (int)evObjNdx[q[j]-1]-1;
      N[oj+op*Nobjects]++;
      N[oj+(Nobjects+1)*Nobjects]++;
    }
    else N[oj+Nobjects*Nobjects]++;

    // Reset histogram of object classes in this image:
    for(int k = 0; k < Nobj; k++) {
      if(k!=oo) Cobj[(int)evObjNdx[k]-1]--;
    }
  }

  // Free temporary memory:
  mxFree(mult);
  mxFree(Cobj);
}
