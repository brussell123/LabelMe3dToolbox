// rayTrace.cpp
// To compile:
// mex rayTrace.cpp

// Two more optimizations to improve speed (if necessary):
// 1. Filter out triangles that do not appear in the viewing frustrum.
// 2. Incorporate better splitting heuristic (e.g. SAH):
// http://www.flipcode.com/archives/Raytracing_Topics_Techniques-Part_7_Kd-Trees_and_More_Speed.shtml

#include <math.h>
#include "mex.h"

#define INF FLT_MAX
// #define INF 1.0f/0.0f

// typedef struct {
//   KDtree *left;
//   KDtree *right;
//   void *value;
// } KDtree;

// Leaf is when both left and right are NULL.
struct KDtree_t {
  struct KDtree_t *left;
  struct KDtree_t *right;
  int dim;
  float split_direction;
  int *leaves;
};

typedef struct KDtree_t KDtree;

float IntersectTriangle(float *Di,int *Uaxis,int *Vaxis,int *Waxis,float *lambda_num,float *bnu,float *bnv,float *cnu,float *cnv,float *Cu,float *Cv,float *nu,float *nv,float *vu,float *vv,int *leaves,int Nleaves,int *jmax);

float IntersectTriangle_KDtree(float *Di,int *Uaxis,int *Vaxis,int *Waxis,float *lambda_num,float *bnu,float *bnv,float *cnu,float *cnv,float *Cu,float *Cv,float *nu,float *nv,float *vu,float *vv,KDtree *kd,float min_lambda,float max_lambda,float *C,int *jmax) {
  // Base case:
  if((kd->left==NULL) && (kd->right==NULL)) {
    // Intersect with leaf node's triangles:
    return IntersectTriangle(Di,Uaxis,Vaxis,Waxis,lambda_num,bnu,bnv,cnu,cnv,Cu,Cv,nu,nv,vu,vv,kd->leaves,kd->dim,jmax);
  }

  float val_min = Di[kd->dim]*min_lambda+C[kd->dim];
  float val_max = Di[kd->dim]*max_lambda+C[kd->dim];

  // Check to see which side of plane (or both) ray lies on:
  if((val_min<=kd->split_direction) && (val_max<=kd->split_direction)) {
    // Go left:
    if(kd->left) {
      return IntersectTriangle_KDtree(Di,Uaxis,Vaxis,Waxis,lambda_num,bnu,bnv,cnu,cnv,Cu,Cv,nu,nv,vu,vv,kd->left,min_lambda,max_lambda,C,jmax);
    }
    else return -INF;
  }
  else if((val_min>kd->split_direction) && (val_max>kd->split_direction)) {
    // Go right:
    if(kd->right) {
      return IntersectTriangle_KDtree(Di,Uaxis,Vaxis,Waxis,lambda_num,bnu,bnv,cnu,cnv,Cu,Cv,nu,nv,vu,vv,kd->right,min_lambda,max_lambda,C,jmax);
    }
    else return -INF;
  }
  else {
    // Intersect ray with splitting direction plane:
    float lambda = (kd->split_direction-C[kd->dim])/Di[kd->dim];

    // Try partition that is closest to camera center:
    float maxLambda = -INF;
    if(C[kd->dim]<=kd->split_direction) {
      if(kd->left) {
	maxLambda = IntersectTriangle_KDtree(Di,Uaxis,Vaxis,Waxis,lambda_num,bnu,bnv,cnu,cnv,Cu,Cv,nu,nv,vu,vv,kd->left,min_lambda,lambda,C,jmax);
      }
    }
    else if(kd->right) {
      maxLambda = IntersectTriangle_KDtree(Di,Uaxis,Vaxis,Waxis,lambda_num,bnu,bnv,cnu,cnv,Cu,Cv,nu,nv,vu,vv,kd->right,min_lambda,lambda,C,jmax);
    }
    if((min_lambda>=maxLambda)&&(maxLambda>=lambda)) return maxLambda;

    // Otherwise, try other partition:
    if(C[kd->dim]<=kd->split_direction) {
      if(kd->right) {
	maxLambda = IntersectTriangle_KDtree(Di,Uaxis,Vaxis,Waxis,lambda_num,bnu,bnv,cnu,cnv,Cu,Cv,nu,nv,vu,vv,kd->right,lambda,max_lambda,C,jmax);
      }
    }
    else if(kd->left) {
      maxLambda = IntersectTriangle_KDtree(Di,Uaxis,Vaxis,Waxis,lambda_num,bnu,bnv,cnu,cnv,Cu,Cv,nu,nv,vu,vv,kd->left,lambda,max_lambda,C,jmax);
    }
    
    return maxLambda;
  }
}

float IntersectTriangle(float *Di,int *Uaxis,int *Vaxis,int *Waxis,float *lambda_num,float *bnu,float *bnv,float *cnu,float *cnv,float *Cu,float *Cv,float *nu,float *nv,float *vu,float *vv,int *leaves,int Nleaves,int *jmax) {
  float Du, Dv, Dw, lambda, Pu, Pv, w1, w2;
  float maxLambda = -INF;
  int j;
  for(int i = 0; i < Nleaves; i++) {
    j = leaves[i];
    Du = Di[Uaxis[j]];
    Dv = Di[Vaxis[j]];
    Dw = Di[Waxis[j]];
    
    // Get distance to triangles:
    lambda = lambda_num[j]/(Dw+nu[j]*Du+nv[j]*Dv);
//     lambda = lambda_num[j]/(Dw+normals[Uaxis[j]+j*3]*Du+normals[Vaxis[j]+j*3]*Dv);
    if(lambda<=0) {
      // Get 3D intersection of ray and plane passing through each face:
      Pu = Cu[j] + lambda*Du - vu[j];
      Pv = Cv[j] + lambda*Dv - vv[j];
//       Pu = C[Uaxis[j]] + lambda*Du - cc[Uaxis[j]+j*3];
//       Pv = C[Vaxis[j]] + lambda*Dv - cc[Vaxis[j]+j*3];
      
      // Get barycentric coordinates:
      w1 = bnu[j]*Pv+bnv[j]*Pu;
      if(w1>=0) {
	w2 = cnu[j]*Pu+cnv[j]*Pv;
	if(w2>=0) {
	  if((w1+w2)<=1) {
	    if(lambda > maxLambda) {
	      maxLambda = lambda;
	      jmax[0] = j;
	    }
	  }
	}
      }
    }
  }
  return maxLambda;
}

void PreCompute(float *vertices,int *faces,int Nfaces,float *C,int *Uaxis,int *Vaxis,int *Waxis,float *lambda_num,float *bnu,float *bnv,float *cnu,float *cnv,float *Cu,float *Cv,float *nu,float *nv,float *vu,float *vv) {
  // Pre-allocate memory:
  float *aa = (float*)mxMalloc(3*Nfaces*sizeof(float));
  float *bb = (float*)mxMalloc(3*Nfaces*sizeof(float));
  float *cc = (float*)mxMalloc(3*Nfaces*sizeof(float));
  float *normals = (float*)mxMalloc(3*Nfaces*sizeof(float));

  // Compute two triangle face directions and triangle origins:
  for(int i = 0; i < Nfaces; i++) {
    for(int j = 0; j < 3; j++) {
      aa[j+i*3] = vertices[j+(faces[0+i*3]-1)*3]-vertices[j+(faces[2+i*3]-1)*3];
    }
  }
  for(int i = 0; i < Nfaces; i++) {
    for(int j = 0; j < 3; j++) {
      bb[j+i*3] = vertices[j+(faces[1+i*3]-1)*3]-vertices[j+(faces[2+i*3]-1)*3];
    }
  }
  for(int i = 0; i < Nfaces; i++) {
    for(int j = 0; j < 3; j++) {
      cc[j+i*3] = vertices[j+(faces[2+i*3]-1)*3];
    }
  }

  // Compute normal vectors for each triangle (b CROSS a):
  for(int i = 0; i < Nfaces; i++) {
    normals[0+i*3] = aa[2+i*3]*bb[1+i*3]-aa[1+i*3]*bb[2+i*3];
    normals[1+i*3] = aa[0+i*3]*bb[2+i*3]-aa[2+i*3]*bb[0+i*3];
    normals[2+i*3] = aa[1+i*3]*bb[0+i*3]-aa[0+i*3]*bb[1+i*3];
  }


  // Pre-compute numerator for lambda (for speed):
  for(int i = 0; i < Nfaces; i++) {
    lambda_num[i] = 0;
    for(int j = 0; j < 3; j++) {
      lambda_num[i] -= normals[j+i*3]*(C[j]-cc[j+i*3]);
    }
  }

  // Get dominant axis for each triangle face:
  for(int i = 0; i < Nfaces; i++) {
    if(fabs(normals[0+i*3])>fabs(normals[1+i*3])) {
      if(fabs(normals[0+i*3])>fabs(normals[2+i*3])) Waxis[i] = 0;
      else Waxis[i] = 2;
    }
    else { 
      if(fabs(normals[1+i*3])>fabs(normals[2+i*3])) Waxis[i] = 1;
      else Waxis[i] = 2;
    }
  }

  // Get minor axes for each triangle face:
  for(int i = 0; i < Nfaces; i++) Uaxis[i] = (Waxis[i]+1)%3;
  for(int i = 0; i < Nfaces; i++) Vaxis[i] = (Waxis[i]+2)%3;

  // Get camera center coordinates along minor axes:
  for(int i = 0; i < Nfaces; i++) Cu[i] = C[Uaxis[i]];
  for(int i = 0; i < Nfaces; i++) Cv[i] = C[Vaxis[i]];

  // Get Barycentric variables:
  float *norm = (float*)mxMalloc(Nfaces*sizeof(float));
  for(int i = 0; i < Nfaces; i++) norm[i] = aa[Uaxis[i]+i*3]*bb[Vaxis[i]+i*3]-aa[Vaxis[i]+i*3]*bb[Uaxis[i]+i*3];
  for(int i = 0; i < Nfaces; i++) bnu[i] = aa[Uaxis[i]+i*3]/norm[i];
  for(int i = 0; i < Nfaces; i++) bnv[i] = -aa[Vaxis[i]+i*3]/norm[i];
  for(int i = 0; i < Nfaces; i++) cnu[i] = bb[Vaxis[i]+i*3]/norm[i];
  for(int i = 0; i < Nfaces; i++) cnv[i] = -bb[Uaxis[i]+i*3]/norm[i];
  mxFree(norm);

  // Normalize normal vectors along dominant axis:
  for(int i = 0; i < Nfaces; i++) {
    normals[Uaxis[i]+i*3] /= normals[Waxis[i]+i*3];
    normals[Vaxis[i]+i*3] /= normals[Waxis[i]+i*3];
  }

  // Normalize numerator along dominant axis:
  for(int i = 0; i < Nfaces; i++) {
    lambda_num[i] /= normals[Waxis[i]+i*3];
  }

  for(int i = 0; i < Nfaces; i++) vu[i] = cc[Uaxis[i]+i*3];
  for(int i = 0; i < Nfaces; i++) vv[i] = cc[Vaxis[i]+i*3];
  for(int i = 0; i < Nfaces; i++) nu[i] = normals[Uaxis[i]+i*3];
  for(int i = 0; i < Nfaces; i++) nv[i] = normals[Vaxis[i]+i*3];

  mxFree(aa);
  mxFree(bb);
  mxFree(cc);
  mxFree(normals);
}

void FreeKDtree(KDtree *kd) {
  if(kd==NULL) return;
  FreeKDtree(kd->left);
  FreeKDtree(kd->right);
  mxFree(kd);
}

KDtree* CreateKDtree(float *vertices,int *faces,float *bound_min,float *bound_max,int dim,int *leaves,int Nleaves,int depth,int maxDepth) {
  // Parameters:
  int maxLeaves = 20;

  // Stopping conditions:
  if(Nleaves==0) {
    return NULL;
  }

  // Allocate memory:
  KDtree *kd = (KDtree*)mxMalloc(sizeof(KDtree));

  if((depth>=maxDepth) || (Nleaves<=maxLeaves)) {
    // Set leaf node:
    kd->left = NULL;
    kd->right = NULL;
    kd->leaves = leaves;
    kd->dim = Nleaves;
    return kd;
  }

  // Step 1: Get split direction
  // FOR SPEED: Use SAH cost function
  float split_direction = (bound_min[dim]+bound_max[dim])/2;

  // Step 2: Split data:
  int *whichBin = (int*)mxMalloc(Nleaves*sizeof(int));
  int Nleft = 0;
  int Nright = 0;
  for(int i = 0; i < Nleaves; i++) {
    if(vertices[dim+(faces[0+leaves[i]*3]-1)*3]<=split_direction) {
      if(vertices[dim+(faces[1+leaves[i]*3]-1)*3]<=split_direction) {
	if(vertices[dim+(faces[2+leaves[i]*3]-1)*3]<=split_direction) {
	  whichBin[i] = 1;
	}
	else whichBin[i] = 3;
      }
      else whichBin[i] = 3;
    }
    else if(vertices[dim+(faces[1+leaves[i]*3]-1)*3]>split_direction) {
      if(vertices[dim+(faces[2+leaves[i]*3]-1)*3]>split_direction) {
	whichBin[i] = 2;
      }
      else whichBin[i] = 3;
    }
    else whichBin[i] = 3;

    switch(whichBin[i]) {
    case 1:
      Nleft++;
      break;
    case 2:
      Nright++;
      break;
    case 3:
      Nleft++;
      Nright++;
      break;
    }
  }
//   fprintf(stdout,"***depth: %d\n",depth);
//   fprintf(stdout,"split: %f\n",split_direction);
//   fprintf(stdout,"Nleaves: %d\n",Nleaves);
//   fprintf(stdout,"Nleft: %d\n",Nleft);
//   fprintf(stdout,"Nright: %d\n",Nright);
//   fflush(stdout);
  
  // Allocate children leaf indices:
  int *leavesLeft = (int*)mxMalloc(Nleft*sizeof(int));
  int *leavesRight = (int*)mxMalloc(Nright*sizeof(int));
  int jLeft = 0;
  int jRight = 0;
  for(int i = 0; i < Nleaves; i++) {
    switch(whichBin[i]) {
    case 1:
      leavesLeft[jLeft++] = leaves[i];
      break;
    case 2:
      leavesRight[jRight++] = leaves[i];
      break;
    case 3:
      leavesLeft[jLeft++] = leaves[i];
      leavesRight[jRight++] = leaves[i];
      break;
    }
  }

  mxFree(whichBin);

  // Step 3: Recurse on children:
  float bound_min_right[3] = {bound_min[0],bound_min[1],bound_min[2]};
  float bound_max_left[3] = {bound_max[0],bound_max[1],bound_max[2]};
  bound_min_right[dim] = split_direction;
  bound_max_left[dim] = split_direction;
  kd->left = CreateKDtree(vertices,faces,bound_min,bound_max_left,(dim+1)%3,leavesLeft,Nleft,depth+1,maxDepth);
  kd->right = CreateKDtree(vertices,faces,bound_min_right,bound_max,(dim+1)%3,leavesRight,Nright,depth+1,maxDepth);
  kd->split_direction = split_direction;
  kd->dim = dim;

//   mxFree(leavesLeft);
//   mxFree(leavesRight);
  mxFree(leaves);

  return kd;
}

// [X,isValid] = rayTrace(D,C,vertices,faces);
//
// Inputs:
// D - 3xM matrix of direction vectors (single)
// C - 3x1 vector for camera center (single)
// vertices - 3xN matrix of vertices (single)
// faces - 3xK matrix of triangle face indices (int32)
//
// Outputs:
// X - 3xM matrix of 3D points that intersects the mesh (single)
// isValid - 1xM vector that indicates whether 3D point is valid intersection (logic)
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
  if(nrhs != 4) {
    mexErrMsgTxt("Error: 4 args. needed.");
    return;
  }

  // Get inputs:
  float *D = (float*)mxGetData(prhs[0]);
  float *C = (float*)mxGetData(prhs[1]);
  float *vertices = (float*)mxGetData(prhs[2]);
  int *faces = (int*)mxGetData(prhs[3]);

  int Nrays = mxGetN(prhs[0]);
  int Nvertices = mxGetN(prhs[2]);
  int Nfaces = mxGetN(prhs[3]);

  fprintf(stdout,"Nfaces: %d\n",Nfaces);
  fflush(stdout);

  // Allocate outputs:
  plhs[0] = mxCreateNumericMatrix(3,Nrays,mxSINGLE_CLASS,mxREAL);
  float *X = (float*)mxGetData(plhs[0]);
  plhs[1] = mxCreateLogicalMatrix(1,Nrays);
  bool *isValid = (bool*)mxGetData(plhs[1]);
  plhs[2] = mxCreateNumericMatrix(1,Nrays,mxINT32_CLASS,mxREAL);
  int *faceNdx = (int*)mxGetData(plhs[2]);

  // Allocate temporary memory:
  float *lambda_num = (float*)mxMalloc(Nfaces*sizeof(float));
  int *Uaxis = (int*)mxMalloc(Nfaces*sizeof(int));
  int *Vaxis = (int*)mxMalloc(Nfaces*sizeof(int));
  int *Waxis = (int*)mxMalloc(Nfaces*sizeof(int));
  float *Cu = (float*)mxMalloc(Nfaces*sizeof(float));
  float *Cv = (float*)mxMalloc(Nfaces*sizeof(float));
  float *bnu = (float*)mxMalloc(Nfaces*sizeof(float));
  float *bnv = (float*)mxMalloc(Nfaces*sizeof(float));
  float *cnu = (float*)mxMalloc(Nfaces*sizeof(float));
  float *cnv = (float*)mxMalloc(Nfaces*sizeof(float));
  float *vu = (float*)mxMalloc(Nfaces*sizeof(float));
  float *vv = (float*)mxMalloc(Nfaces*sizeof(float));
  float *nu = (float*)mxMalloc(Nfaces*sizeof(float));
  float *nv = (float*)mxMalloc(Nfaces*sizeof(float));

  // Pre-compute quantities on mesh:
  PreCompute(vertices,faces,Nfaces,C,Uaxis,Vaxis,Waxis,lambda_num,bnu,bnv,cnu,cnv,Cu,Cv,nu,nv,vu,vv);

  // Compute KD-tree:
  fprintf(stdout,"Creating KDtree...\n");
  fflush(stdout);
  int maxDepth = 20;//3;
  float bound_min[3] = {INF,INF,INF};
  float bound_max[3] = {-INF,-INF,-INF};
  int *leaves = (int*)mxMalloc(Nfaces*sizeof(int));
  for(int i = 0; i < Nfaces; i++) leaves[i] = i;
  for(int i = 0; i < Nvertices; i++) {
    for(int j = 0; j < 3; j++) {
      if(vertices[j+i*3]<bound_min[j]) bound_min[j] = vertices[j+i*3];
      if(vertices[j+i*3]>bound_max[j]) bound_max[j] = vertices[j+i*3];
    }
  }
//   for(int i = 0; i < 3; i++) fprintf(stdout,"Bounds %d: %f %f\n",i,bound_min[i],bound_max[i]);
  KDtree *root = CreateKDtree(vertices,faces,bound_min,bound_max,0,leaves,Nfaces,0,maxDepth);
//   mxFree(leaves);
  fprintf(stdout,"Finished creating KDtree...\n");
  fflush(stdout);

  int Nleaves = Nfaces;
  leaves = (int*)mxMalloc(Nleaves*sizeof(int));
  for(int i = 0; i < Nleaves; i++) leaves[i] = i;

  // Get 3D points:
  float maxLambda;
  float min_lambda,max_lambda;
//   float *Di;
  float Di[3];
  int j;
  float ll;
  int jmax;
  for(int i = 0; i < Nrays; i++) {
//     Di = D+i*3;
    // This was originally written for LH coordinates; change for RH coordinates
    Di[0] = -D[i*3];
    Di[1] = -D[1+i*3];
    Di[2] = -D[2+i*3];

    // Get bounds on ray:
//     min_lambda = -INF;
    max_lambda = -INF;
    for(j = 0; j < 3; j++) {
      ll = (bound_min[j]-C[j])/Di[j];
      if(ll<=0) {
	if((C[j]>=bound_min[j])&&(ll>max_lambda)) max_lambda = ll;
// 	if((C[j]<bound_min[j])&&(ll>min_lambda)) min_lambda = ll;
      }
      ll = (bound_max[j]-C[j])/Di[j];
      if(ll<=0) {
// 	if((C[j]>=bound_max[j])&&(ll>min_lambda)) min_lambda = ll;
	if((C[j]<bound_max[j])&&(ll>max_lambda)) max_lambda = ll;
      }
    }
    min_lambda = 0.0f;

//     fprintf(stdout,"Min/max lambda: (%f,%f)\n",min_lambda,max_lambda);
//     fflush(stdout);
    
    // Intersect triangles:
    maxLambda = IntersectTriangle_KDtree(Di,Uaxis,Vaxis,Waxis,lambda_num,bnu,bnv,cnu,cnv,Cu,Cv,nu,nv,vu,vv,root,min_lambda,max_lambda,C,&jmax);
//     maxLambda = IntersectTriangle(Di,Uaxis,Vaxis,Waxis,lambda_num,bnu,bnv,cnu,cnv,Cu,Cv,nu,nv,vu,vv,leaves,Nleaves);
    if(maxLambda != -INF) {
      isValid[i] = 1;
      faceNdx[i] = jmax+1;
      for(j = 0; j < 3; j++) X[j+i*3] = -D[j+i*3]*maxLambda+C[j];
    }
  }

  mxFree(leaves);

  // Free memory:
  mxFree(lambda_num);
  mxFree(Uaxis);
  mxFree(Vaxis);
  mxFree(Waxis);
  mxFree(Cu);
  mxFree(Cv);
  mxFree(bnu);
  mxFree(bnv);
  mxFree(cnu);
  mxFree(cnv);
  mxFree(vu);
  mxFree(vv);
  mxFree(nu);
  mxFree(nv);

  // Deallocate KDtree:
  FreeKDtree(root);
}
