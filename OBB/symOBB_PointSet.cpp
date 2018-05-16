//--------------------------------------------------------------------------
//Developed by : Diego Gonzalez Becerra(*)
//(*)Prof.Oliver van Kaick's masters student at Carleton University. (2016-2017)
//---------------------------------------------------------------------------
// This file contains the code to compute an "Oriented Bounding Box" (OBB) for a set of points in 3D
// This is a MEX-file to be called from Matlab. The calling syntax is:
//			[mbox, boxes] = symOBB_PointSet(samples, coh, calcWithSym);
// Input:
//	samples:		set of 3D samples or points (a point cloud): [n*3] matrix
//	coh:			[m*3] matrix representing the Convex Hull of the points given by the parameter "samples"
//	calcWithSym:	scalar value indicating if "symmetry planes" will be used as part of the computation, (see Reference a)
//					if calcWithSym == 1, then, the "symmetry" will be used, otherwise Not
//---------------------------------------------------------------------------
// *****	*****
//IMPORTANT NOTE: 
//This function is only a modified version of the MEX function (file) "sym_bb" 
//originally created by Noa Fish, who kindly gave us the code to compute the OBB for a Mesh.
//This code is only an adaptation for the case of a Point Set, (a "point cloud" representation)
// *****	*****
//References:
//	(a) Fish, N., Averkiou, M., van Kaick, O., et al. "Meta-representation of shape families". (2014)
//---------------------------------------------------------------------------
#include "Wm5Vector3.h"
#include "Wm5ContSymBox3.h"
#include "Wm5Memory.h"
#include "Wm5Query.h"
#include "mex.h"

using namespace Wm5;

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]) {

	//const mxArray *mVertices;
	//-----
	//Matlab's FUNCTION THAT IS CALLED FROM THIS MEX REQUIRES 3 INPUTS: 
	//(1) The set of 3D samples or points
	//(2) The convex-hull of the set of points, computed with Matlab's function convhull (or convhulln)
	const mxArray *mSamples;	//Set of 3d points
	const mxArray *mFaces;		//Set of indices of the triangular faces of the Convex Hull (CH) of the points
	const mxArray *mUseSymm;	//Do use "symmetry planes" (1) or Not (any other value <> 1)
	//-----
	//const mxArray *mOrigFaces;

	//double *vertices, *samples, *faces, *orig_faces;
	double *samples, *faces;
	double *out_box, *out_planes;
	double useSymm;
	//int nV, nS, nF, nOF;
	int nS, nF;

	/* Check number of input and output parameters */
	//if (nrhs != 4) {
	if (nrhs != 3) {
		mexErrMsgTxt("The function 'symOBB_PointSet' requires 3 input arguments");
	}
	if (nlhs != 2) {
		mexErrMsgTxt("The function 'symOBB_PointSet' requires 2 outputs");
	}

	/* Get matlab inputs */
	mSamples = prhs[0];
	//mVertices = prhs[1];
	mFaces = prhs[1];
	//mOrigFaces = prhs[3];
	mUseSymm = prhs[2];

	/* Get dimensions of input data */
	//nV = mxGetM(mVertices);
	nS = mxGetM(mSamples);
	nF = mxGetM(mFaces);
	//nOF = mxGetM(mOrigFaces);

	/* Get pointer to input data */
	//vertices = mxGetPr(mVertices);
	samples = mxGetPr(mSamples);
	faces = mxGetPr(mFaces);
	useSymm = mxGetScalar(mUseSymm);
	//orig_faces = mxGetPr(mOrigFaces);
	
	//Vector3<double> *points, *samps, *facs, *orfacs;
	Vector3<double> *samps, *facs;

	/* Transfer the samples and faces of the CH from one structure to the other */
	samps = new1<Vector3d>(nS);
	for (int i = 0; i < nS; i++) {
		samps[i] = Vector3<double>(samples[i], samples[nS + i], samples[2 * nS + i]);
	}

	facs = new1<Vector3d>(nF);
	for (int i = 0; i < nF; i++) {
		facs[i] = Vector3<double>(faces[i] - 1, faces[nF + i] - 1, faces[2 * nF + i] - 1);
	}
		
	/* Compute box */
	double epsilon = 0.00001;
	Query::Type queryType = Query::QT_REAL;
	//Query::Type queryType = Query::QT_RATIONAL;
	SymBox3<double> bx = SymBox3<double>(nS, nF, samps, facs, epsilon, queryType, (useSymm == 1 ? true : false));
	//delete1(points);
	delete1(samps);
	delete1(facs);
	//delete1(orfacs);

	int prob = bx.getProblem();
	if (prob)
	{
		plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
		return;
	}

	/* Create output data */
	plhs[0] = mxCreateDoubleMatrix(3, 5, mxREAL);

	/* Get pointer to output data */
	out_box = mxGetPr(plhs[0]);

	Box3<double> box = (Box3<double>) bx;

	/* Transfer box to output */
	out_box[0] = box.Center[0];
	out_box[1] = box.Center[1];
	out_box[2] = box.Center[2];
	out_box[3] = box.Axis[0][0];
	out_box[4] = box.Axis[0][1];
	out_box[5] = box.Axis[0][2];
	out_box[6] = box.Axis[1][0];
	out_box[7] = box.Axis[1][1];
	out_box[8] = box.Axis[1][2];
	out_box[9] = box.Axis[2][0];
	out_box[10] = box.Axis[2][1];
	out_box[11] = box.Axis[2][2];
	out_box[12] = box.Extent[0];
	out_box[13] = box.Extent[1];
	out_box[14] = box.Extent[2];

	int numBoxes = bx.getNumBoxes();

	if (numBoxes == 0)
	{
		plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
		return;
	}

	plhs[1] = mxCreateDoubleMatrix(18, numBoxes, mxREAL);
	out_planes = mxGetPr(plhs[1]);

	//double* normals = bx.getNormals();
	//double* ppoints = bx.getPoints();
	double* extents = bx.getExtents();
	double* axes = bx.getAxes();
	double* centers = bx.getCenters();
	double* scores = bx.getScores();

	for (int i = 0; i < numBoxes; i++)
	{

		out_planes[18 * i] = centers[3 * i];
		out_planes[18 * i + 1] = centers[3 * i + 1];
		out_planes[18 * i + 2] = centers[3 * i + 2];

		out_planes[18 * i + 3] = axes[9 * i];
		out_planes[18 * i + 4] = axes[9 * i + 1];
		out_planes[18 * i + 5] = axes[9 * i + 2];
		out_planes[18 * i + 6] = axes[9 * i + 3];
		out_planes[18 * i + 7] = axes[9 * i + 4];
		out_planes[18 * i + 8] = axes[9 * i + 5];
		out_planes[18 * i + 9] = axes[9 * i + 6];
		out_planes[18 * i + 10] = axes[9 * i + 7];
		out_planes[18 * i + 11] = axes[9 * i + 8];

		out_planes[18 * i + 12] = extents[3 * i];
		out_planes[18 * i + 13] = extents[3 * i + 1];
		out_planes[18 * i + 14] = extents[3 * i + 2];

		out_planes[18 * i + 15] = scores[3 * i];
		out_planes[18 * i + 16] = scores[3 * i + 1];
		out_planes[18 * i + 17] = scores[3 * i + 2];

	}

	bx.freeFields();

}
