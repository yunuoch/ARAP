#ifndef ARAPDEFORMER_H
#define ARAPDEFORMER_H

#include <set>
#include <map>
#include <vector>

#include<Eigen\Dense>

#include"Mesh3D.h";

using namespace Eigen;



/*
  As-rigid-as-possible model for obj meshes, tet meshes and cubic meshes.

  1. For obj meshes, this class implements the model for fast geometric shape deformation described in the following paper:

	Olga Sorkine and Marc Alexa: As-Rigid-As-Possible Surface Modeling
	In Proc. of Eurographics Symposium on Geometry Processing (2007), Vol. 4, p. 30

  2. For tet and cubic meshes, the class implements the extension for volumetric meshes described in the same paper.
*/

class ARAPModel
{
public:
	ARAPModel(Mesh3D*tet = nullptr, const double * vtxWeights = nullptr);
	virtual ~ARAPModel();

	void updateRotations(const double * disp);

	void buildDispRHS(std::vector<double> rhs[3]);
	void buildPosRHS(std::vector<double> rhs[3]);


	const Vector3d & getRestVertex(int i) const { return restVertices[i]; }
	const std::vector<Vector3d> & getRestVertices() const { return restVertices; }

	const Matrix3d & getRotation(int i) const { return R[i]; }
	const std::vector<Matrix3d> & getRotations() const { return R; }

	void setRotations(const Matrix3d * rotations);

	// we omit the derivative of rotation Ri, so energy hessian is always the same: a weighted Laplacian matrix
	// the hessian is of size 3*#vtx
	//void buildEnergyHessian(SparseMatrix ** sparseMatrix) const;

	// compute energy and its gradient; the derivative of rotation is omitted in calculation
	// either energy or gradient can be nullptr, dim is of size 3*#vtx
	void getEnergyAndGradient(const double * disp, double * energy, double * gradient = nullptr);

	void buildPosRHSOnOneRotation(int vtxIdx, std::vector<double> rhs[3]);

	static void assembleDimension(const std::vector<double> input[3], double * output);

protected:
	void initialize();

	int n;

	std::vector<Vector3d> restVertices;
	std::vector<double> Lp[3];
	std::vector<Matrix3d> R;

public:
	MatrixXd Ltest;
	MatrixXd Wtest;
};


class Deform
{
public:

	Deform(Mesh3D*tet, int numFixedVertices_, const int * fixedVertices_, int numThreads_);
	virtual ~Deform();

	void setNumThreads(int threads);

	virtual void updateHandles(const std::map<int, Vector3d> & newHandles);

	// do one iteration of rotation and position optimization, return a rest shape (disp = 0) if no constraints
	void deformOneIter(double * disp); //disp serves as input and output
	virtual void deformOneIter(const double * dispLast, double * disp);

	// do rotation and position optimization until relative error < epsilon or maxIteration is reached
	void deform(const double * dispLast, double * disp, double epsilon, unsigned int maxIteration);

	vector<double> deformtest(const double * dispLast, double * disp, double epsilon, unsigned int maxIteration);

protected:
	// update handle displacements. return true if handle vertex indices don't match
	bool checkHandles(const std::map<int, Vector3d> & newHandles);

	void rebuildSolver(const std::map<int, Vector3d> & newHandles);

	//Solve 3 linear systems (n * n)
	//call Lagrange multiplier solver SolveLinearSystem()
	//store the position to sceneObjectDeformable (it uses displacement. instead of position)
	
	void optimizePositions(double * disp);
	vector<double> optimizePositionstest(/*double * disp*/);

	ARAPModel arapModel;

	std::map<int, Vector3d> handles;

	std::vector<int> fixedVertices;
	std::set<int> fixedVerticesSet;

	int numThreads;

	//per-edge Laplacian matrix: n * n

	int n, n3;

	//for testing purpose..
	void SolveLinearSystemEigen(double* x, double* rhs);

	vector<Vector3d> result;
};
#endif
