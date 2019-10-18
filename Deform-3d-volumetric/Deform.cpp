#include "Deform.h"
#include <float.h>
#include "triple.h"
#include <algorithm>
#include <iostream>

using namespace std;



// ARAP Energy: E = \sum_i \sum_j\in N(i) wij | (pi'-pj') - Ri (pi-pj) |^2,
//                  i: vtx index, j: nbr index in i's neighborhood N(i), pi: rest position of vtx i, pi': deformed position of i
//                  Ri: the current rotation around vtx i
//                  wij = wji = 0.5 * \sum_k cot(angle_ikj), the weight on edge (i,j) is the 0.5 * sum of each angle /_ikj which
//                  opposites the edge (i,j) on a triangle (i,j,k)
// partial E/ partial pi' = 4 \sum_j\in N(i) wij ( (pi'-pj') - 0.5 * (Ri+Rj)(pi-pj)
//                        = 4 \sum_j wij [ (pi'-pj') - \sum_j 0.5 * wij (Ri+Rj)(pi-pj) ]
// To set pE = 0, we have: L p = b, where
//               L_ij = wij is a laplacian matrix, b_i = \sum_j 0.5 * wij (Ri+Rj)(pi-pj)


void ARAPModel::initialize()
{
	R.resize(n, Matrix3d::Identity());
	VectorXd p0; //Lp = L * p0
	p0.resize(n);
	
	for (int i = 0; i < 3; i++)
	{
		Lp[i].resize(n);
		for (int j = 0; j < n; j++)
			p0[j] = restVertices[j][i];
		VectorXd Lptemp = Ltest * p0;
		for (int j = 0; j < n; j++)
		{
			Lp[i][j] = Lptemp[j];
		}
	}
}

double clamp(double a, double b, double c)
{
	if (a < b) return b;
	if (a > c) return c;
	return a;
}

#define M_PI 3.1415926

double getTriangleAngleRobust(const Vector3d & v0, const Vector3d & v1, const Vector3d & v2)
{
	Vector3d e1 = v1 - v0;
	Vector3d e2 = v2 - v0;
	double l2e1 = pow((e1).norm(), 2), l2e2 = pow((e2).norm(), 2);
	double alpha = 0.0;
	if (l2e1 > 0 && l2e2 > 0)
	{
		double cosAlpha = e1.dot(e2) / sqrt(l2e1 * l2e2);
		cosAlpha = clamp(cosAlpha, -1.0, 1.0);
		alpha = acos(cosAlpha);
	}
	else if (l2e1 == 0 && l2e2 == 0)
		alpha = M_PI / 3;
	else
		alpha = M_PI / 2;
	return alpha;
}

ARAPModel::ARAPModel(Mesh3D*tet, const double * vtxWeights)
{
	if (!tet)
	{
		tet = new Mesh3D();
	}


	tet->SimpleMesh();
	n = tet->num_of_vertices;

	vector<Vector3d> vertex = { Vector3d(0,0,0),Vector3d(1,0,0),Vector3d(0,1,0),Vector3d(0,0,1) };
	for (int i = 0; i < n; i++)
		restVertices.push_back(tet->vertices[i]);


	Ltest = MatrixXd::Zero(n, n);
	Wtest = MatrixXd::Zero(n, n);

	const double min_angle = acos(0.9999);
	const double max_angle = acos(-0.9999);


	vector<triple<int, int, double>> weightBuffer;


	for (int i = 0; i < tet->faces.size(); i++)
	{
		int v0 = tet->faces[i][0];
		int v1 = tet->faces[i][1];
		int v2 = tet->faces[i][2];
		Vector3d p0 = tet->vertices[v0];
		Vector3d p1 = tet->vertices[v1];
		Vector3d p2 = tet->vertices[v2];


		//compute 0.5 * cot (angle(p1p0p2))
		double w = 0.5 / tan(clamp(getTriangleAngleRobust(p0, p1, p2), min_angle, max_angle));
		// cout << "computing vtx " << v0 << " " << v1 << " " << v2 << " " << w << endl;
		//clamp
		if (w < 0) w = 0;

		// formula from the paper: wij = 0.5 (cot alpha_ij + cot beta_ij)
		// alpha_ij and beta_ij are the angles opposite of the mesh edge (i, j)

		weightBuffer.emplace_back(v1, v2, w);

		//cout << "a " << p0 << "b " << p1 << "c " << p2 << endl;
		//cout << endl << v0 << " " << v1 << " " << v2 << " " << w << endl;
	}

	{
		for (const auto & t : weightBuffer)
		{
			int v1 = t.first, v2 = t.second;
			double w = t.third;

			Wtest(v1, v2) += w;
			Wtest(v2, v1) += w;

			if (vtxWeights) // default vtxWeight is 1.0 for each vtx.
				w *= (vtxWeights[v1] + vtxWeights[v2]) / 2.0;

			Ltest(v1, v2) += (-w);
			Ltest(v2, v1) += (-w);
			Ltest(v1, v1) += w;
			Ltest(v2, v2) += w;

		}
	}

	initialize();
}

ARAPModel::~ARAPModel()
{

}

Vector3d tovec(const double a[])
{
	return Vector3d(a[0], a[1], a[2]);
}

void ARAPModel::updateRotations(const double * disp)
{
	double Si[9], U[9], Ri[9];

	for (int vtxIdx = 0; vtxIdx < n; vtxIdx++)
	{
		Vector3d pos_i_def, pos_j_def, pos_i_rest, pos_j_rest;
		pos_i_rest = restVertices[vtxIdx];
		pos_i_def = pos_i_rest + tovec(disp + vtxIdx * 3);


		memset(Si, 0, sizeof(double) * 9);
		Matrix3d cov;
		Matrix3d si = Matrix3d::Zero();


		for (int col = 0; col < n; col++)
		{
			if (col == vtxIdx) continue;
			pos_j_rest = restVertices[col];
			pos_j_def = pos_j_rest + tovec(disp + col * 3);

			double w = -1 * (Ltest(vtxIdx, col));
			Vector3d e_rest = pos_i_rest - pos_j_rest;
			Vector3d e_def = pos_i_def - pos_j_def;
			cov = w * e_rest* e_def.transpose();
			si += cov;
		}

		double tolerance = 1E-6;
		int forceRotation = 1;

		JacobiSVD<Eigen::MatrixXd>svd(si, ComputeThinU | ComputeThinV);

		Vector3d sigma = svd.singularValues();

		Matrix3d U = svd.matrixU();
		Matrix3d V = svd.matrixV();

		if (U.determinant() < 0) {
			U.col(2) *= -1;
			sigma(2) *= -1;
		}
		if (V.determinant() < 0) {
			V.col(2) *= -1;
			sigma(2) *= -1;
		}
		Matrix3d RiM = U * V.transpose();

		R[vtxIdx] = RiM.transpose();

	}
}


void ARAPModel::buildPosRHS(std::vector<double> rhs[3])
{
	for (int i = 0; i < 3; i++)
		rhs[i].assign(n, 0.0);

	for (int vtxIdx = 0; vtxIdx < n; vtxIdx++)
	{
		//int rowLength = L->GetRowLength(vtxIdx);
		Vector3d pos_i = restVertices[vtxIdx];
		// iterate all the neighbors of current vertex
		for (int entryIndex = 0; entryIndex < n; entryIndex++)
		{
			// w_ij
			double w = -1 * (Ltest(vtxIdx, entryIndex)); // L stores negative wij
			if (entryIndex == vtxIdx)
				continue;

			// R_i + R_j
			Matrix3d Rsum = R[vtxIdx] + R[entryIndex];
			Vector3d pos_j = restVertices[entryIndex];
			// p_i - p_j
			Vector3d pos_diff = pos_i - pos_j;
			// (w_ij / 2) * (R_i + R_j) * (p_i - p_j)
			Vector3d res = (w / 2) * Rsum * pos_diff;

			//seperate the dofs
			for (int d = 0; d < 3; d++)
				rhs[d][vtxIdx] += res[d];
		}
	}

}

void ARAPModel::buildDispRHS(std::vector<double> rhs[3])
{
	buildPosRHS(rhs);
	for (int d = 0; d < 3; d++)
		for (int i = 0; i < n; i++)
			rhs[d][i] -= Lp[d][i];
}

void ARAPModel::buildPosRHSOnOneRotation(int vtxIdx, std::vector<double> rhs[3])
{
}





void addtoarray(Vector3d& a, double b[])
{
	a[0] += b[0], a[1] += b[1], a[2] += b[2];
}






// E = \sum_i \sum_j\in N(i) wij | (pi'-pj') - Ri (pi-pj) |^2,
// partial E/ partial pi' = 4 \sum_j\in N(i) wij ( (pi'-pj') - 0.5 * (Ri+Rj)(pi-pj) )
void ARAPModel::getEnergyAndGradient(const double * disp, double * energy, double * gradient)
{
	if (gradient)
		memset(gradient, 0, n * 3);
	updateRotations(disp);

	if (energy)
		*energy = 0.;

	for (int vtxIdx = 0; vtxIdx < n; vtxIdx++)
	{
		//int rowLength = L->GetRowLength(vtxIdx);
		Vector3d resPi = restVertices[vtxIdx];
		Vector3d defPi = resPi + Vector3d(disp + 3 * vtxIdx);
		// iterate all the neighbors of current vertex
		for (int j = 0; j < n; j++)
		{
			// w_ij
			double w = -1 * (Ltest(vtxIdx, j)); // L stores negative wij
			if (j == vtxIdx)
				continue;
			Vector3d resPj = restVertices[j];
			Vector3d defPj = resPj + Vector3d(disp + 3 * j);

			// p_i - p_j
			Vector3d resDif = resPi - resPj;
			Vector3d defDif = defPi - defPj;

			// (R_i) * (p_i - p_j)
			Vector3d rotResDif = R[vtxIdx] * resDif;

			if (energy)
				*energy += w * pow((defDif - rotResDif).norm(), 2);////////////

			if (gradient)
			{
				Vector3d g = defDif - 0.5 * rotResDif;
				g *= 4 * w;
				//g.addToArray(gradient + 3 * vtxIdx);
				addtoarray(g, gradient + 3 * vtxIdx);
			}
		}
	}
}

void ARAPModel::setRotations(const Matrix3d * rotations)
{
	for (int i = 0; i < n; i++)
		R[i] = rotations[i];
	//  memcpy(R.data(), rotations, sizeof(Matrix3d) * n);
}

void ARAPModel::assembleDimension(const std::vector<double> input[3], double * output)
{
	int n = input[0].size();
	for (int i = 0; i < n; i++)
	{
		output[3 * i + 0] = input[0][i];
		output[3 * i + 1] = input[1][i];
		output[3 * i + 2] = input[2][i];
	}
}

Deform::Deform(Mesh3D*tet, int numFixedVertices_, const int * fixedVertices_, int numThreads_) :
	arapModel(tet), numThreads(numThreads_)
{
	n = tet->num_of_vertices;
	n3 = 3 * n;

	fixedVertices.resize(numFixedVertices_);
	//memcpy(fixedVertices.data(), fixedVertices_, sizeof(int) * numFixedVertices_);
	fixedVertices = vector<int>();
	fixedVerticesSet.insert(fixedVertices.begin(), fixedVertices.end());
}

Deform::~Deform()
{
}

bool Deform::checkHandles(const std::map<int, Vector3d> & newHandles)
{
	bool rebuild = false;
	if (newHandles.size() != handles.size())
		rebuild = true;
	else
	{
		// checking handle vertex indices and updating handle displacements
		map<int, Vector3d>::iterator oldIt = handles.begin();
		map<int, Vector3d>::const_iterator newIt = newHandles.begin();
		for (; newIt != newHandles.end(); oldIt++, newIt++)
			if (oldIt->first != newIt->first)
			{
				rebuild = true;
				break;
			}
			else
			{
				oldIt->second = newIt->second;
				//cout << "move handle " << oldIt->first  << " to " << oldIt->second << endl;
			}
	}
	return rebuild;
}

void Deform::updateHandles(const std::map<int, Vector3d> & newHandles)
{
	bool rebuild = checkHandles(newHandles);
	if (rebuild)
		rebuildSolver(newHandles);
}

void Deform::rebuildSolver(const std::map<int, Vector3d> & newHandles)
{
	//	for (int i = 0; i < 3; i++)
	//		delete (solver[i]);
	handles = newHandles;

	// no constraints, so return
	if (newHandles.size() == 0)
	{
		//		for (int i = 0; i < 3; i++)
		//			solver[i] = nullptr;
		return;
	}

	////////////////////cout << "Build ARAPDeformer solver with " << handles.size() << " handles" << endl;

	// form constraints

	set<int> constrainedVertices(fixedVertices.begin(), fixedVertices.end());
	for (map<int, Vector3d>::const_iterator it = newHandles.begin(); it != newHandles.end(); it++)
	{
		constrainedVertices.insert(it->first);
		assert(it->first < n);
	}
	vector<int> constrainedVerticesVec(constrainedVertices.begin(), constrainedVertices.end());
	for (size_t i = 0; i < constrainedVerticesVec.size(); i++)
		assert(constrainedVerticesVec[i] >= 0 && constrainedVerticesVec[i] < n);

	for (int i = 0; i < 3; i++)
	{
		//reconstruct the solver and do Cholesky decomposition
		int updatable = 0, addDirichlet = 1;

	}

}

void Deform::deform(const double * dispLast, double * disp, double epsilon, unsigned int maxIteration)
{
	if (handles.size() == 0)
	{
		memset(disp, 0, sizeof(double) * 3 * n);
		return;
	}

	double ep2 = epsilon * epsilon;


	vector<double> buffer(n3);
	memcpy(buffer.data(), dispLast, sizeof(double) * n3);

	double error2 = 0;
	unsigned int iter = 0;
	double dispLastLen2 = 0;
	do
	{
		dispLastLen2 = 0;
		for (int i = 0; i < n3; i++)
			dispLastLen2 += buffer[i] * buffer[i];

		deformOneIter(&buffer[0], disp);
		error2 = 0;
		for (int i = 0; i < n3; i++)
		{
			error2 += (buffer[i] - disp[i]) * (buffer[i] - disp[i]);
		}
		//////////////////////cout << sqrt(error2 / dispLastLen2) << " ";
		iter++;
		if (iter >= maxIteration)
			break;
		memcpy(buffer.data(), disp, sizeof(double) * n3);
	} while (error2 / dispLastLen2 > ep2);
	////////////////////cout << iter << " " << endl;
}

vector<double> Deform::deformtest(const double * dispLast, double * disp, double epsilon, unsigned int maxIteration)
{
	if (handles.size() == 0)
	{
		memset(disp, 0, sizeof(double) * 3 * n);
		//return;
	}

	double ep2 = epsilon * epsilon;


	vector<double> buffer(n3);//last deform config
	vector<double> buffer2(n3);//temp deform config
	memcpy(buffer.data(), dispLast, sizeof(double) * n3);

	double error2 = 0;
	unsigned int iter = 0;
	double dispLastLen2 = 0;
	do
	{
		dispLastLen2 = 0;
		for (int i = 0; i < n3; i++)
			dispLastLen2 += buffer[i] * buffer[i];

		//deformOneIter(&buffer[0], disp);
		arapModel.updateRotations(&buffer[0]);
		vector<double> temp=optimizePositionstest(/*disp*/);
		error2 = 0;
		for (int i = 0; i < n3; i++)
		{
			error2 += (buffer[i] - temp[i]) * (buffer[i] - temp[i]);
		}
		iter++;
		if (iter >= maxIteration)
			break;
		//memcpy(buffer.data(), disp, sizeof(double) * n3);		
		buffer = temp;
		buffer2 = temp;
	} while (error2 / dispLastLen2 > ep2);

	return buffer2;
}

void Deform::deformOneIter(const double * dispLast, double * disp)
{
	if (handles.size() == 0) //if there is no constraints, set to rest configuration
	{
		memset(disp, 0, sizeof(double) * 3 * n);
		return;
	}

	arapModel.updateRotations(dispLast);
	optimizePositions(disp);

}

void Deform::deformOneIter(double * disp)
{
	deformOneIter(disp, disp);
}

vector<double> Deform::optimizePositionstest(/*double * disp*/)
{
	vector<double> rhs[3], x[3];
	for (int i = 0; i < 3; i++)
	{
		rhs[i].resize(n);
		x[i].resize(n);
	}

	arapModel.buildDispRHS(rhs);

	for (std::map<int, Vector3d>::iterator it = handles.begin(); it != handles.end(); it++)
	{
		Vector3d temp = it->second;

		int vtx = it->first;
		for (int d = 0; d < 3; d++)
			x[d][vtx] = temp[d];
	}

	//arapModel.getL()->Print();

	// 3 solves
	for (int i = 0; i < 3; i++)
		SolveLinearSystemEigen(&x[i][0], &rhs[i][0]);

	vector<double> out;
	for (int k = 0; k < n; k++)
		for (int i = 0; i < 3; i++)
		{
			out.push_back(x[i][k]);
		}
	//arapModel.assembleDimension(x, disp);

	//for (int i = 0; i < n; i++)
	//{
	//	//out.push_back(Vector3d(x[0][i], x[1][i], x[2][i]));
	//}
	return out;
}

void Deform::optimizePositions(double * disp)
{
	vector<double> rhs[3], x[3];
	for (int i = 0; i < 3; i++)
	{
		rhs[i].resize(n);
		x[i].resize(n);
	}

	arapModel.buildDispRHS(rhs);

	for (std::map<int, Vector3d>::iterator it = handles.begin(); it != handles.end(); it++)
	{
		Vector3d disp = it->second;

		int vtx = it->first;
		for (int d = 0; d < 3; d++)
			x[d][vtx] = disp[d];
	}

	//arapModel.getL()->Print();

	// 3 solves
	for (int i = 0; i < 3; i++)
		SolveLinearSystemEigen(&x[i][0], &rhs[i][0]);

}

void Deform::SolveLinearSystemEigen(double* x, double* rhs)
{
	vector<int> fix_id;
	for (auto it = handles.begin(); it != handles.end(); it++)
	{
		fix_id.push_back(it->first);
	}

	int num = n;
	int num_fix = fix_id.size();
	int num_free = num - num_fix;



	MatrixXd LL, LL_cons;
	MatrixXd ADirichlet;
	VectorXd x_cons, rhs_cons;

	LL.resize(num, num);
	LL = arapModel.Ltest;
	LL_cons.resize(num - num_fix, num - num_fix);


	rhs_cons.resize(num_free);
	ADirichlet.resize(num_free, num_fix);



	vector<int> old_id;
	old_id.clear();
	for (int l = 0; l < fix_id[0]; l++)
	{
		old_id.push_back(l);
	}
	for (int k = 0; k < num_fix - 1; k++)
	{
		for (int l = fix_id[k] + 1; l < fix_id[k + 1]; l++)
		{
			//cout << "aaaa" << l << endl;
			old_id.push_back(l);
		}
	}
	for (int l = fix_id[num_fix - 1] + 1; l < num; l++)
	{
		old_id.push_back(l);
	}
	for (int k = 0; k < old_id.size(); k++)
	{
		rhs_cons[k] = rhs[old_id[k]];
		//cout << "aaaa" << old_id[k] << endl;
	}
	for (int k = 0; k < old_id.size(); k++)
		for (int l = 0; l < old_id.size(); l++)
		{
			LL_cons(k, l) = LL(old_id[k], old_id[l]);
		}
	for (int k = 0; k < old_id.size(); k++)
		for (int l = 0; l < fix_id.size(); l++)
		{
			//cout << old_id[k] << fix_id[l] << endl;
			ADirichlet(k, l) = LL(old_id[k], fix_id[l]);
		}

	VectorXd DirichletBuffer, DirichletFreeBuffer;
	DirichletBuffer.resize(num_fix);
	for (int i = 0; i < num_fix; i++)
		DirichletBuffer[i] = x[fix_id[i]];
	rhs_cons -= ADirichlet * DirichletBuffer;

	//solve linear system...
	//x_cons = LL_cons.inverse()*rhs_cons;
	x_cons = LL_cons.ldlt().solve(rhs_cons);

	for (int i = 0; i < num_fix; i++)
	{
		x[fix_id[i]] = DirichletBuffer[i];
	}
	for (int i = 0; i < num_free; i++)
	{
		x[old_id[i]] = x_cons[i];
	}
	//cout << endl;
}

void Deform::setNumThreads(int threads)
{
}
