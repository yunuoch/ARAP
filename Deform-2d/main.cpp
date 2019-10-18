#include <iostream>
#include <Eigen/Dense>
#include <set>

#include "Deform.h"
#include "common.h"

using namespace std;
using namespace zzs;

using Eigen::MatrixXd;

int main()
{
	//a test case for deform
	Deform* arap = new Deform();

	Vertex a(0, 0), b(0, 1), c(1, 0), d(1, 1),e(2,2);//undeformed vertices position
	arap->addPoint(a);
	arap->addPoint(b);
	arap->addPoint(c);
	arap->addPoint(d);
	arap->addPoint(e);

	vector<Vector3i> tri = { Vector3i(0,1,2),Vector3i(1,2,3),Vector3i(0,2,4) };//triangle list
	arap->buildMesh2(tri);
	//arap->buildMesh();

	const std::set<size_t> con = { 0,1,2};//deformed vertices index(constraints)
	arap->updateConstraints(con);

	arap->setVertex(0, Point(-1,0));//deformed vertices position
	arap->setVertex(1, Point(-1, 2));
	arap->setVertex(2, Point(-1, 0));

	arap->updateMesh(false);
	for (int i = 0; i < arap->getDeformedVerts().size(); i++)
	{
		//deformed vertices position
		cout <<i<<": "<< arap->getDeformedVerts()[i].vPosition[0]<<" "<<arap->getDeformedVerts()[i].vPosition[1] << endl;
	}
	

	

	system("pause");
	return 0;
}
