#include<iostream>
#include <map>
#include"Deform.h"
#include "Mesh3D.h"

using namespace std;


int main()
{


	Deform* deformer;

	Mesh3D* tet = new Mesh3D;

	int fixv[] = { 3 };
	deformer = new Deform(tet,0, fixv, 1);

	static vector<double> disp;
	disp.resize(tet->num_of_vertices * 3,1);

	map<int, Vector3d> handle = { {0, Vector3d(2, 0, 0)}, { 2,Vector3d(-1,0,0) } };
	//map<int, Vector3d> handle = { {0,Vector3d(0,0,0)},{1,Vector3d(-2,0,0)} };
	//map<int, Vector3d> handle = { {1,Vector3d(1,2,3)},{2,Vector3d(1,2,3)},{3,Vector3d(1,2,3)} };
	deformer->updateHandles(handle);
	vector<double> temp=deformer->deformtest(&disp[0], &disp[0], 1e-5, 20);
	for (int i = 0; i < tet->num_of_vertices; i++)
	{
		cout <<i<<": "<< temp[3 * i] << " " << temp[3 * i + 1] << " " << temp[3 * i + 2] << endl;
	}


	system("pause");
	return 0;
}