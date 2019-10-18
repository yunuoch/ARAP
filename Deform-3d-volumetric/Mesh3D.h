#pragma once

#include"help.h"
#include"triple.h"

class Facet
{
public:
	Facet(int x, int y, int z, int x_id, int y_id, int z_id)
		:x(x), y(y), z(z), is_surf(true), flag(false), core_id(Vector3i(x_id, y_id, z_id)),
		is_break(0, 0, 0)
	{};
	int x, y, z;
	bool is_surf;
	bool flag;

	Vector3i core_id;
	Vector3i is_break;//initial as 0 0 0

	//some added data structure...
	//Edget* yz, zx, zy;
	Vector3i edge_id;//edge_index of y-z/z-x/x-y
};

class Edget
{
public:
	Edget(int x, int y)
		:x(x),y(y),max_strech(0),is_break(0),orig_dist(0),is_surfedge(false)
	{};
	Edget(int x, int y, double dist, int x_id, int y_id)
		:x(x), y(y), max_strech(0), is_break(0), orig_dist(dist), core_id(Vector2i(x_id, y_id)), is_surfedge(false)
	{};

	int x, y;
	
	double max_strech;
	double orig_dist;

	bool is_break;//initial as 0

	Vector2i core_id;

	bool is_surfedge;

};


class Mesh3D
{
public:
	Mesh3D();
	//Mesh3D(tetgenio & out);
	Mesh3D(string strv, string strt);//read vert&tet directly from file
	//Mesh3D(vector<Vector3d> &vert);//read vert and generate tets
	Mesh3D(vector<Vector3d> &vert, vector<Vector4i> &tet);
	Mesh3D(string filename, int filetype);
	//Mesh3D(string strv);//read directly from 0.txt..
	~Mesh3D();

	int num_of_vertices;
	int num_of_tets;

	int num_of_faces=0;
	

	vector<Vector3d> vertices;
	vector<Vector4i> tets;

	vector<Vector3i> edges;//(x,y,i)
	vector<Vector4i> faces;//(x,y,z,i)

	//map<pair<int, int>, triple<int,double,double>> edges_;
	map<pair<int, int>, pair<int,Edget*>> edges_;

	map<triple<int, int, int>, pair<int,Facet*>> faces_;
	
	vector<Facet*> FacetList;
	vector<Edget*> EdgetList;

	int num_of_edges;

	vector<int> intface, surface;
	vector<int> intvec, survec;

	vector<vector<int>> crck_surf;//crack surface stored at each particle..//id in core triple...
	vector<Vector3i> crck_bound;

	vector<Matrix3d> R, F;


	bool is_survec(int vec);

	vector<double> distance_original;
	vector<double> max_strech;
	//map<pair<int, int>, double> distance_original, max_strech;

	vector<vector<int>> core_;//e1,f1,t1,e2,f2,t2.....

	int FindEdge(int p, int q);
	int FindEdgeId(int p, int q);
	Edget* FindEdgePtr(int p, int q);
	int FindFace(int x, int y, int z);
	//vector<int> FindFaceId(int x, int y);//find face like (u,x,y) or (u,y,x)....
	Facet* FindFacePtr(int x, int y,int z);//find face like(x,y,z) or (x,z,y)

	void PrintParticle();
	void PrintTet();
	void PrintEdge();
	void PrintFace();
	void PrintCore();
	


	void Test();

	void PreProcess();


	Mesh3D copy()
	{
		Mesh3D cpy(*this);
		cpy.FacetList.clear();
		cpy.EdgetList.clear();
		for (int i = 0; i < FacetList.size(); i++)
		{
			Facet* pf = new Facet(*(FacetList[i]));
			cpy.FacetList.push_back(pf);
		}
		for (int i = 0; i < EdgetList.size(); i++)
		{
			Edget* pe = new Edget(*(EdgetList[i]));
			cpy.EdgetList.push_back(pe);
		}
		cout << cpy.FacetList.size() << FacetList.size() << endl;
		return cpy;
	};

	vector<vector<Edget*>> vert_link;


	//temporary usage for arap
	//to be modified!
	void SimpleMesh();
	int FindFaceOld(int x, int y, int z);

};

