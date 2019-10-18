#include "Mesh3D.h"



Mesh3D::Mesh3D()
{
	vertices.clear();
	tets.clear();


	//tets = {Vector4i()}

	vertices = { Vector3d(0,0,0),Vector3d(1,0,0),Vector3d(0,1,0),Vector3d(0,0,1) };
	num_of_vertices = 4;
	tets = { Vector4i(0,1,2,3) };
	num_of_tets = 1;

	/*vertices = { Vector3d(0,0,0),Vector3d(1,0,0),Vector3d(0,1,0),Vector3d(0,0,1),
		Vector3d(0.5,0,0.5),Vector3d(0.5,0.5,0),Vector3d(0,0.5,0.5) };

	num_of_vertices = 7;
	tets = { Vector4i(0,1,4,5),Vector4i(0,2,5,6),Vector4i(0,3,4,6),Vector4i(0,4,5,6) };
	num_of_tets = 4;*/

	//cout << num_of_vertices << endl;
	//cout << num_of_tets << endl;
	crck_surf.resize(num_of_vertices);
	crck_bound.clear();
}


Mesh3D::Mesh3D(string strv, string strt)
{
	vertices.clear(); tets.clear(); 

	ifstream infile;

	infile.open(strv);
	while (!infile.eof())
	{
		double px, py, pz;
		infile >> px >> py >> pz;
		vertices.push_back(Vector3d(px, py, pz));
		//cout << px << py << endl;
	}
	infile.close();
	num_of_vertices = vertices.size();


	infile.open(strt);
	while (!infile.eof())
	{
		int px, py, pz, pw;
		infile >> px >> py >> pz >> pw;
		tets.push_back(Vector4i(px, py, pz, pw));
	}
	infile.close();
	num_of_tets = tets.size();

	//cout << num_of_vertices << endl;
	//cout << num_of_tets << endl;

	crck_surf.resize(num_of_vertices);
	crck_bound.clear();

	R.resize(num_of_vertices, Matrix3d::Identity());
	F.resize(num_of_vertices, Matrix3d::Identity());

	cout << num_of_vertices << endl;
	cout << num_of_tets << endl;
}

void readTetmeshVtk(std::istream& in, vector<Vector3d> & X, vector<Vector4i>& indices)
{
	auto initial_X_size = X.size();
	auto initial_indices_size = indices.size();
	std::string line;
	Vector3d position;
	bool reading_points = false;
	bool reading_tets = false;
	size_t n_points = 0;
	size_t n_tets = 0;
	while (std::getline(in, line)) {
		std::stringstream ss(line);
		if (line.size() == (size_t)(0)) {
		}
		else if (line.substr(0, 6) == "POINTS") {
			reading_points = true;
			reading_tets = false;
			ss.ignore(128, ' '); // ignore "POINTS"
			ss >> n_points;
		}
		else if (line.substr(0, 5) == "CELLS") {
			reading_points = false;
			reading_tets = true;
			ss.ignore(128, ' '); // ignore "CELLS"
			ss >> n_tets;
		}
		else if (line.substr(0, 10) == "CELL_TYPES") {
			reading_points = false;
			reading_tets = false;
		}
		else if (reading_points) {
			for (size_t i = 0; i < 3; i++)
				ss >> position(i);

			position(0) += 5;
			position(1) += 5;
			position(2) += 5;
			X.emplace_back(position);
		}
		else if (reading_tets) {
			ss.ignore(128, ' '); // ignore "4"
			Vector4i tet;
			for (size_t i = 0; i < 4; i++)
				ss >> tet(i);
			indices.emplace_back(tet);
		}
	}
}
Mesh3D::Mesh3D(string filename, int filetype)
{
	vertices.clear(); tets.clear();
	std::ifstream fs;
	fs.open(filename);
	//ZIRAN_ASSERT(fs, "could not open ", vtk_file);
	readTetmeshVtk(fs, vertices, tets);
	num_of_vertices = vertices.size();
	num_of_tets = tets.size();
	cout << num_of_vertices << endl;
	cout << num_of_tets << endl;
	crck_surf.resize(num_of_vertices);
	crck_bound.clear();
	if (filetype == 0)
		//read vtk file...
	{
	}

	R.resize(num_of_vertices, Matrix3d::Identity());
	F.resize(num_of_vertices, Matrix3d::Identity());
}


Mesh3D::Mesh3D(vector<Vector3d> &vert, vector<Vector4i> &tet)
{
	vertices = vert;
	tets = tet;
	num_of_vertices = vert.size();
	num_of_tets = tet.size();
}

Mesh3D::~Mesh3D()
{
}

bool Mesh3D::is_survec(int vec)
{
	for (int i = 0; i < survec.size(); i++)
	{
		if (survec[i]==vec)
		{
			return true;
		}
	}
	return false;
}

int Mesh3D::FindEdge(int p, int q)
{
	auto it = edges_.find(make_pair(p, q));
	if (it != edges_.end())
	{
		auto pe = it->second.second; int id = it->second.first;
		return pe->core_id(id);
	}

	return -1;
}

int Mesh3D::FindEdgeId(int p, int q)
{
	for (int i = 0; i < EdgetList.size(); i++)
	{
		auto pe = EdgetList[i];
		if (pe->x == p && pe->y == q || pe->y == p && pe->x == q)
		{
			return i;
		}
	}
	return -1;
}

Edget* Mesh3D::FindEdgePtr(int p, int q)
{
	auto it = edges_.find(make_pair(p, q));
	if (it != edges_.end())
	{
		return it->second.second;
	}
	return nullptr;
}

int Mesh3D::FindFace(int x, int y, int z)
{
	auto it = faces_.find(make_triple(x, y, z));
	if (it != faces_.end()) return it->second.second->core_id(it->second.first);

	it = faces_.find(make_triple(x, z, y));
	if (it != faces_.end()) return it->second.second->core_id(it->second.first);

	return -1;

}



Facet* Mesh3D::FindFacePtr(int x, int y, int z)
{
	auto it = faces_.find(make_triple(x, y, z));
	if (it != faces_.end()) return it->second.second;

	it = faces_.find(make_triple(x, z, y));
	if (it != faces_.end()) return it->second.second;

	return nullptr;
	//return faces_[make_triple(x, y, z)]+faces_[make_triple(x, z, y)];*/
}


void Mesh3D::PrintParticle()
{
	for (int i = 0; i < 10; i++)
	{
		cout << "Particle " << i << ": " << vertices[i].x() << " " 
			<< vertices[i].y() << " " << vertices[i].z() << endl;
	}
}

void Mesh3D::PrintTet()
{
	for (int i = 0; i < num_of_tets; i++)
	{
		cout << "Tet " << i << ": " << tets[i].x() << " " << tets[i].y() << " " <<
			tets[i].z() << " " << tets[i].w() << " " << endl;
	}

}

void Mesh3D::PrintEdge()
{
	for (int i = 0; i < edges.size(); i++)
	{
		cout << "Edge " << i << ": " << edges[i].x() << " "
			<< edges[i].y() << " " << edges[i].z() << endl;
	}
}

void Mesh3D::PrintFace()
{
	cout << FacetList.size() << " " << faces.size() / 3 << endl;

	for (int i = 0; i < FacetList.size(); i++)
	{
		cout << "Face " << i << ": " << faces[3*i].x() << " " << faces[3*i].y() << " " <<
			faces[3*i].z() << " " << faces[3*i].w() << " "<< faces[3 * i+1].w() << " " << faces[3 * i + 2].w() << endl;
		cout << "Facet " << i << ": " << FacetList[i]->x << " " << FacetList[i]->y << " " <<
			FacetList[i]->z << " " << FacetList[i]->core_id(0) << " " << FacetList[i]->core_id(1) << " "
			<< FacetList[i]->core_id(2) << " " << endl;
	}
}

void Mesh3D::PrintCore()
{

}


void Mesh3D::Test()
{

}

Vector3d aver(Vector3d x, Vector3d y)
{
	return y + (x - y) / 2;
}

void Mesh3D::PreProcess()
//PreProcess for future convenience
//to be modified..
{
	//add particle neighbors..
	vert_link.resize(num_of_vertices);
	for (auto pe : EdgetList)
	{
		vert_link[pe->x].push_back(pe);
		vert_link[pe->y].push_back(pe);
	}


}



int Mesh3D::FindFaceOld(int x, int y, int z)
{
	for (int i = 0; i < faces.size(); i++)
	{
		if ((faces[i][0] == x && faces[i][1] == y && faces[i][2] == z) ||
			(faces[i][0] == x && faces[i][1] == z && faces[i][2] == y))
		{
			return faces[i][3];
		}
	}
	return -1;
}

void Mesh3D::SimpleMesh()
//TODO:
{
	int id = 0;
	//tet center
	for (int i = 0; i < num_of_tets; i++)
	{
		int x = tets[i].x(), y = tets[i].y(), z = tets[i].z(), w = tets[i].w();
	}

	//face center
	//order:yzw--zwx--wxy--xyz
	id = num_of_vertices + 4 * num_of_tets;
	for (int i = 0; i < num_of_tets; i++)
	{
		int x = tets[i].x(), y = tets[i].y(), z = tets[i].z(), w = tets[i].w();
		Vector3d yzw = (vertices[y] + vertices[z] + vertices[w]) / 3,
			zwx = (vertices[z] + vertices[w] + vertices[x]) / 3,
			wxy = (vertices[w] + vertices[x] + vertices[y]) / 3,
			xyz = (vertices[x] + vertices[y] + vertices[z]) / 3;

		if (FindFaceOld(y, z, w) == -1)
		{
			int id = vertices.size();
			faces.push_back(Vector4i(y, z, w, id));
			faces.push_back(Vector4i(z, w, y, id + 1));
			faces.push_back(Vector4i(w, y, z, id + 2));

		}
		if (FindFaceOld(z, w, x) == -1)
		{
			int id = vertices.size();
			faces.push_back(Vector4i(z, w, x, id));
			faces.push_back(Vector4i(w, x, z, id + 1));
			faces.push_back(Vector4i(x, z, w, id + 2));

		}
		if (FindFaceOld(w, x, y) == -1)
		{
			int id = vertices.size();
			faces.push_back(Vector4i(w, x, y, id));
			faces.push_back(Vector4i(x, y, w, id + 1));
			faces.push_back(Vector4i(y, w, x, id + 2));

		}
		if (FindFaceOld(x, y, z) == -1)
		{
			int id = vertices.size();
			faces.push_back(Vector4i(x, y, z, id));
			faces.push_back(Vector4i(y, z, x, id + 1));
			faces.push_back(Vector4i(z, x, y, id + 2));

		}
	}

}
