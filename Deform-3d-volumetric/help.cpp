#include "help.h"

vector<int> ArrangeCore(vector<int> a)
{
	vector<int> b;
	b.resize(a.size());
	for (int i = 0; i < b.size(); i++)
	{
		b[i] = i;
	}

	for (int i = 0; i < a.size(); i++)
	{
		for (int j = 0; j < a.size(); j++)
			if (i != j && a[i] == a[j])
			{
				b[i] = j;
				b[j] = i;
			}
	}

	int s = 0;
	for (int j = 0; j < a.size(); j++)
	{
		if (b[j] == j) s = j;
	}

	vector<int> c;
	while (c.size() < a.size() / 2 + 1)
	{
		c.push_back(s);
		int temp;
		if (s % 2 == 0) temp = s + 1;
		else temp = s - 1;
		s = b[temp];
	}
	for (int i = 0; i < c.size(); i++)
		c[i] = a[c[i]];
	//cout << endl;
	return c;
}

void WriteParticle(vector<Vector3d>& points,int num)
{
	ofstream outfile;
	//string str = "E:/output/meshing/test/20/particle.obj";

	string buffer= "E:/output/meshing/test/20/particle";
	string over = ".obj";
	stringstream ss;
	string str;

	ss <<num;
	ss >> str;
	str += over;
	str = buffer + str;

	outfile.open(str);

	for (int i = 0; i < points.size(); i++)
	{
		outfile <<"v "<< points[i].x() << " " << points[i].y() << " " << points[i].z() << endl;
		//cout<< points[i].x() << " " << points[i].y() << " " << points[i].z() << endl;
	}
	outfile.close();
}
