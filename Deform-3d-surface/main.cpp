#include<iostream>
#include"Deform.h"
#include<vector>

using namespace std;
using namespace Eigen;

int main()
{
	
	
	float p[] = { 0,0,0,1,0,0,0,1,0,0,0,1 };
	int p_num = 4;

	vector<vector<int>> tri = { {0,1,2},{0,2,3},{0,1,3},{1,2,3} };
	vector<vector<int>> adj = { {1,2,3},{0,2,3},{0,1,3},{0,1,2} };

	vector<float> pp = { 1,0,0,2,0,0,1,1,0 };
	vector<int> vv = { 0,1,2 };



	Deform* MyDeform=new Deform(p, p_num, adj, tri);



	MyDeform->set_hard_ctrs(pp, vv);

	MyDeform->do_Deform();

	float* pprime = MyDeform->get_P_Prime();
	for (int i = 0; i < p_num; i++)
	{
		cout << pprime[3 * i] << " " << pprime[3 * i + 1] << " " << pprime[3 * i + 2] << endl;
	}

	system("pause");
	return 0;
}