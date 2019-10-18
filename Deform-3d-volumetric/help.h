#pragma once

#include <iostream>
#include <vector>
#include <Eigen\Dense>

using namespace std;
using namespace Eigen;

//#include "tetgen/tetgen-master/tetgen.h"

#include<fstream>  
#include<sstream>
#include<string.h>

#include<map>
#include<set>


vector<int> ArrangeCore(vector<int> a);

void WriteParticle(vector<Vector3d>& points,int num=0);