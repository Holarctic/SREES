#ifndef GSSOLVER_H
#define GSSOLVER_H
#include "Node.h"
#include "Line.h"
#include <vector>
#include <algorithm>    // std::max

using namespace std;
class GSSolver
{
	int numberOfIter_GS;

public:
	GSSolver();
	int solveGS(int numberOfNodes, vector<Node *> &nodes, vector<Line> lines, int maxNumberOfIterations, double epsilon);
	int getNumberOfIter() { return numberOfIter_GS; }
	void setNumberOfIter(int iter);
};
#endif