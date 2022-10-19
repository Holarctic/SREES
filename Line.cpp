#include "Line.h"


Line:: Line(int n_startNode, int n_endNode, std::complex<double> n_z) {
	setLine(n_startNode, n_endNode, n_z);
}

void Line::setLine(int n_startNode, int n_endNode, std::complex<double> n_z) {
	Line::startNode = n_startNode;
	Line::endNode = n_endNode;
	Line::z = n_z;
}
