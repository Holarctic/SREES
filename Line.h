#ifndef LINE_H
#define LINE_H


#include <complex>

class Line {
    int startNode, endNode;
    std::complex<double> z;
public:

    Line(int n_startNode, int n_endNode, std::complex<double> n_z);

    void setLine(int n_startNode, int n_endNode, std::complex<double> n_z);

    int getStartNode() { return startNode; }
    int getEndNode() { return endNode; }
    std::complex<double> getZ() { return z; }

};
#endif
