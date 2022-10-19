// #include <gtk/gtk.h>
#include <Eigen/Dense>
#include <complex>

#include <iostream>
#include <string.h>
#include <math.h>

#include "Line.h"
#include "Node.h"


class NRSSolver
{
    int numberOfIter;

public:
    NRSSolver();
    int solveNewtonRaphson(int numberOfNodes, std::vector<Node *> nodes, std::vector<Line> lines, int maxNumberOfIterations, double epsilon);
    int getNumberOfIter() { return numberOfIter; }
    void setNumberOfIter(int iter);
    Eigen::MatrixXcd getAdmittanceMatrix(int numberOfNodes, std::vector<Line> lines);
};