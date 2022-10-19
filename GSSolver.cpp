#include "GSSolver.h"
#include <Eigen/Dense>
#include <complex>

using namespace std;

GSSolver::GSSolver() : numberOfIter_GS(0) {}


int GSSolver::solveGS(int numberOfNodes, vector<Node *> &nodes, vector<Line> lines, int maxNumberOfIterations, double epsilon)
{

    Eigen::MatrixXcd y(numberOfNodes+1, numberOfNodes+1);
    for (int i = 0; i < numberOfNodes; i++)
    {
        for (int j = 0; j < numberOfNodes; j++)
        {
            y(i, j) = 0;
        }
    }
    for (int n = 0; n < lines.size(); n++)
    {
        if (lines[n].getStartNode() == 0)
        {
            y(lines[n].getEndNode() - 1, lines[n].getEndNode() - 1) += 1. / lines[n].getZ();
        }
        else if (lines[n].getEndNode() == 0)
        {
            y(lines[n].getStartNode() - 1, lines[n].getStartNode() - 1) += 1. / lines[n].getZ();
        }
        else
        {
            y(lines[n].getStartNode() - 1, lines[n].getStartNode() - 1) += 1. / lines[n].getZ();
            y(lines[n].getEndNode() - 1, lines[n].getEndNode() - 1) += 1. / lines[n].getZ();
            y(lines[n].getStartNode() - 1, lines[n].getEndNode() - 1) -= 1. / lines[n].getZ();
            y(lines[n].getEndNode() - 1, lines[n].getStartNode() - 1) -= 1. / lines[n].getZ();
        }
    }
    bool isDone = false;
    double maxDV = 0;
    while (!isDone && numberOfIter_GS < maxNumberOfIterations)
    {
        maxDV = 0;

        for (int j = 0; j < numberOfNodes; j++)
        {
            std::complex<double> complexV(0, 0);
            std::complex<double> complexVn(0, 0);
            std::complex<double> sum(0, 0);
            std::complex<double> S_conj(0, 0);
            double DV = 0;
            std::complex<double> new_complexV(0, 0);
            double absDV = 0;
            switch (nodes[j]->getType())
            {
            case 0:       // slack
                continue; // poznato je V i d
            case 1:       // pq

                complexV = std::complex<double>(nodes[j]->getV() * cos(nodes[j]->getD()), nodes[j]->getV() * sin(nodes[j]->getD()));

                sum = (0, 0);
                for (int n = 0; n < numberOfNodes; n++)
                {
                    if (n != j)
                    {
                        complexVn = std::complex<double>(nodes[n]->getV() * cos(nodes[n]->getD()), nodes[n]->getV() * sin(nodes[n]->getD()));
                        sum += y.coeff(j, n) * complexVn;
                    }
                }
                S_conj = std::complex<double>(nodes[j]->getP(), (-1) * nodes[j]->getQ());
                new_complexV = 1. / y.coeff(j, j) * (S_conj / std::conj(complexV) - sum);
                DV = std::abs(new_complexV) - std::abs(complexV);
                absDV = abs(DV);
                nodes[j]->setV(std::abs(new_complexV));
                nodes[j]->setD(std::arg(new_complexV));
                if (absDV > maxDV)
                {
                    maxDV = absDV;
                }
                continue;
            case 2: // pv
                std::complex<double> new_complexV_corr(0, 0);
                complexV = std::complex<double>(nodes[j]->getV() * cos(nodes[j]->getD()), nodes[j]->getV() * sin(nodes[j]->getD()));
                double new_Q = nodes[j]->getQ();
                sum = (0, 0);
                for (int n = 0; n < numberOfNodes; n++)
                {
                    if (n != j)
                    {
                        complexVn = std::complex<double>(nodes[n]->getV() * cos(nodes[n]->getD()), nodes[n]->getV() * sin(nodes[n]->getD()));
                        sum += y.coeff(j, n) * complexVn;
                    }
                }

                new_Q = -imag(conj(complexV) * (sum + y.coeff(j, j) * complexV));
                
                S_conj = std::complex<double>(nodes[j]->getP(), -new_Q);
                
                new_complexV = 1. / y.coeff(j, j) * (S_conj / std::conj(complexV) - sum);
                new_complexV_corr = nodes[j]->getV() * new_complexV / std::abs(new_complexV);
                
                DV = std::abs(new_complexV_corr) - std::abs(complexV);
                absDV = std::abs(DV);
                
                nodes[j]->setD(std::arg(new_complexV_corr));
                nodes[j]->setQ(new_Q);
                
                if (absDV > maxDV)
                {
                    maxDV = absDV;
                }
            }
        }

        if (maxDV < epsilon)
        {
            isDone = true;
            // računamo p i q za slack čvor
            for (int i = 0; i < numberOfNodes; i++)
            {
                if (nodes[i]->getType() == 0)
                {
                    std::complex<double> sum(0, 0);
                    std::complex<double> complexV(0, 0);

                    for (int k = 0; k < numberOfNodes; k++)
                    {
                        complexV = std::complex<double>(nodes[k]->getV() * cos(nodes[k]->getD()), nodes[k]->getV() * sin(nodes[k]->getD()));
                        sum += y.coeff(i, k) * complexV;
                    }
                    complexV = std::complex<double>(nodes[i]->getV() * cos(nodes[i]->getD()), nodes[i]->getV() * sin(nodes[i]->getD()));
                    std::complex<double> S(0, 0);
                    S = complexV * sum;
                    nodes[i]->setP(std::real(S));
                    nodes[i]->setQ(std::imag(S));
                }
            }
        }
        numberOfIter_GS++;
    }
    return numberOfIter_GS;
}
void GSSolver::setNumberOfIter(int iter)
{
    numberOfIter_GS = iter;
}