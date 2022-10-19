#include "NRSolver.h"   
using namespace std;

NRSSolver::NRSSolver()
{
    NRSSolver::numberOfIter = 0;
}

struct PQ
{
    vector<double> p, q;
};
struct XElement
{
    int index;
    double value;
};


struct PQ getNodePQ(Node *node, vector<Node *> nodes, int numberOfNodes, Eigen::MatrixXcd admittanceMatrix)
{
    double p = 0.0;
    double q = 0.0;
    int i = node->getNumber() - 1;
    double v_i = node->getV();
    double d_i = node->getD();
    for (int j = 0; j < numberOfNodes; j++)
    {
        complex<double> line_data = admittanceMatrix(i, j);
        double line_abs = abs(line_data);
        double line_arg = arg(line_data);
        double p_tmp = v_i * nodes[j]->getV() * line_abs * cos(line_arg - d_i + nodes[j]->getD());
        p += p_tmp;
        double q_tmp = v_i * nodes[j]->getV() * line_abs * sin(line_arg - d_i + nodes[j]->getD());
        q += q_tmp;

    }
    q *= -1;
    vector<double> vector_p{p};
    vector<double> vector_q{q};
    return PQ{vector_p, vector_q};
}

double get_j1_diag_element(Node *node, int numberOfNodes, vector<Node *> nodes, const Eigen::MatrixXcd &admittanceMatrix)
{
    double sum = 0;
    int i = node->getNumber() - 1;
    for(int j = 0;j<numberOfNodes;j++) 
    {
        if (node->getNumber() == nodes[j]->getNumber())
        {
            continue;
        }
        complex<double> line_data = admittanceMatrix(i, j);
        double product = abs(node->getV()) * abs(nodes[j]->getV()) * abs(line_data);
        double angle = arg(line_data) - node->getD() + nodes[j]->getD();
        sum += (product * sin(angle));
    }
    return sum;
}

double get_j1_off_diag_element(Node *node_i, Node *node_j, const Eigen::MatrixXcd &admittanceMatrix)
{
    complex<double> line_data = admittanceMatrix(node_i->getNumber() - 1, node_j->getNumber() - 1);
    double product = abs(node_i->getV()) * abs(node_j->getV()) * abs(line_data);
    double angle = arg(line_data) - node_i->getD() + node_j->getD();

    return -1.0 * product * sin(angle);
}

vector<vector<double>> getJ1(vector<Node *> nodes, int numberOfNodes, const Eigen::MatrixXcd &admittanceMatrix)
{
    vector<vector<double>> jacobian_1;

    for (int i = 0; i < numberOfNodes; i++)
    {
        if (nodes[i]->getType() == 0)
        {
            continue;
        }
        vector<double> tmp;
        for (int j = 0; j < numberOfNodes; j++)
        {
            if (nodes[j]->getType() == 0)
            {
                continue;
            }
            if (i == j)
            {
                tmp.push_back(get_j1_diag_element(nodes[i], numberOfNodes, nodes, admittanceMatrix));
            }
            else
            {
                tmp.push_back(get_j1_off_diag_element(nodes[i], nodes[j], admittanceMatrix));
            }
        }
        jacobian_1.push_back(tmp);
    }
    return jacobian_1;
}

double get_j2_diag_element(Node *node, int numberOfNodes, vector<Node *> nodes, const Eigen::MatrixXcd &admittanceMatrix)
{
    int i = node->getNumber() - 1;
    complex<double> node_admittance = admittanceMatrix(i, i);
    double sum = 2 * abs(node->getV()) * abs(node_admittance) * cos(arg(node_admittance));
    for (int j = 0; j < numberOfNodes; j++)
    {
        if (node->getNumber() == nodes[j]->getNumber())
        {
            continue;
        }
        complex<double> line_data = admittanceMatrix(i, j);
        double product = abs(nodes[j]->getV()) * abs(line_data);
        double angle = arg(line_data) - node->getD() + nodes[j]->getD();
        sum += (product * cos(angle));
    }
    return sum;
}

double get_j2_off_diag_element(Node *node_i, Node *node_j, const Eigen::MatrixXcd &admittanceMatrix)
{
    complex<double> line_data = admittanceMatrix(node_i->getNumber() - 1, node_j->getNumber() - 1);
    double product = abs(node_i->getV()) * abs(line_data);
    double angle = arg(line_data) - node_i->getD() + node_j->getD();

    return product * cos(angle);
}

vector<vector<double>> getJ2(vector<Node *> nodes, int numberOfNodes, const Eigen::MatrixXcd &admittanceMatrix)
{
    vector<vector<double>> jacobian_2;

    for (int i = 0; i < numberOfNodes; i++)
    {
        if (nodes[i]->getType() == 0)
        {
            continue;
        }
        vector<double> tmp;
        for (int j = 0; j < numberOfNodes; j++)
        {
            if (nodes[j]->getType() == 0)
            {
                continue;
            }
            if (nodes[j]->getType() == 2)
            {
                continue;
            }
            if (i == j)
            {
                tmp.push_back(get_j2_diag_element(nodes[i], numberOfNodes, nodes, admittanceMatrix));
            }
            else
            {
                tmp.push_back(get_j2_off_diag_element(nodes[i], nodes[j], admittanceMatrix));
            }
        }
        jacobian_2.push_back(tmp);
    }
    return jacobian_2;
}

double get_j3_diag_element(Node *node, int numberOfNodes, vector<Node *> nodes, const Eigen::MatrixXcd &admittanceMatrix)
{
    int i = node->getNumber() - 1;
    complex<double> node_admittance = admittanceMatrix(i, i);
    double sum = 0;
    for (int j = 0; j < numberOfNodes; j++)
    {
        if (node->getNumber() == nodes[j]->getNumber())
        {
            continue;
        }
        complex<double> line_data = admittanceMatrix(i, j);
        double product = abs(node->getV()) * abs(nodes[j]->getV()) * abs(line_data);
        double angle = arg(line_data) - node->getD() + nodes[j]->getD();
        sum += (product * cos(angle));
    }
    return sum;
}

double get_j3_off_diag_element(Node *node_i, Node *node_j, const Eigen::MatrixXcd &admittanceMatrix)
{
    complex<double> line_data = admittanceMatrix(node_i->getNumber() - 1, node_j->getNumber() - 1);
    double product = abs(node_i->getV()) * abs(node_j->getV()) * abs(line_data);
    double angle = arg(line_data) - node_i->getD() + node_j->getD();

    return -1.0 * product * cos(angle);
}

vector<vector<double>> getJ3(vector<Node *> nodes, int numberOfNodes, const Eigen::MatrixXcd &admittanceMatrix)
{
    vector<vector<double>> jacobian_3;

    for (int i = 0; i < numberOfNodes; i++)
    {
        if (nodes[i]->getType() == 0)
        {
            continue;
        }
        if (nodes[i]->getType() == 2)
        {
            continue;
        }
        vector<double> tmp;
        for (int j = 0; j < numberOfNodes; j++)
        {
            if (nodes[j]->getType() == 0)
            {
                continue;
            }
            if (i == j)
            {
                tmp.push_back(get_j3_diag_element(nodes[i], numberOfNodes, nodes, admittanceMatrix));
            }
            else
            {
                tmp.push_back(get_j3_off_diag_element(nodes[i], nodes[j], admittanceMatrix));
            }
        }
        jacobian_3.push_back(tmp);
    }
    return jacobian_3;
}

double get_j4_diag_element(Node *node, int numberOfNodes, vector<Node *> nodes, const Eigen::MatrixXcd &admittanceMatrix)
{
    int i = node->getNumber() - 1;
    complex<double> node_admittance = admittanceMatrix(i, i);
    double sum = 0;
    for (int j = 0; j < numberOfNodes; j++)
    {
        if (node->getNumber() == nodes[j]->getNumber())
        {
            continue;
        }
        complex<double> line_data = admittanceMatrix(i, j);
        double product = abs(nodes[j]->getV()) * abs(line_data);
        double angle = arg(line_data) - node->getD() + nodes[j]->getD();
        sum += (product * sin(angle));
    }
    double prefix = -2.0 * abs(node->getV()) * abs(node_admittance);
    prefix *= sin(arg(node_admittance));
    return prefix - sum;
}

double get_j4_off_diag_element(Node *node_i, Node *node_j, const Eigen::MatrixXcd &admittanceMatrix)
{
    complex<double> line_data = admittanceMatrix(node_i->getNumber() - 1, node_j->getNumber() - 1);
    double product = abs(node_i->getV()) * abs(line_data);
    double angle = arg(line_data) - node_i->getD() + node_j->getD();

    return -1.0 * product * sin(angle);
}

vector<vector<double>> getJ4(vector<Node *> nodes, int numberOfNodes, const Eigen::MatrixXcd &admittanceMatrix)
{
    vector<vector<double>> jacobian_4;

    for (int i = 0; i < numberOfNodes; i++)
    {
        if (nodes[i]->getType() == 0)
        {
            continue;
        }
        if (nodes[i]->getType() == 2)
        {
            continue;
        }
        vector<double> tmp;
        for (int j = 0; j < numberOfNodes; j++)
        {
            if (nodes[j]->getType() == 2)
            {
                continue;
            }
            if (nodes[j]->getType() == 0)
            {
                continue;
            }
            if (i == j)
            {
                tmp.push_back(get_j4_diag_element(nodes[i], numberOfNodes, nodes, admittanceMatrix));
            }
            else
            {
                tmp.push_back(get_j4_off_diag_element(nodes[i], nodes[j], admittanceMatrix));
            }
        }
        jacobian_4.push_back(tmp);
    }
    return jacobian_4;
}

// start SO copy
double getDeterminant(const vector<vector<double>> vect)
{
    if (vect.size() != vect[0].size())
    {
        throw runtime_error("Matrix is not quadratic");
    }
    int dimension = vect.size();

    if (dimension == 0)
    {
        return 1;
    }

    if (dimension == 1)
    {
        return vect[0][0];
    }

    // Formula for 2x2-matrix
    if (dimension == 2)
    {
        return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
    }

    double result = 0;
    int sign = 1;
    for (int i = 0; i < dimension; i++)
    {

        // Submatrix
        vector<vector<double>> subVect(dimension - 1, vector<double>(dimension - 1));
        for (int m = 1; m < dimension; m++)
        {
            int z = 0;
            for (int n = 0; n < dimension; n++)
            {
                if (n != i)
                {
                    subVect[m - 1][z] = vect[m][n];
                    z++;
                }
            }
        }

        // recursive call
        result = result + sign * vect[0][i] * getDeterminant(subVect);
        sign = -sign;
    }

    return result;
}

vector<vector<double>> getTranspose(const vector<vector<double>> matrix1)
{

    // Transpose-matrix: height = width(matrix), width = height(matrix)
    vector<vector<double>> solution(matrix1[0].size(), vector<double>(matrix1.size()));

    // Filling solution-matrix
    for (size_t i = 0; i < matrix1.size(); i++)
    {
        for (size_t j = 0; j < matrix1[0].size(); j++)
        {
            solution[j][i] = matrix1[i][j];
        }
    }
    return solution;
}

vector<vector<double>> getCofactor(const vector<vector<double>> vect)
{
    if (vect.size() != vect[0].size())
    {
        throw runtime_error("Matrix is not quadratic");
    }

    vector<vector<double>> solution(vect.size(), vector<double>(vect.size()));
    vector<vector<double>> subVect(vect.size() - 1, vector<double>(vect.size() - 1));

    for (size_t i = 0; i < vect.size(); i++)
    {
        for (size_t j = 0; j < vect[0].size(); j++)
        {

            int p = 0;
            for (size_t x = 0; x < vect.size(); x++)
            {
                if (x == i)
                {
                    continue;
                }
                int q = 0;

                for (size_t y = 0; y < vect.size(); y++)
                {
                    if (y == j)
                    {
                        continue;
                    }

                    subVect[p][q] = vect[x][y];
                    q++;
                }
                p++;
            }
            solution[i][j] = pow(-1, i + j) * getDeterminant(subVect);
        }
    }
    return solution;
}

vector<vector<double>> getInverse(const vector<vector<double>> vect)
{
    double d = 1.0 / getDeterminant(vect);
    vector<vector<double>> solution(vect.size(), vector<double>(vect.size()));

    for (size_t i = 0; i < vect.size(); i++)
    {
        for (size_t j = 0; j < vect.size(); j++)
        {
            solution[i][j] = vect[i][j];
        }
    }

    solution = getTranspose(getCofactor(solution));

    for (size_t i = 0; i < vect.size(); i++)
    {
        for (size_t j = 0; j < vect.size(); j++)
        {
            solution[i][j] *= d;
        }
    }

    return solution;
}

// end SO copy

vector<vector<double>> matrix_mult(
    vector<vector<double>> a, vector<vector<double>> b)
{
    vector<vector<double>> mult;
    int r1 = a.size();
    int r2 = b.size();
    int c1 = a[0].size();
    int c2 = b[0].size();
    for (int i = 0; i < r1; ++i)
    {
        vector<double> tmp(c2, 0.0);
        mult.push_back(tmp);
    }

    for (int i = 0; i < r1; ++i)
        for (int j = 0; j < c2; ++j)
            for (int k = 0; k < c1; ++k)
            {
                mult[i][j] += a[i][k] * b[k][j];
            }
    return mult;
}

int get_nth_element_neq(vector<Node *> nodes, int n, int filter_type)
{
    int found = 0;
    for (int i = 0; i < nodes.size(); i++)
    {
        if (nodes[i]->getType() != filter_type)
        {
            found++;
        }
        if (found == n)
        {
            return i;
        }
    }
    return -1;
}
int get_nth_element_eq(vector<Node *> nodes, int n, int select_type)
{
    int found = 0;
    for (int i = 0; i < nodes.size(); i++)
    {
        if (nodes[i]->getType() == select_type)
        {
            found++;
        }
        if (found == n)
        {
            return i;
        }
    }
    return -1;
}
vector<vector<double>> delete_first_row_and_col(vector<vector<double>> matrix)
{
    vector<vector<double>> result;
    for (int i = 0; i < matrix.size(); i++)
    {
        if (i == 0)
        {
            continue;
        }
        vector<double> tmp;
        for (int j = 0; j < matrix[0].size(); j++)
        {
            if (j != 0)
            {
                tmp.push_back(matrix[i][j]);
            }
        }
        result.push_back(tmp);
    }
    return result;
}

int NRSSolver::solveNewtonRaphson(
    int numberOfNodes, vector<Node *> nodes, vector<Line> lines, int maxNumberOfIterations, double epsilon)
{

    int num_PV_nodes = 0;
    int num_PQ_nodes = 0;
    for (int i = 0; i < nodes.size(); i++)
    {
        if (nodes[i]->getType() == 2)
        {
            num_PV_nodes++;
        }
        if (nodes[i]->getType() == 1)
        {
            num_PQ_nodes++;
        }
    }
    int iteration = 0;
    vector<XElement> x;
    for (int i = 0; i < nodes.size(); i++)
    {
        if (nodes[i]->getType() == 0)
        {
            for (int j = 1; j < nodes.size(); j++)
            {
                x.push_back(XElement{j, nodes[i]->getD()});
            }
            for (int j = 0; j < nodes.size(); j++)
            {
                if (nodes[j]->getType() == 2)
                {
                    x.push_back(XElement{j, 1.0});
                }
                if (nodes[j]->getType() == 1)
                {
                    nodes[j]->setV(1.0);
                }
            }
            break;
        }
    }
    Eigen::MatrixXcd admittanceMatrix(numberOfNodes, numberOfNodes);
    for (int i = 0; i < numberOfNodes; i++)
    {
        for (int j = 0; j < numberOfNodes; j++)
        {
            admittanceMatrix(i, j) = 0;
        }
    }
    for (int n = 0; n < lines.size(); n++)
    {
        if (lines[n].getStartNode() == 0)
        {
            admittanceMatrix(lines[n].getEndNode() - 1, lines[n].getEndNode() - 1) += 1. / lines[n].getZ();
        }
        else if (lines[n].getEndNode() == 0)
        {
            admittanceMatrix(lines[n].getStartNode() - 1, lines[n].getStartNode() - 1) += 1. / lines[n].getZ();
        }
        else
        {
            admittanceMatrix(lines[n].getStartNode() - 1, lines[n].getStartNode() - 1) += 1. / lines[n].getZ();
            admittanceMatrix(lines[n].getEndNode() - 1, lines[n].getEndNode() - 1) += 1. / lines[n].getZ();
            admittanceMatrix(lines[n].getStartNode() - 1, lines[n].getEndNode() - 1) -= 1. / lines[n].getZ();
            admittanceMatrix(lines[n].getEndNode() - 1, lines[n].getStartNode() - 1) -= 1. / lines[n].getZ();
        }
    }
    vector<double> p_sch, q_sch;
    for (int j = 0; j < numberOfNodes; j++)
    {
        if (nodes[j]->getType() == 2)
        {
            p_sch.push_back(nodes[j]->getP());
        }
        else if (nodes[j]->getType() == 1)
        {
            p_sch.push_back(nodes[j]->getP());
            q_sch.push_back(nodes[j]->getQ());
        }
    }
    while (iteration < maxNumberOfIterations)
    {
        vector<double> P, delta_P;
        vector<double> Q, delta_Q;
        for (int j = 0; j < numberOfNodes; j++)
        {
            struct PQ tmp = getNodePQ(nodes[j], nodes, numberOfNodes, admittanceMatrix);
            if (nodes[j]->getType() == 0)
            {
                nodes[j]->setP(tmp.p[0]);
                nodes[j]->setQ(tmp.q[0]);
                continue;
            }
            nodes[j]->setP(tmp.p[0]);
            nodes[j]->setQ(tmp.q[0]);
            for (int k = 0; k < tmp.p.size(); k++)
            {
                P.push_back(tmp.p[k]);
            }
            for (int k = 0; k < tmp.q.size(); k++)
            {
                Q.push_back(tmp.q[k]);
            }
        }
        for (int i = 0; i < p_sch.size(); i++)
        {
            delta_P.push_back(p_sch[i] - P[i]);
        }
        for (int i = 0; i < q_sch.size(); i++)
        {
            delta_Q.push_back(q_sch[i] - Q[i]);
        }
        
        vector<vector<double>> J1 = getJ1(nodes, numberOfNodes, admittanceMatrix);
        vector<vector<double>> J2 = getJ2(nodes, numberOfNodes, admittanceMatrix);
        vector<vector<double>> J3 = getJ3(nodes, numberOfNodes, admittanceMatrix);
        vector<vector<double>> J4 = getJ4(nodes, numberOfNodes, admittanceMatrix);

        vector<vector<double>> J(J1.size() + J3.size());
        for (int i = 0; i < J1.size(); i++)
        {
            J[i].insert(J[i].end(), J1[i].begin(), J1[i].end());
            J[i].insert(J[i].end(), J2[i].begin(), J2[i].end());
        }
        int tmp = J1.size();
        for (int i = 0; i < J3.size(); i++)
        {
            J[tmp + i].insert(J[tmp + i].end(), J3[i].begin(), J3[i].end());
            J[tmp + i].insert(J[tmp + i].end(), J4[i].begin(), J4[i].end());
        }
        

        vector<vector<double>> inv_j = getInverse(J);
        vector<double> delta_vector_row;
        for (int i = 0; i < delta_P.size(); i++)
        {
            delta_vector_row.push_back(delta_P[i]);
        }
        for (int i = 0; i < delta_Q.size(); i++)
        {
            delta_vector_row.push_back(delta_Q[i]);
        }
        vector<vector<double>> delta_vector;
        for (int i = 0; i < delta_vector_row.size(); i++)
        {
            vector<double> tmp{delta_vector_row[i]};
            delta_vector.push_back(tmp);
        }
        
        vector<vector<double>> result = matrix_mult(inv_j, delta_vector);

        // // update results
        for (int i = 0; i < numberOfNodes; i++)
        {
            x[i].value += result[i][0];
        }
        int n = 1; // get 1st element
        for (int i = 0; i < num_PQ_nodes; i++)
        {
            int update_index = get_nth_element_eq(nodes, n, 1);
            if (update_index >= 0)
            {
                nodes[update_index]->setD(x[i].value);
                n++;
            }
        }
        n = 1;
        for (int i = num_PQ_nodes; i < num_PV_nodes + num_PQ_nodes; i++)
        {
            int update_index = get_nth_element_eq(nodes, n, 2);
            if (update_index >= 0)
            {
                nodes[update_index]->setD(x[i].value);
                n++;
            }
        }
        n = 1;
        for (int i = num_PV_nodes + num_PQ_nodes; i < x.size(); i++)
        {
            int update_index = get_nth_element_eq(nodes, n, 1);
            if (update_index >= 0)
            {
                nodes[update_index]->setV(x[i].value);
                n++;
            }
        }

        double p_residual = 0;
        double q_residual = 0;
        int num_nok = 0;
        for (int i = 0; i < delta_P.size(); i++)
        {
            if (delta_P[i] > epsilon)
            {
                num_nok++;
                break;
            }
        }
        for (int i = 0; i < delta_Q.size(); i++)
        {
            if (delta_Q[i] > epsilon)
            {
                num_nok++;
                break;
            }
        }
        if (num_nok == 0)
        {
            break;
        }
        iteration++;
    }
    this->numberOfIter = iteration;
    return getNumberOfIter();
}


void NRSSolver::setNumberOfIter(int iter)
{
    NRSSolver::numberOfIter = iter;
}