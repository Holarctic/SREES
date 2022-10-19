#include "Node.h"

Node::Node(int n_number, int n_type, double n_V, double n_d, double n_P, double n_Q)
{
    SetNode(n_number, n_type, n_V, n_d, n_P, n_Q);
}

void Node::SetNode(int n_number, int n_type, double n_V, double n_d, double n_P, double n_Q){
    this->number = n_number;
    this->type = n_type;
    this->V = n_V;
    this->d = n_d;
    this->P = n_P;
    this->Q = n_Q;
}

Node::Node(int number=0, int type=0) : 
    number(number), 
    type(type), 
    V(1),
    d(0),
    P(0),
    Q(0),
    node_row(Gtk::ListBoxRow())
{
    build_row();
}

void Node::setType(int type) {
    this->type = type;
    build_row();
}


void Node::setImpedance(Node* node, std::complex<double> impedance)
{
    this->impedance[node->getNumber()] = impedance;
    node->impedance[this->getNumber()] = impedance;
}

void Node::build_row() {
    auto label = Gtk::Label("Cvor " + std::to_string(number) + "\t\tTip: " + tipovi[type]);
    node_row.set_child(*Gtk::manage(&label));
}


Gtk::ListBoxRow & Node::getConnectionRow(Node* node)
{
    if (node->getNumber() == this->getNumber()) {
        throw "Ne moze se napraviti veza sa samim sobom";
    }
    if (connection_row.find(node->getNumber()) == connection_row.end()) {
        auto label = Gtk::Label("Veza " + to_string(getNumber()) + "-" + std::to_string(node->getNumber()) + ":    ");
        connection_entry[node->getNumber()] = Gtk::Entry();
        auto box = Gtk::Box();
        box.set_size_request(200);
        box.append(label);
        box.append(connection_entry[node->getNumber()]);
        connection_row[node->getNumber()] = Gtk::ListBoxRow();
        connection_row[node->getNumber()].set_child(*Gtk::manage(&box));
    }
    connection_entry[node->getNumber()].set_text(node->connection_entry[getNumber()].get_text());
    return connection_row[node->getNumber()];
}

void Node::save_impedance(Node* node)
{
    try {
        complex<double> z = string_to_complex(connection_entry[node->getNumber()].get_text());
        setImpedance(node, z);
        if (DEBUG) cout << endl;
    } catch (const char* msg) {
        if (DEBUG) cout << msg << endl;
        return;
    }
}

complex<double> Node::string_to_complex(string s) {
    if (DEBUG) cout << "-" << s << "- "; 
    if (s.empty()) throw "Prazan string";
    if (s.at(0)  != '(') s = "(" + s;
    if (s.back() != ')') s = s + ")";
    replace(s.begin(), s.end(), '.', ','); // cuz of the locale, ocekuje brojeve sa ","
    istringstream is(s);
    
    complex<double> c;  
    is >> c;
    if (!is) throw "Neispravan format kompleksnog broja";
    return c;
}

string Node::complex_to_string(complex<double> c, int precision=2) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << c;
    return ss.str();
}

std::complex<double> Node::getImpedance(int node) {
    if (impedance.find(node) == impedance.end()) {
        throw "None";
    }
    return impedance[node];
}


Node::Node(Node& orig) {
    this->number = orig.number+1;
    this->type = orig.type;
    this->V = orig.V;
    this->d = orig.d;
    this->P = orig.P;
    this->Q = orig.Q;
}

string Node::print_node() {
    string s = "Cvor " + to_string(number-1) + ":\t\t";
    s += "V = " + to_string(V) + " p.u.\t";
    s += "d = " + to_string(d) + " p.u.\t";
    s += "P = " + to_string(P) + " p.u.\t";
    s += "Q = " + to_string(Q) + " p.u.\t";
    return s;
}