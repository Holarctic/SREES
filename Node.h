#ifndef NODE_H
#define NODE_H

#ifndef DEBUG
#define DEBUG 0
#endif // DEBUG

#include <map>
#include <complex>
#include "gtkmm.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

using namespace std;
class Node {
private:
    int number, type;
    /*
    * types:
      SLACK = 0,
	  LOAD/PQ = 1,
	  GEN/PV = 2,
    */
    double V, d, P, Q;
    const string tipovi[3] = {"Slack", "PQ", "PV"};
    std::complex<double> Z;
    /*
     P = Net Real Power on the Bus
     Q = Net Reactive Power on the Bus
     V = Bus Voltage
     d = Bus Voltage angle
    */
   std::map<int, std::complex<double>> impedance;
   std::map<int, Gtk::ListBoxRow> connection_row;
   std::map<int, Gtk::Entry> connection_entry;
   /*
    * impedanse[Node] = impedansa izmedju ovog i drugog cvora
    *                   veza je ostvarena pomocu pointera, treba napraviti i u drugom Node jer je impedansa simetricna izmedju
    */
   Gtk::ListBoxRow node_row;
   void build_row();

   std::complex<double> complex_V;
   
public:
    // Node();
    Node(int number, int type);
    Node(int n_number, int n_type, double n_V, double n_d, double n_P, double n_Q);
    Node(Node&);
    void SetNode(int n_number, int n_type, double n_V, double n_d, double n_P, double n_Q);

    int getNumber() { return number; }
    int getType(){ return type;}
    double getV() { return V; }
    double getD() { return d; }
    double getP() { return P; }
    double getQ() { return Q; }
    std::complex<double> getZ() { return Z; }

    void setV(double V) { this->V = V; }
    void setD(double d) { this->d = d; }
    void setP(double P) { this->P = P; }
    void setQ(double Q) { this->Q = Q; };
    void setZ(std::complex<double> Z) { this->Z = Z; }
    void setType(int type);
    void setImpedance(Node*, std::complex<double>);
    std::complex<double> getImpedance(int node);

    Gtk::ListBoxRow & getRow() { return node_row; }
    Gtk::ListBoxRow & getConnectionRow(Node* node);

    void save_impedance(Node*);
    string double_to_string(double);
    complex<double> string_to_complex(string );
    string complex_to_string(complex<double>, int );

    string print_node();

};
#endif