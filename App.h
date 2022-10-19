#ifndef GTKMM_EXAMPLE_APP_H
#define GTKMM_EXAMPLE_APP_H

#ifndef DEBUG
#define DEBUG 0
#endif // DEBUG

#include <gtkmm.h>
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <chrono>

#include "Node.h"
#include "GSSolver.h"
#include "Line.h"
#include "NRSolver.h"


using namespace std;

class App : public Gtk::Window
{

public:
    App();
    virtual ~App();

protected:
    // Signal handlers:

    // Member widgets:
    Gtk::Grid window_grid;

    Gtk::Fixed mesh;
    Gtk::Grid node_info;
    Gtk::Grid conection_info;
private:
    const int slack{0}, pq{1}, pv{2};  // types of nodes
    const string tipovi[3] = {"Slack", "PQ", "PV"};

    Gtk::Button btn_add_node;
    Gtk::Button btn_delete_node;

    Gtk::ListBoxRow row_label_node_number;
    void new_node(int tip=1);
    void delete_node();

    std::vector<Node*> nodes;
    Node * selected_node;

    void build_node_list(Gtk::Grid &);
    Gtk::Grid node_list;
    Gtk::ListBox *listbox_nodes;

    Gtk::ListBox *listbox_connections;

    void build_node_info(Gtk::Grid &);
    
    Gtk::Grid info_grid;
    Gtk::Grid data_grid;
    Gtk::CheckButton type[3];
    Gtk::Entry *input_napon, *input_ugao, *input_P, *input_Q;
    Gtk::Button btn_solve;
    void build_conection_info(Gtk::Grid &);

    void change_selected_node(Gtk::ListBoxRow *);

    void change_selection();
    void populate_conections();
    void save_node_param();
    void save_conections();
    string double_to_string(double);
    void solve();
    void save_state();
};

#endif