#include "App.h"
using namespace std;

void print_node_list(vector<Node *> &nodes, string line_start="")
{
    for (int i = 0; i < nodes.size(); i++)
    {
        cout << line_start << nodes[i]->print_node() << endl;
    }
}
void print_node_list(vector<Node> &nodes, string line_start="")
{
    for (int i = 0; i < nodes.size(); i++)
    {
        cout << line_start << nodes[i].print_node() << endl;
    }
}
void print_line_list(vector<Line> &lines, string line_start="")
{
    for (int i = 0; i < lines.size(); i++)
    {
        cout << line_start << lines[i].getStartNode() << "->" << lines[i].getEndNode() << ":  " << lines[i].getZ() << endl;
    }
    cout << endl;
}

App::App()
    : window_grid(), mesh(), node_info(), conection_info(), nodes(), selected_node(nullptr), node_list(), info_grid(), data_grid()
{
    // Gtk::Fixed mesh;
    // mesh.set_hexpand(true);
    // app->set_accel_for_action("example.new", "<Primary>n");
    Gtk::Box * ToolBar = new Gtk::Box();
    ToolBar->set_style(Gtk::STYLE_CLASS_PRIMARY_TOOLBAR);
    // ToolBar->append("New", "New", "New", "New");
    // window_grid.attach(*ToolBar, 0, 0, 1, 1);
    build_node_info(node_info);
    build_conection_info(conection_info);

    listbox_nodes = new Gtk::ListBox();
    listbox_nodes->signal_row_selected().connect([this](Gtk::ListBoxRow *row)
                                                 { this->change_selected_node(row); });

    new_node(slack);
    listbox_nodes->set_hexpand(true);
    node_info.set_hexpand(true);
    conection_info.set_hexpand(true);

    // window_grid.attach(mesh, 0, 0, 1, 1);
    window_grid.attach(*listbox_nodes, 0, 0, 1, 1);
    window_grid.attach(node_info, 1, 0, 1, 1);
    window_grid.attach(conection_info, 2, 0, 1, 1);
    set_child(window_grid);
}

App::~App()
{
}

void App::build_node_list(Gtk::Grid &node_list)
{
}

void App::build_node_info(Gtk::Grid &node_info)
{
    node_info.set_vexpand(true);
    node_info.set_row_spacing(10);
    node_info.set_column_spacing(10);
    node_info.set_margin(10);

    Gtk::Box node_god;
    btn_add_node = Gtk::Button("Dodaj novi cvor");
    btn_add_node.signal_clicked().connect([this]()
                                          { this->new_node(); });
    node_god.append(btn_add_node);

    btn_delete_node = Gtk::Button("Obrisi zadnji cvor");
    btn_delete_node.signal_clicked().connect([this]()
                                             { this->delete_node(); });
    node_god.append(btn_delete_node);
    node_info.attach(node_god, 0, 0, 1, 1);

    row_label_node_number = Gtk::ListBoxRow();
    node_info.attach(row_label_node_number, 0, 1, 1, 1);

    // Node type selector
    Gtk::Box box_type;
    for (int i = 0; i < 3; i++)
    {
        type[i] = Gtk::CheckButton(tipovi[i]);
        if(i!=0) type[i].set_group(type[0]);
        type[i].signal_toggled().connect([this]()
                                         { this->change_selection(); });
        box_type.append(type[i]);
    }
    type[1].set_active(true);
    node_info.attach(box_type, 0, 2, 1, 1);

    auto napon = new Gtk::Label("Napon: ");
    input_napon = new Gtk::Entry();
    data_grid.attach(*napon, 0, 1, 1, 1);
    data_grid.attach(*input_napon, 1, 1, 1, 1);

    auto ugao = new Gtk::Label("Ugao: ");
    input_ugao = new Gtk::Entry();
    data_grid.attach(*ugao, 0, 2, 1, 1);
    data_grid.attach(*input_ugao, 1, 2, 1, 1);

    auto P = new Gtk::Label("P: ");
    input_P = new Gtk::Entry();
    data_grid.attach(*P, 0, 3, 1, 1);
    data_grid.attach(*input_P, 1, 3, 1, 1);

    auto Q = new Gtk::Label("Q: ");
    input_Q = new Gtk::Entry();
    data_grid.attach(*Q, 0, 4, 1, 1);
    data_grid.attach(*input_Q, 1, 4, 1, 1);

    btn_solve = Gtk::Button("Izracunaj");
    btn_solve.signal_clicked().connect([this]()
                                       { this->solve(); });
    data_grid.attach(btn_solve, 0, 5, 2, 1);

    node_info.attach(data_grid, 0, 3, 1, 1);
}

void App::change_selection()
{
    for (int i = 0; i < 3; i++)
    {
        if (type[i].get_active())
        {
            if (nodes.empty())
                break;
            selected_node->setType(i);
            if (DEBUG)
                cout << "Tip: " << i << endl;
            break;
        }
    }
}

void App::new_node(int tip)
{
    nodes.push_back(new Node(nodes.size(), tip));
    if (DEBUG)
        cout << "Dodan cvor " << nodes.back()->getNumber() << " tip: " << nodes.back()->getType() << endl;
    listbox_nodes->append(nodes.back()->getRow());
    if (nodes.size() == 1)
    {
        listbox_nodes->select_row(nodes.back()->getRow());
        change_selected_node(&nodes.back()->getRow());
    }
    save_conections();
    populate_conections();
}

void App::delete_node()
{
    // nema nista da se obrise, prazan queue
    if (DEBUG && nodes.empty())
        cout << "Nema cvorova za brisanje" << endl;
    if (nodes.empty())
        return;
    // hocemo da ostane barem jedan cvor
    if (DEBUG && nodes.size() == 1)
        cout << "Zadnji cvor se ne moze obrisati" << endl;
    if (nodes.size() == 1)
        return;

    if (DEBUG)
        cout << "Obrisan cvor " << nodes.back()->getNumber();
    listbox_nodes->remove(nodes.back()->getRow());
    delete nodes.back();
    nodes.pop_back();
    if (DEBUG)
        cout << " Broj cvorova: " << nodes.size() << endl;
    if (nodes.empty())
        selected_node = nullptr;
    if (listbox_nodes->get_selected_row() == nullptr)
    {
        listbox_nodes->select_row(nodes.back()->getRow());
        change_selected_node(&nodes.back()->getRow());
    }
    save_conections();
    populate_conections();
}

// potrebno je spasiti unesene vrijednosti u selected_node
// zatim promijeniti selected_node na novi cvor
// zetim ucitati vrijednosti iz selected_node u input_napon, input_ugao, input_P, input_Q
void App::change_selected_node(Gtk::ListBoxRow *row)
{
    if (DEBUG && listbox_nodes->get_selected_row() == nullptr)
        cout << "Nema selektovanog cvora" << endl;
    if (listbox_nodes->get_selected_row() == nullptr)
        return;

    // pocetak programa
    if (selected_node == nullptr)
    {
        if (DEBUG)
            cout << "Pocetak programa, prvi cvor" << endl;
        selected_node = nodes.back();
        return;
    }

    save_state();
    selected_node = nodes[listbox_nodes->get_selected_row()->get_index()];
    if (DEBUG)
        cout << "Selektovan cvor: " << selected_node->getNumber() << endl;
    row_label_node_number.set_child(*Gtk::manage(new Gtk::Label("Cvor " + to_string(selected_node->getNumber()) + " Info")));
    if (DEBUG)
        cout << "Ucitavam vrijednosti iz cvora: " << selected_node->getNumber() << endl;
    input_napon->set_text(double_to_string(selected_node->getV()));
    input_ugao->set_text(double_to_string(selected_node->getD()));
    input_P->set_text(double_to_string(selected_node->getP()));
    input_Q->set_text(double_to_string(selected_node->getQ()));
    type[selected_node->getType()].set_active(true);

    populate_conections();
}

string App::double_to_string(double d)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(3) << d;
    return ss.str();
}

void App::build_conection_info(Gtk::Grid &conection_info)
{
    conection_info.set_vexpand(true);
    conection_info.set_row_spacing(10);
    conection_info.set_column_spacing(10);
    conection_info.set_margin(10);
    Gtk::Label label("Impedanse izmedju cvorova");
    conection_info.attach(label, 0, 0, 1, 1);

    listbox_connections = new Gtk::ListBox();
    gtk_list_box_set_selection_mode(listbox_connections->gobj(), GTK_SELECTION_NONE);
    conection_info.attach(*listbox_connections, 0, 1, 1, 1);
}

void App::populate_conections()
{
    if (DEBUG)
        cout << endl;
    while (listbox_connections->get_last_child() != nullptr)
    {
        listbox_connections->remove(*listbox_connections->get_last_child());
    }

    if (DEBUG)
        cout << "Obrisane konekcije" << endl;

    if (DEBUG && nodes.empty())
        cout << "Nema cvorova za konekcije" << endl;
    if (nodes.empty())
        return;
    if (DEBUG && nodes.size() == 1)
        cout << "Zadnji cvor se ne moze konektovati" << endl;
    if (nodes.size() == 1)
        return;
    if (DEBUG)
        cout << "Popunjavam konekcije" << endl;
    for (int i = 0; i < nodes.size(); i++)
    {
        if (DEBUG)
            cout << "Konektujem cvor " << selected_node->getNumber() << " sa cvorom " << nodes[i]->getNumber() << "; i = " << i << endl;

        if (selected_node->getNumber() == nodes[i]->getNumber())
        {
            if (DEBUG)
                cout << "Ne moze se konektovati sa samim sobom" << endl;
            continue;
        }

        listbox_connections->append(selected_node->getConnectionRow(nodes[i]));
    }
    if (DEBUG)
        cout << endl;
}
void App::save_node_param()
{
    if (DEBUG)
        cout << "Spasavam vrijednosti u cvor: " << selected_node->getNumber() << endl;
    selected_node->setV(atof(input_napon->get_text().data()));
    selected_node->setD(atof(input_ugao->get_text().data()));
    selected_node->setP(atof(input_P->get_text().data()));
    selected_node->setQ(atof(input_Q->get_text().data()));
}
void App::save_conections()
{
    if (DEBUG)
        cout << endl
             << "Spasavam unesene vrijednosti" << endl;
    for (int i = 0; i < nodes.size(); i++)
    {
        if (selected_node->getNumber() == nodes[i]->getNumber())
        {
            continue;
        }
        if (DEBUG)
            cout << selected_node->getNumber() << " -> " << nodes[i]->getNumber() << ": ";
        selected_node->save_impedance(nodes[i]);
    }
    if (DEBUG)
        cout << endl;
}

void App::save_state() {
    if (DEBUG)
        cout << "Spasavam stanje" << endl;
    save_node_param();
    save_conections();
}

void App::solve()
{
    save_state();
    
    if (DEBUG)
        cout << "Solver pokrenut" << endl;
    GSSolver solver;
    NRSSolver solver2;
    vector<Node *> GS_nodes;
    vector<Node *> NR_nodes;
    double epsilon = pow(10, -5);
    int max_iter = 100;
    vector<Line> lines;
    

    for (int i = 0; i < nodes.size(); i++)
    {
        GS_nodes.push_back(new Node(*nodes[i]));
        NR_nodes.push_back(new Node(*nodes[i]));
    }
    cout << "input: " << endl << endl;
    print_node_list(GS_nodes);
    cout << endl;

    // create lines from node connections
    for (int i = 0; i < nodes.size(); i++)
    {
        for (int j = 0; j < i && j < nodes.size(); j++)
        {
            if (i == j)
                continue;
            try
            {
                lines.push_back(Line(nodes[i]->getNumber()+1, nodes[j]->getNumber()+1,
                                     nodes[i]->getImpedance(nodes[j]->getNumber())));
            }
            catch (const char *msg)
            {
                continue;
            }
        }
    }

    print_line_list(lines);

    if (DEBUG)
        cout << "Pozivam GSsolver" << endl;

    auto start_GS = chrono::high_resolution_clock::now();
    int iter_num_GS = solver.solveGS(GS_nodes.size(), GS_nodes, lines, max_iter, epsilon);
    auto end_GS = chrono::high_resolution_clock::now();

    cout << "output: " << endl;
    cout << "\tGS metod: " << endl;
    cout << "\t\tBroj iteracija: " << iter_num_GS << endl;
    cout << "\t\tVrijeme izvrsavanja: " << chrono::duration_cast<chrono::microseconds>(end_GS - start_GS).count() << "us" << endl;
    print_node_list(GS_nodes, "\t\t");
    cout << endl;
   
    if (DEBUG) 
        cout << "Pozivam NRsolver" << endl;

    auto start_NR = chrono::high_resolution_clock::now();
    int iter_num_NR = solver2.solveNewtonRaphson(NR_nodes.size(), NR_nodes, lines, max_iter, epsilon);
    auto end_NR = chrono::high_resolution_clock::now();
    // drugi solver..
    cout << "\tNR metod: " << endl;
    cout << "\t\tBroj iteracija: " << iter_num_NR << endl;
    cout << "\t\tTrajanje: " << chrono::duration_cast<chrono::microseconds>(end_NR - start_NR).count() << "us" << endl;
    print_node_list(NR_nodes, "\t\t");
    cout << endl;


    cout << "\n----------------------------------------\n" << endl;
    
}

