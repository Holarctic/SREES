#include "App.h"
#include <gtkmm/application.h>
#include <iostream>

int main(int argc, char *argv[])
{
  auto app = Gtk::Application::create("org.gtkmm.example");
  // app->set_accel_for_action("example.new", "<Primary>n");
  // Shows the window and returns when it is closed.
  return app->make_window_and_run<App>(argc, argv);
}