#include <QApplication>

#include "mainwindow.h"

int main(int argc, char* argv[])
{
  QApplication application(argc, argv);
  QApplication::setApplicationName("Ruptura Lab");
  QApplication::setOrganizationName("Ruptura");

  MainWindow window;
  window.show();

  return QApplication::exec();
}
