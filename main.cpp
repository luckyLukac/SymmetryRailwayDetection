#include <QtWidgets/QApplication>

#include "MainWindow.hpp"


int main(int argc, char *argv[]) {
    // Context and OpenGL initialization.
    // Version: OpenGL 3.3
    // Profile: Core
    QSurfaceFormat glFormat;
    glFormat.setVersion(3, 3);
    glFormat.setProfile(QSurfaceFormat::CoreProfile);
    glFormat.setSamples(8);
    QSurfaceFormat::setDefaultFormat(glFormat);

    // Creating the main window of the application.
    QApplication::setAttribute(Qt::AA_UseDesktopOpenGL);
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
