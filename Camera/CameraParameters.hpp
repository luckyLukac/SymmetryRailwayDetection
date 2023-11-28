#pragma once

#include <QDialog>


namespace Ui {
    class CameraParameters;
}


/// <summary>
/// Camera parameters structure.
/// </summary>
struct Camera {
    double FOV;
    std::vector<double> position;
    std::vector<double> rotation;

    Camera(const double FOV, const std::vector<double>& position, const std::vector<double>& rotation) : FOV(FOV), position(position), rotation(rotation) {};
};


/// <summary>
/// Class for manipulating with camera parameters.
/// </summary>
class CameraParameters : public QDialog {
    Q_OBJECT
private:
    Ui::CameraParameters* ui;            // Form data structure.
    double FOV;                          // Field of View.
    std::vector<double> cameraPosition;  // Camera position by X, Y and Z axes.
    std::vector<double> cameraRotation;  // Camera rotation by X, Y and Z axes.

public:
    explicit CameraParameters(const Camera& camera, QWidget* parent = nullptr);  // Constructor of the form.
    ~CameraParameters();                                                         // Destructor of the form.

    Camera getParams();  // Getting all parameters of the camera.
};