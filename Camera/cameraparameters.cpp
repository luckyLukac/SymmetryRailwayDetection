#include "CameraParameters.hpp"
#include "ui_CameraParameters.h"


// Constructor of the form.
CameraParameters::CameraParameters(const Camera& camera, QWidget* parent) :
    QDialog(parent),
    ui(new Ui::CameraParameters),
    FOV(camera.FOV),
    cameraPosition(camera.position),
    cameraRotation(camera.rotation)
{
    ui->setupUi(this);

    // Setting all parameters in the GUI.
    ui->sbx_FOV->setValue(camera.FOV);
    ui->sbx_CameraPositionX->setValue(-cameraPosition[0]);
    ui->sbx_CameraPositionY->setValue(cameraPosition[2]);
    ui->sbx_CameraPositionZ->setValue(-cameraPosition[1]);
    ui->sbx_CameraRotationX->setValue(-cameraRotation[1]);
    ui->sbx_CameraRotationY->setValue(cameraRotation[0]);
    ui->sbx_CameraRotationZ->setValue(cameraRotation[2]);
}

// Getting all parameters of the camera.
Camera CameraParameters::getParams() {
    // Retrieving new values.
    FOV = ui->sbx_FOV->value();
    cameraPosition[0] = ui->sbx_CameraPositionX->value();
    cameraPosition[1] = ui->sbx_CameraPositionY->value();
    cameraPosition[2] = ui->sbx_CameraPositionZ->value();
    cameraRotation[0] = ui->sbx_CameraRotationX->value();
    cameraRotation[1] = ui->sbx_CameraRotationY->value();
    cameraRotation[2] = ui->sbx_CameraRotationZ->value();

    return Camera(FOV, cameraPosition, cameraRotation);
}

// Destructor of the form.
CameraParameters::~CameraParameters() {
    delete ui;
}
