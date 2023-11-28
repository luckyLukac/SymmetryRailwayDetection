#include <iostream>
#include <QFileDialog>
#include <QMainWindow>
#include <QMessageBox>

#include "GeometricRailwayExtraction.hpp"
#include "MainWindow.hpp"
#include "PointReader/PointReader.hpp"
#include "RailObstacleDetector.hpp"


MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    m_Ui.setupUi(this);
    
    // SIGNAL-SLOT CONNECTIONS
    connect(m_Ui.actionAbout, SIGNAL(triggered()), this, SLOT(displayBasicInfo()));
    connect(m_Ui.actionResetCamera, SIGNAL(triggered()), this, SLOT(resetCamera()));
    connect(m_Ui.actionTop_down_view, SIGNAL(triggered()), this, SLOT(topDownView()));
    connect(m_Ui.actionCameraParameters, SIGNAL(triggered()), this, SLOT(openCameraParameters()));
    connect(m_Ui.btn_FindFileAndLoad, SIGNAL(clicked()), this, SLOT(findPointCloud()));
    connect(m_Ui.btn_LoadPointCloud, SIGNAL(clicked()), this, SLOT(loadPointCloud()));
    connect(m_Ui.btn_GeometricPointFilter, SIGNAL(clicked()), this, SLOT(applyGeometricPointFilter()));
    connect(m_Ui.btn_DetectRailways, SIGNAL(clicked()), this, SLOT(detectRailways()));

    connect(m_Ui.rb_PerspectiveProjection, SIGNAL(clicked()), this, SLOT(setPerspectiveProjection()));
    connect(m_Ui.rb_OrthoProjection, SIGNAL(clicked()), this, SLOT(setOrthogonalProjection()));
    connect(m_Ui.cbx_RenderPoints, SIGNAL(toggled(bool)), this, SLOT(renderPoints(bool)));
    connect(m_Ui.rb_RenderVoxelEdges, SIGNAL(clicked()), this, SLOT(renderVoxelsAsEdges()));
    connect(m_Ui.rb_RenderVoxelCubes, SIGNAL(clicked()), this, SLOT(renderVoxelsAsCubes()));
    connect(m_Ui.cbx_RenderVoxels, SIGNAL(toggled(bool)), this, SLOT(renderVoxels(bool)));
    connect(m_Ui.cbx_RenderSymmetryVoxels, SIGNAL(toggled(bool)), this, SLOT(renderSymmetryVoxels(bool)));
    connect(m_Ui.cbx_RenderInterestingVoxels, SIGNAL(toggled(bool)), this, SLOT(renderInterestingVoxels(bool)));
    connect(m_Ui.cbx_RenderMaterialVoxels, SIGNAL(toggled(bool)), this, SLOT(renderMaterialVoxels(bool)));
    connect(m_Ui.cbx_RenderOrdinaryVoxels, SIGNAL(toggled(bool)), this, SLOT(renderOrdinaryVoxels(bool)));
    connect(m_Ui.sbx_ReflectionSymmetries, SIGNAL(currentIndexChanged(int)), this, SLOT(renderSymmetry(int)));
    connect(m_Ui.cbx_RenderSymmetryPlane, SIGNAL(toggled(bool)), this, SLOT(renderSymmetryPlane(bool)));
    connect(m_Ui.btn_TopDownView, SIGNAL(clicked()), this, SLOT(railwayTopDownView()));
}

MainWindow::~MainWindow()
{}



// PRIVATE FUNCTIONS
void MainWindow::print(const std::string& title, const std::string& body) const {
    m_Ui.tbx_Console->moveCursor(QTextCursor::End);
    m_Ui.tbx_Console->setTextColor(QColor(0, 100, 0));
    m_Ui.tbx_Console->insertPlainText(QString::fromStdString(title) + "\n");
    m_Ui.tbx_Console->setTextColor(QColor(0, 0, 0));
    m_Ui.tbx_Console->insertPlainText(QString::fromStdString(body) + "\n\n");
    m_Ui.tbx_Console->moveCursor(QTextCursor::End);
}

std::string MainWindow::pointCloudDataToString(const u128 pointCount, const VoxelMesh& voxelMesh) const {
    std::stringstream ss;

    ss << "Number of points: " << pointCount << "\n";
    ss << "Δx = " << voxelMesh.deltaX / 100 << " m \n";
    ss << "Δy = " << voxelMesh.deltaY / 100 << " m \n";
    ss << "Δz = " << voxelMesh.deltaZ / 100 << " m";

    return ss.str();
}

std::string MainWindow::railwayDetectionString(const u128 time, const uint numberOfSymmetries, const VoxelMesh& voxelMesh) const {
    std::stringstream ss;

    ss << "Voxel mesh: " << voxelMesh.xSize << " × " << voxelMesh.ySize << " × " << voxelMesh.zSize << " = " << voxelMesh.count << "\n";
    ss << "Number of railways: " << numberOfSymmetries << "\n";
    ss << "Time [s]: " << time / 1000.0 << " s\n";

    return ss.str();
}


// SLOTS
void MainWindow::displayBasicInfo() const {
    // Creating and displaying a message
    // box with the basic data.
    QMessageBox msg(
        QMessageBox::Information,
        "Authors",
        "GeMMA\n2023\n\nLuka Lukač\nDavid Podgorelec",
        QMessageBox::Close
    );
    msg.exec();
}

void MainWindow::resetCamera() const {
    m_Ui.openGLWidget->resetCamera();
}

void MainWindow::topDownView() const {
    // If there are no points, we should do exactly nothing.
    if (m_LocalSymmetry.getPoints().empty()) {
        return;
    }

    const VoxelMesh& vm = m_LocalSymmetry.getVoxelMesh();  // Getting the voxel mesh.

    // Setting camera parameters.
    m_Ui.openGLWidget->resetCamera();
    m_Ui.openGLWidget->moveCameraToPoint((vm.minX + vm.maxX) / 2, (vm.minY + vm.maxY) / 2, -vm.maxZ - 150);
    m_Ui.openGLWidget->rotateCamera(-90, 0, 0);
}

void MainWindow::openCameraParameters() const {
    auto params = m_Ui.openGLWidget->getCameraParameters();  // Getting current camera parameters.

    // Execution of the form with camera parameters.
    // If clicked OK, the new camera parameters are set.
    CameraParameters c(params);
    if (c.exec()) {
        // Getting new camera parameters.
        const Camera newCamera = c.getParams();
        const double newFOV = newCamera.FOV;
        const std::vector<double> newPositions = newCamera.position;
        const std::vector<double> newRotations = newCamera.rotation;

        // Setting new camera parameters.
        m_Ui.openGLWidget->resetCamera();
        m_Ui.openGLWidget->setFOV(newFOV);
        m_Ui.openGLWidget->moveCameraToPoint(newPositions[0], newPositions[1], -newPositions[2]);
        m_Ui.openGLWidget->rotateCamera(newRotations[1], -newRotations[0], 0);
    }
}

void MainWindow::findPointCloud() {
    // Creating a QFileDialog that accepts LAS files only.
    QFileDialog dialog;
    dialog.setFileMode(QFileDialog::ExistingFile);
    //dialog.setNameFilter(tr("LAS files (*.las)"));

    // Dialog start.
    if (dialog.exec()) {
        // Reading the path to the file from the dialog.
        if (!dialog.selectedFiles().empty()) {
            m_Ui.tbx_FileName->setText(dialog.selectedFiles().at(0));
            m_Ui.btn_LoadPointCloud->setEnabled(true);
        }
    }

    loadPointCloud();
}

void MainWindow::loadPointCloud() {
    // Reading points.
    const std::string path = m_Ui.tbx_FileName->text().toStdString();                                        // Getting the path to the LAS file.
    const std::string suffix = path.substr(path.find_last_of('.') + 1, path.size() - path.find_last_of('.'));

    if (suffix == "pc" || suffix == "csv") {
        auto readPoints = PointReader::readPointsFromPointCloud(path);
        m_LocalSymmetry.setPoints(readPoints);
        m_Points = m_LocalSymmetry.getPoints();
    }

    // Calculating a voxel mesh and rendering points.
    const VoxelMesh& voxelMesh = m_LocalSymmetry.getVoxelMesh();                                   // Getting voxel mesh to find extreme points.
    const std::vector<PointPosition> positions(m_Points.size(), PointPosition::notInSymmetry);     // Positions from a non-existent plane in the beginning are undefined.
    m_Ui.openGLWidget->setPoints(m_Points, positions);                                             // Rendering points in OpenGL.
    m_Ui.openGLWidget->setVoxelMesh(voxelMesh);                                                    // Voxel mesh setting in OpenGL.
    m_Ui.openGLWidget->setRenderSymmetryAxis(false);                                               // Disabling symmetry axis rendering in OpenGL.
    m_Ui.openGLWidget->setRenderSymmetryPlane(false);                                              // Disabling symmetry plane rendering in OpenGL.

    // Basic description dump in the console and GUI manipulation.
    print("POINT CLOUD", pointCloudDataToString(m_Points.size(), voxelMesh));
    m_Ui.btn_GeometricPointFilter->setEnabled(true);
}

void MainWindow::applyGeometricPointFilter() {
    GeometricRailwayExtraction geometricFilter = GeometricRailwayExtraction(m_LocalSymmetry.getPoints(), m_LocalSymmetry.getVoxelMesh());
    m_FilteredPoints = geometricFilter.extractCandidateRailwayPoints(10, 20);

    m_LocalSymmetry.setPoints(m_FilteredPoints, false);
    m_Ui.openGLWidget->setPointPositions(m_FilteredPoints, PointPosition::railwayCandidate);

    m_Ui.btn_GeometricPointFilter->setEnabled(false);
    m_Ui.btn_DetectRailways->setEnabled(true);
}

void MainWindow::detectRailways() {
    // Basic algorithm parameters.
    const uint voxelSideX = m_Ui.sbx_VoxelSideX->value();
    const uint voxelSideY = m_Ui.sbx_VoxelSideY->value();
    const uint voxelSideZ = m_Ui.sbx_VoxelSideZ->value();
    const double tolerance = m_Ui.sbx_Tolerance->value();
    const double angleTolerance = m_Ui.sbx_AngleTolerance->value();
    const uint minimumLineSegmentLength = m_Ui.sbx_MinimumLineSegmentLength->value();
    const uint minimumVoxelCluster = m_Ui.sbx_MinimumVoxelCluster->value();

    // Railway detection.
    auto start = std::chrono::high_resolution_clock::now();
    m_Symmetries = m_LocalSymmetry.detectReflectionSymmetries(voxelSideX, voxelSideY, voxelSideZ, tolerance, angleTolerance, minimumVoxelCluster, minimumLineSegmentLength, false);
    auto end = std::chrono::high_resolution_clock::now();

    if (m_Symmetries.empty()) {
        QMessageBox(QMessageBox::Information, "*bad program*", "No railways found on scene. Try different parameters and pray to the Lord!", QMessageBox::Ok).exec();
        return;
    }

    RailObstacleDetector obstacleDetector(m_Symmetries, m_Points, m_FilteredPoints);
    m_Symmetries = obstacleDetector.detectRails();

    // Visualizator update.
    const std::vector<PointPosition> positions(m_Points.size(), PointPosition::notInSymmetry);     // Positions from a non-existent plane in the beginning are undefined.
    m_Ui.openGLWidget->setPoints(m_Points, positions);                                             // Rendering points in OpenGL.
    m_Ui.openGLWidget->setVoxels(m_Symmetries[0].getVoxelVector(), m_Symmetries[0].getVoxelMesh(), std::vector<PointPosition>(m_Symmetries[0].getVoxelMesh().count, PointPosition::notInSymmetry));  // Setting points and rendering in OpenGL.

    // Adding symmetries to the combo box.
    m_Ui.sbx_ReflectionSymmetries->clear();
    int i = 0;
    for (const ReflectionSymmetry& symmetry : m_Symmetries) {
        m_Ui.sbx_ReflectionSymmetries->addItem("(" + QString::number(i++) + ")  " + QString::fromStdString(symmetry.getPlane().toString()));
    }

    // Console output.
    const u128 milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    print("RAILWAY DETECTION", railwayDetectionString(milliseconds, m_Symmetries.size(), m_Symmetries[0].getVoxelMesh()));
}


void MainWindow::setPerspectiveProjection() const {
    m_Ui.openGLWidget->setRenderPerspectiveProjection(true);
}

void MainWindow::setOrthogonalProjection() const {
    m_Ui.openGLWidget->setRenderPerspectiveProjection(false);
}

void MainWindow::renderPoints(const bool render) const {
    m_Ui.openGLWidget->setRenderPoints(render);
}

void MainWindow::renderVoxelsAsEdges() const {
    m_Ui.openGLWidget->setRenderVoxelsAsEdges(true);
}

void MainWindow::renderVoxelsAsCubes() const {
    m_Ui.openGLWidget->setRenderVoxelsAsEdges(false);
}

void MainWindow::renderVoxels(const bool render) const {
    m_Ui.openGLWidget->setRenderVoxels(render);
}

void MainWindow::renderSymmetryVoxels(const bool render) const {
    m_Ui.openGLWidget->setRenderSymmetryVoxels(render);
}

void MainWindow::renderInterestingVoxels(const bool render) const {
    m_Ui.openGLWidget->setRenderSuperInterestingVoxels(render);
}

void MainWindow::renderMaterialVoxels(const bool render) const {
    m_Ui.openGLWidget->setRenderInterestingVoxels(render);
}

void MainWindow::renderOrdinaryVoxels(const bool render) const {
    m_Ui.openGLWidget->setRenderOrdinaryVoxels(render);
}

void MainWindow::renderSymmetry(const int index) const {
    if (index != -1) {
        VoxelVector voxelVector = m_LocalSymmetry.getVoxelVector();
        const VoxelMesh voxelMesh = m_LocalSymmetry.getVoxelMesh();
        const uint uindex = static_cast<uint>(index);

        // Getting line segments and the plane.
        std::vector<LineSegment> lineSegments = m_Symmetries[uindex].getLineSegments();
        const Plane plane = m_Symmetries[uindex].getPlane();

        // Voxel resetting.
        std::vector<Voxel> voxels = m_Symmetries[uindex].getSymmetryVoxels();
        std::vector<Voxel> avoxels = m_Symmetries[uindex].getObstacleVoxels();
        std::vector<Voxel> avoxelsAcross = m_Symmetries[uindex].getVoxelsAcrossObstacle();
        std::vector<Voxel>& openGLvoxels = m_Ui.openGLWidget->getVoxels();
        for (Voxel& voxel : openGLvoxels) {
            voxel.inSymmetry = false;
            voxel.inAsymmetry = false;
            voxel.inAsymmetryAcross = false;
        }

        // Visualization of a voxel in the symmetry.
        for (const Voxel& v : voxels) {
            const auto [x, y, z] = v.getCoordinates();

            uint index = (
                z * voxelMesh.ySize * voxelMesh.xSize +
                y * voxelMesh.xSize +
                x
            );

            try {
                openGLvoxels[index].inSymmetry = true;
            }
            catch (const std::exception&) {

            }
        }
        for (const Voxel& v : avoxels) {
            const auto [x, y, z] = v.getCoordinates();

            uint index = (
                z * voxelMesh.ySize * voxelMesh.xSize +
                y * voxelMesh.xSize +
                x
            );

            try {
                openGLvoxels[index].inSymmetry = false;
                openGLvoxels[index].interesting = false;
                openGLvoxels[index].material = false;
                openGLvoxels[index].inAsymmetry = true;
            }
            catch (const std::exception&) {

            }

            if (z < voxelVector.size() && y < voxelVector[0].size() && x < voxelVector[0][0].size()) {
                voxelVector[z][y][x].inAsymmetry = true;
            }
        }
        for (const Voxel& v : avoxelsAcross) {
            const auto [x, y, z] = v.getCoordinates();

            uint index = (
                z * voxelMesh.ySize * voxelMesh.xSize +
                y * voxelMesh.xSize +
                x
            );


            try {
                openGLvoxels[index].inSymmetry = false;
                openGLvoxels[index].interesting = false;
                openGLvoxels[index].material = false;
                openGLvoxels[index].inAsymmetryAcross = true;
            }
            catch (const std::exception&) {

            }

            if (z < voxelVector.size() && y < voxelVector[0].size() && x < voxelVector[0][0].size()) {
                voxelVector[z][y][x].inAsymmetryAcross = true;
            }
        }

        // Setting and rendering points in OpenGL.
        const std::vector<Point>& points = m_Points;                                                                                                   // Getting LAS points.
        const std::vector<PointPosition> positions = m_Symmetries[uindex].calculatePositionsFromPlane(points, voxelVector, voxelMesh, plane, true, true);  // Positions in the beginning are neutral.
        m_Ui.openGLWidget->setPoints(points, positions);
        m_Ui.openGLWidget->setLineSegments(lineSegments, plane);
        m_Ui.openGLWidget->setRenderLineSegments(false);

        // Calculating the voxel positions.
        std::vector<Point> voxelPoints = m_LocalSymmetry.pointsFromMaterialVoxels();
        std::vector<PointPosition> voxelPositions = m_Symmetries[uindex].calculatePositionsFromPlane(voxelPoints, voxelVector, voxelMesh, plane, true, false);
        m_Ui.openGLWidget->setVoxelPositions(voxelPositions);
        m_Ui.openGLWidget->update();
    }
}

void MainWindow::renderSymmetryPlane(const bool render) const {
    m_Ui.openGLWidget->setRenderSymmetryPlane(render);
}

void MainWindow::railwayTopDownView() const {
    // Getting the basic data.
    const VoxelMesh& voxelMesh = m_LocalSymmetry.getVoxelMesh();
    const uint index = static_cast<uint>(m_Ui.sbx_ReflectionSymmetries->currentIndex());
    const Plane plane = m_Symmetries[index].getPlane();
    std::tuple<Vector3d, double, double, double> planeParams = plane.calculateStartPoint(voxelMesh);
    Vector3d start = std::get<0>(planeParams);
    double lengthX = std::get<2>(planeParams);
    double lengthY = std::get<3>(planeParams);

    // Moving camera to the symmetry.
    m_Ui.openGLWidget->resetCamera();
    if (plane.parallelVector().y < 0) {
        m_Ui.openGLWidget->moveCameraToPoint((start.x + lengthX) / 10.0, (start.y - lengthY) / 10.0, (-voxelMesh.maxZ - 150) / 10.0);
    }
    else {
        m_Ui.openGLWidget->moveCameraToPoint((start.x + lengthX) / 10.0, (start.y + lengthY) / 10.0, (-voxelMesh.maxZ - 150) / 10.0);
    }
    m_Ui.openGLWidget->rotateCamera(-90.0f, 0.0f, 0.0f);
    m_Ui.openGLWidget->update();
}