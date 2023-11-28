#pragma once

#include <LocalSymmetryDetector/LocalSymmetry.hpp>
#include <LocalSymmetryDetector/ReflectionSymmetry.hpp>
#include <QtWidgets/QMainWindow>
#include <vector>

#include "ui_MainWindow.h"


class MainWindow : public QMainWindow {
    Q_OBJECT
private:
    Ui::MainWindow m_Ui;
    std::vector<Point> m_Points;
    std::vector<Point> m_FilteredPoints;
    ReflectionSymmetry m_LocalSymmetry;
    std::vector<ReflectionSymmetry> m_Symmetries;


// PRIVATE FUNCTIONS
private:
    /// <summary>
    /// Function for printing the text in Qt console.
    /// </summary>
    /// <param name="title">: text title</param>
    /// <param name="body">: text body</param>
    void print(const std::string& title, const std::string& body) const;

    /// <summary>
    /// Stringifying basic point cloud data.
    /// </summary>
    /// <param name="pointCount">: number of points</param>
    /// <param name="voxelMesh">: voxel mesh</param>
    /// <returns>Stringified point cloud data</returns>
    std::string pointCloudDataToString(const u128 pointCount, const VoxelMesh& voxelMesh) const;

    /// <summary>
    /// Stringifying Railway Detection Algorithm (RDA).
    /// <param name="time">: elapsed time in milliseconds</param>
    /// <param name="numberOfSymmetries">: number of detected symmetries</param>
    /// <param name="voxelMesh">: voxel mesh</param>
    /// </summary>
    std::string railwayDetectionString(const u128 time, const uint numberOfSymmetries, const VoxelMesh& voxelMesh) const;

// SLOTS
private slots:
    /// <summary>
    /// Displaying the basic information about the app.
    /// </summary>
    void displayBasicInfo() const;

    /// <summary>
    /// Resetting camera position and rotation to initial values.
    /// </summary>
    void resetCamera() const;

    /// <summary>
    /// Setting camera view to top-down.
    /// </summary>
    void topDownView() const;

    /// <summary>
    /// Opening of the window with camera parameters.
    /// </summary>
    void openCameraParameters() const;

    /// <summary>
    /// Slot for finding a point cloud.
    /// </summary>
    void findPointCloud();

    /// <summary>
    /// Loading a chosen point cloud.
    /// </summary>
    void loadPointCloud();

    /// <summary>
    /// Applying the geometric point filter on the point cloud.
    /// </summary>
    void applyGeometricPointFilter();

    /// <summary>
    /// Function for railway detection.
    /// </summary>
    void detectRailways();


    /// <summary>
    /// Perspective projection setter.
    /// </summary>
    void setPerspectiveProjection() const;

    /// <summary>
    /// Orthogonal projection setter.
    /// </summary>
    void setOrthogonalProjection() const;

    /// <summary>
    /// Points rendering setter.
    /// <param name="render">: rendering of points if true</param>
    /// </summary>
    void renderPoints(const bool render) const;

    /// <summary>
    /// Voxel rendering setter (as edge).
    /// </summary>
    void renderVoxelsAsEdges() const;

    /// <summary>
    /// Voxel rendering setter (as cube).
    /// </summary>
    void renderVoxelsAsCubes() const;

    /// <summary>
    /// Voxel rendering setter.
    /// </summary>
    /// <param name="render">: rendering of voxels if true</param>
    void renderVoxels(const bool render) const;

    /// <summary>
    /// Symmetry voxel rendering setter.
    /// </summary>
    /// <param name="render">: rendering of symmetry voxels if true</param>
    void renderSymmetryVoxels(const bool render) const;

    /// <summary>
    /// Interesting voxel rendering setter.
    /// </summary>
    /// <param name="render">: rendering of interesting voxels if true</param>
    void renderInterestingVoxels(const bool render) const;

    /// <summary>
    /// Material voxel rendering setter.
    /// </summary>
    /// <param name="render">: rendering of material voxels if true</param>
    void renderMaterialVoxels(const bool render) const;

    /// <summary>
    /// Ordinary voxel rendering setter.
    /// </summary>
    /// <param name="render">: rendering of ordinary voxels if true</param>
    void renderOrdinaryVoxels(const bool render) const;

    /// <summary>
    /// Rendering of the chosen symmetry.
    /// </summary>
    /// <param name="index">: index of the symmetry</param>
    void renderSymmetry(const int index) const;

    /// <summary>
    /// Symmetry plane rendering setter.
    /// </summary>
    /// <param name="render">: rendering of symmetry plane if true</param>
    void renderSymmetryPlane(const bool render) const;

    /// <summary>
    /// Top-down view of the railway.
    /// </summary>
    void railwayTopDownView() const;

public:
    MainWindow(QWidget* parent = nullptr);
    ~MainWindow();
};
