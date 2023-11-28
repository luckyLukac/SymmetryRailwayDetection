#ifndef OPENGLWIDGET_H
#define OPENGLWIDGET_H

#include <glm/ext.hpp>
#include <glm/glm.hpp>
#include <LAS/Data/nVector.h>
#include <QOpenGLExtraFunctions>
#include <QtOpenGLWidgets/QOpenGLWidget>
#include <tuple>
#include <vector>

#include "Camera/CameraParameters.hpp"
#include "LocalSymmetryDetector/Structs/LineSegment.tpp"
#include "LocalSymmetryDetector/Structs/Plane.tpp"
#include "LocalSymmetryDetector/Structs/Point.tpp"
#include "LocalSymmetryDetector/Structs/Voxel.hpp"

using namespace Symmetry;


// Widget class where OpenGL is rendered.
class OpenGLWidget : public QOpenGLWidget {
private:
    // OBJECT VARIABLES
    // OpenGL variables
    QOpenGLExtraFunctions* gl = nullptr;          // OpenGL functions.
    GLuint id_Program = 0;                        // OpenGL render object.
    GLuint id_VAO_cube = 0;                       // Vertex Array Object for cube (additional vertex data).
    GLuint id_VBO_cube = 0;                       // Vertex Buffer Object for cube (basic vertex data).
    GLuint id_cube_normals = 0;                       // Vertex Buffer Object for cube (basic vertex data).
    GLuint id_VAO_line = 0;                       // Vertex Array Object for line (additional vertex data).
    GLuint id_VBO_line = 0;                       // Vertex Buffer Object for line (basic vertex data).
    GLuint id_VAO_axis = 0;                       // Vertex Array Object for symmetry axis (additional vertex data).
    GLuint id_VBO_axis = 0;                       // Vertex Buffer Object for symmetry axis (basic vertex data).
    GLuint id_VBO_instance = 0;
    GLuint id_VBO_instance_color = 0;

    // Camera variables.
    float FOV = 60.0;                             // Field of view (in degrees).
    float cameraPosition[3] = {0.0, 0.0, -15.0};  // Camera position coordinates (X, Y, Z).
    float cameraRotation[3] = {0.0, 0.0, 0.0};    // Camera rotation in degrees (X, Y, Z).

    // Render variables.
    std::vector<std::tuple<Point, PointPosition, bool>> points;                        // Tuples of points (from a LAS file), positions from the plane and rendering.
    std::vector<Voxel> voxels;                                                         // Voxels from a voxelized scene.
    std::vector<PointPosition> voxelPositions;                                         // Voxel positions according to a symmetry plane.
    std::pair<glm::vec3, bool> symmetryAxis = std::make_pair(glm::vec3(), false);      // Symmetry axis to be rendered and its validity.
    std::vector<glm::vec3> lineSegmentsVectors;                                        // Line segments to be rendered.
    std::vector<glm::vec3> lineSegmentsPos;                                            // Line segments to be rendered.
    Plane plane;                                                                       // Plane position.
    VoxelMesh voxelMesh;                                                               // Voxel mesh for data about voxels.

    // Render booleans.
    bool perspectiveProjection = true;                                                 // Rendering in perspective projection if true.
    bool renderBoundingBox = false;                                                    // Rendering bounding box of points.
    bool renderPoints = true;                                                          // Rendering all points if true.
    bool renderVoxels = true;                                                          // Rendering voxels if true.
    bool renderVoxelsAsEdges = false;                                                   // Rendering voxels as edges if true.
    bool renderOrdinaryVoxels = false;                                                  // Rendering ordinary voxels if true.
    bool renderSuperInterestingVoxels = false;                                          // Rendering super-interesting voxels if true.
    bool renderInterestingVoxels = false;                                               // Rendering interesting voxels if true.
    bool renderSymmetryVoxels = true;                                                  // Rendering symmetry voxels if true.
    bool renderLineSegments = false;                                                   // Rendering line segments if true.
    bool renderSymmetryAxis = false;                                                    // Rendering symmetry axis if true.
    bool renderSymmetryPlane = false;                                                  // Rendering symmetry plane if true.
    int renderLayer = -1;                                                              // Index of the layer to be rendered (all layers if -1).


    // OBJECT METHODS
    void compileShaders();                                                                         // Compiling and adding shaders to the program.
    glm::mat4 getProjectionMatrix();                                                               // Calculating the projection matrix.
    glm::mat4 getViewMatrix();                                                                     // Calculating the view matrix.
    void setPosition(const glm::mat4& P, const glm::mat4& V, const glm::mat4& M = glm::mat4(1));   // Setting the position of an object (Projection-View-Model).
    void setColor(const glm::vec4& color);                                                         // Setting a color of objects to be rendered.
    void setColor(const LAS::Data::Vector4d& color);                                               // Setting a color of objects to be rendered.
    void paintBoundingBox(const glm::mat4& P, const glm::mat4& V);                                 // Bounding box render procedure.
    void paintPoints(const glm::mat4& P, const glm::mat4& V);                                      // Point render procedure.
    void paintVoxels(const glm::mat4& P, const glm::mat4& V);                                      // Voxel render procedure.
    void paintVoxelsAsEdges(const glm::mat4& P, const glm::mat4& V);                               // Voxel edge render procedure.
    void paintVoxelsAsCubes(const glm::mat4& P, const glm::mat4& V);                               // Voxel cube render procedure.
    void paintLineSegments(const glm::mat4& P, const glm::mat4& V);                                // Line segment render procedure.
    void paintSymmetryAxis(const glm::mat4& P, const glm::mat4& V);                                // Symmetry axis render procedure.
    void paintSymmetryPlane(const glm::mat4& P, const glm::mat4& V);                               // Symmetry plane render procedure.
    void displayError();                                                                           // Displaying a potential OpenGL error.

    // SLOTS AND HELPER VARIABLES
    float lastClick[3] = {-1.0, -1.0, -1.0};           // Last click position on the widget.
    void mouseMoveEvent(QMouseEvent* event) override;  // Camera translation or rotation.
    void wheelEvent(QWheelEvent* event) override;      // Field of view change in the interval of [15°, 140°].


protected:
    // OPENGL OBJECT METHODS
    void initializeGL() override;          // OpenGL initialization procedure.
    void paintGL() override;               // OpenGL render procedure.
    void resizeGL(int w, int h) override;  // OpenGL resize procedure.

public:
    // CONSTRUCTOR AND DESTRUCTOR
    OpenGLWidget(QWidget* parent);  // Constructor of the widget.
    ~OpenGLWidget();                // Destructor of the widget.

    // GETTERS AND SETTERS
    void setPoints(const std::vector<Point>& pointVector, const std::vector<PointPosition>& positions);  // Setter for points to be rendered.
    void setPointPositions(const std::vector<Point>& points, const PointPosition& position);
    void setVoxelMesh(const VoxelMesh& vm);                                           // Voxel mesh setter.
    void setVoxels(const std::vector<std::vector<std::vector<Voxel>>>& voxels, const VoxelMesh& vm, const std::vector<PointPosition>& positions);  // Setter for voxels to be rendered.
    void setLineSegments(const std::vector<LineSegment>& lineSegments, const Plane& plane);  // Setter for line segments to be rendered.
    void setSymmetryAxis(const Point& axis, bool active);                             // Setter for symmetry axis.
    std::vector<Voxel>& getVoxels();                                                  // Voxels getter.
    VoxelMesh& getVoxelMesh();                                                        // Voxel mesh getter.
    void setRenderPerspectiveProjection(const bool b);                                // Perspective projection boolean setter.
    void setRenderBoundingBox(const bool b);                                          // Bounding box boolean setter.
    void setRenderPoints(const bool b);                                               // Point rendering setter.
    void setRenderVoxels(const bool b);                                               // Voxel rendering setter.
    void setRenderVoxelsAsEdges(const bool b);                                        // Voxel rendering as edges setter.
    void setRenderOrdinaryVoxels(const bool b);                                       // Ordinary voxels rendering setter.
    void setRenderSuperInterestingVoxels(const bool b);                               // Super-interesting voxels rendering setter.
    void setRenderInterestingVoxels(const bool b);                                    // Interesting voxels rendering setter.
    void setRenderSymmetryVoxels(const bool b);                                       // Symmetry voxels rendering setter.
    void setRenderLineSegments(const bool b);                                         // Line segment rendering setter.
    void setRenderSymmetryAxis(const bool b);                                         // Symmetry axis render procedure.
    void setRenderSymmetryPlane(const bool b);                                        // Symmetry plane rendering setter.
    void setRenderLayer(const int layer);                                             // Layer rendering setter.
    void setVoxelPositions(const std::vector<PointPosition>& positions);          // Setting voxel positions according to a symmetry plane.
    std::vector<PointPosition>& getVoxelPositions();                              // Getting voxel positions.

    // OBJECT METHODS
    Camera getCameraParameters();                                                     // Getting camera parameters.
    void setFOV(const float FOV);                                                     // Setting camera Field of View (FOV).
    void increaseFOV(const float delta);                                              // Camera Field of View (FOV) increasement.
    void moveCameraForward(const float forward);                                      // Forward camera movement.
    void moveCameraUp(const float up);                                                // Upward camera movement.
    void moveCameraToSide(const float side);                                          // Side camera movement.
    void moveCameraToPoint(const float x, const float y, const float z);              // Camera movement to a certain point.
    void rotateCamera(const float x, const float y, const float z);                   // Camera rotation.
    void resetCamera();                                                               // Camera reset to initial values.
};

#endif // OPENGLWIDGET_H
