#include <fstream>
#include <glm/ext.hpp>
#include <iostream>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMessageBox>
#include <QMouseEvent>

#include "LocalSymmetryDetector/Structs/Color.hpp"
#include "LocalSymmetryDetector/Structs/Functions.hpp"
#include "LocalSymmetryDetector/Structs/Plane.tpp"
#include "LocalSymmetryDetector/Structs/Point.tpp"
#include "PlanePicker.hpp"


using namespace Symmetry;



void PlanePicker::compileShaders() {
    // Creating an object for OpenGL rendering.
    if (m_gl == nullptr) {
        return;
    }

    // Creating a new OpenGL instance.
    m_program = m_gl->glCreateProgram();

    // Vertex shader.
    {
        // Creating a vertex shader.
        GLint vShader = m_gl->glCreateShader(GL_VERTEX_SHADER);

        // Reading vertex shader.
        std::ifstream in("./OpenGL/vshader.vsh");
        std::stringstream buffer;
        buffer << in.rdbuf();
        std::string shaderContent = buffer.str();
        const GLchar* shaderContentArr = shaderContent.c_str();

        // Compiling and adding to the program.
        m_gl->glShaderSource(vShader, 1, &shaderContentArr, nullptr);
        m_gl->glCompileShader(vShader);
        m_gl->glAttachShader(m_program, vShader);
    }
    // Fragment shader.
    {
        // Creating a fragment shader.
        GLuint fShader = m_gl->glCreateShader(GL_FRAGMENT_SHADER);

        // Reading fragment shader.
        std::ifstream in("./OpenGL/fshader.fsh");
        std::stringstream buffer;
        buffer << in.rdbuf();
        std::string shaderContent = buffer.str();
        const GLchar* shaderContentArr = shaderContent.c_str();

        // Compiling and adding to the program.
        m_gl->glShaderSource(fShader, 1, &shaderContentArr, nullptr);
        m_gl->glCompileShader(fShader);
        m_gl->glAttachShader(m_program, fShader);
    }

    // Linking the object to the program.
    m_gl->glLinkProgram(m_program);
}

void PlanePicker::setPosition(const glm::mat4& P, const glm::mat4& V, const glm::mat4& M) {
    // Setting projection-view-model (PVM) matrix.
    GLint P_matrix = m_gl->glGetUniformLocation(m_program, "projection");  // Getting P matrix position in vertex shader.
    m_gl->glUniformMatrix4fv(P_matrix, 1, GL_FALSE, glm::value_ptr(P));     // Setting P matrix position in vertex shader.

    GLint V_matrix = m_gl->glGetUniformLocation(m_program, "view");        // Getting P matrix position in vertex shader.
    m_gl->glUniformMatrix4fv(V_matrix, 1, GL_FALSE, glm::value_ptr(V));     // Setting P matrix position in vertex shader.

    GLint M_matrix = m_gl->glGetUniformLocation(m_program, "model");       // Getting P matrix position in vertex shader.
    m_gl->glUniformMatrix4fv(M_matrix, 1, GL_FALSE, glm::value_ptr(M));     // Setting P matrix position in vertex shader.
}

Point Symmetry::PlanePicker::currentWorldPosition(const uint x, const uint y) {   
    glm::mat4 P = calculateProjectionMatrix(true);
    const glm::mat4 V = calculateViewMatrix();
    const int sign = x < static_cast<uint>(width() / 2) ? -1 : 1;
    const glm::vec3 position = glm::unProject(glm::vec3(x, y, 5000.0), V, P, glm::vec4(0.0, 0.0, width(), height()));
    const Point transposedPosition = Point(position.x + sign * (2 * m_aspectRatio), position.z + m_voxelMesh.deltaY, 0);

    return Point(transposedPosition.x, transposedPosition.y, 0);
}

glm::mat4 PlanePicker::calculateProjectionMatrix(const bool ratio) {
    // If the projection is set to orthogonal, the orthogonal matrix is returned.
    GLint viewport[4]{};
    m_gl->glGetIntegerv(GL_VIEWPORT, viewport);                                                                     // Getting width and height of the viewport.
    m_aspectRatio = static_cast<float>(viewport[2]) / viewport[3];                                                  // Calculating the ratio between the width and the height.
    const float width = 0.1 * (m_voxelMesh.deltaX > m_voxelMesh.deltaY ? m_voxelMesh.deltaX : m_voxelMesh.deltaY);  // Calculating the width of the projection.
    const float height = width;                                                                                     // Calculating the height of the projection.
    
    if (ratio) {
        return glm::ortho(-width * m_aspectRatio, width * m_aspectRatio, -height, height, -100.0f, 1000.0f);
    }
     
    return glm::ortho(-width, width, -height, height, -100.0f, 1000.0f);
}

glm::mat4 PlanePicker::calculateViewMatrix() const {
    glm::mat4 V = glm::mat4(1);
    V = glm::translate(V, glm::vec3(-m_voxelMesh.deltaX / 2, 0, m_voxelMesh.deltaY / 2));
    V = glm::inverse(V);
    V = glm::rotate_slow(V, glm::radians(0.0f), glm::vec3(0, 1, 0));
    V = glm::rotate_slow(V, glm::radians(-90.0f), glm::vec3(1, 0, 0));
    V = glm::inverse(V);

    return V;
}

void PlanePicker::displayError() const {
    // Potential error display.
    const GLuint err = m_gl->glGetError();
    if (err) {
        QMessageBox msgBox(
            QMessageBox::Icon::Critical,
            "Error",
            "OpenGL error: " + QString::number(err)
        );
        msgBox.exec();

        QApplication::exit(2);
    }
}

void PlanePicker::setColor(const LAS::Data::Vector4d& color) {
    GLint matrix = m_gl->glGetUniformLocation(m_program, "Color");  // Getting color matrix position in fragment shader.
    m_gl->glUniform4f(matrix, color.x, color.y, color.z, color.w);   // Setting color matrix in fragment shader.
}

void PlanePicker::paintPoints(const glm::mat4& P, const glm::mat4& V) {
    std::vector<glm::mat4> modelMatrices(m_points.size());

    // Setting the position of the point.
    setPosition(P, V);
    setColor(Color::black);

    // Drawing points.
    const double side = 0.3;  // Point side length given in meters.
    for (u128 i = 0; i < m_points.size(); i++) {
        const Point point = m_points[i];

        // Model matrix.
        glm::mat4 M = glm::translate(glm::mat4(1), glm::vec3(point.x - side / 2, point.y - side / 2, point.z - side / 2));
        M = glm::scale_slow(M, glm::vec3(side, side, side));
        modelMatrices[i] = M;

        setPosition(P, V, M);
        m_gl->glBindVertexArray(m_VAO_Cube);
        m_gl->glBindBuffer(GL_ARRAY_BUFFER, m_VBO_Cube);
        m_gl->glDrawArrays(GL_TRIANGLES, 0, 36);
    }

    //if (!m_points.empty()) {
    //    // Passing model matrices to the vertex shader.
    //    m_gl->glBindBuffer(GL_ARRAY_BUFFER, m_VBO_Instance);
    //    m_gl->glBufferData(GL_ARRAY_BUFFER, static_cast<qopengl_GLsizeiptr>(m_points.size() * sizeof(glm::mat4)), &modelMatrices[0], GL_STATIC_DRAW);
    //    m_gl->glEnableVertexAttribArray(3);
    //    m_gl->glEnableVertexAttribArray(4);
    //    m_gl->glEnableVertexAttribArray(5);
    //    m_gl->glEnableVertexAttribArray(6);
    //    m_gl->glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)0);
    //    m_gl->glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(1 * sizeof(glm::vec4)));
    //    m_gl->glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(2 * sizeof(glm::vec4)));
    //    m_gl->glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(3 * sizeof(glm::vec4)));
    //    m_gl->glVertexAttribDivisor(3, 1);
    //    m_gl->glVertexAttribDivisor(4, 1);
    //    m_gl->glVertexAttribDivisor(5, 1);
    //    m_gl->glVertexAttribDivisor(6, 1);

    //    // Instanced draw of the cubes.
    //    m_gl->glBindVertexArray(m_VAO_Cube);
    //    m_gl->glUniform1i(m_gl->glGetUniformLocation(m_program, "instanced"), true);
    //    m_gl->glDrawArraysInstanced(GL_TRIANGLES, 0, 36, static_cast<GLsizei>(m_points.size()));
    //    m_gl->glUniform1i(m_gl->glGetUniformLocation(m_program, "instanced"), false);
    //}
}

void PlanePicker::paintBoundingBox(const glm::mat4& P, const glm::mat4& V) {
    // Setting the color of the symmetry axis.
    setColor(Color::black);

    // Drawing lower bounding box line segment.
    glm::mat4 M = glm::mat4(1);
    M = glm::translate(M, glm::vec3(0, 0, 0));
    M = glm::scale(M, glm::vec3(m_voxelMesh.deltaX, 1, 1));
    // Setting the position of the symmetry axis.
    setPosition(P, V, M);
    m_gl->glBindVertexArray(m_VAO_Line);
    m_gl->glBindBuffer(GL_ARRAY_BUFFER, m_VBO_Line);
    m_gl->glDrawArrays(GL_LINES, 0, 2);


    // Drawing lower bounding box line segment.
    M = glm::mat4(1);
    M = glm::translate(M, glm::vec3(0, 0, 0 - m_voxelMesh.deltaY));
    M = glm::scale(M, glm::vec3(m_voxelMesh.deltaX, 1, 1));
    setPosition(P, V, M);
    m_gl->glBindVertexArray(m_VAO_Line);
    m_gl->glBindBuffer(GL_ARRAY_BUFFER, m_VBO_Line);
    m_gl->glDrawArrays(GL_LINES, 0, 2);

    // Drawing left bounding box line segment.
    M = glm::mat4(1);
    M = glm::translate(M, glm::vec3(0, 0, 0));
    M = glm::rotate(M, glm::radians(90.0f), glm::vec3(0, 1, 0));
    M = glm::scale(M, glm::vec3(m_voxelMesh.deltaY, 1, 1));
    setPosition(P, V, M);
    m_gl->glBindVertexArray(m_VAO_Line);
    m_gl->glBindBuffer(GL_ARRAY_BUFFER, m_VBO_Line);
    m_gl->glDrawArrays(GL_LINES, 0, 2);

    // Drawing right bounding box line segment.
    M = glm::mat4(1);
    M = glm::translate(M, glm::vec3(m_voxelMesh.deltaX, 0, 0));
    M = glm::rotate(M, glm::radians(90.0f), glm::vec3(0, 1, 0));
    M = glm::scale(M, glm::vec3(m_voxelMesh.deltaY, 1, 1));
    setPosition(P, V, M);
    m_gl->glBindVertexArray(m_VAO_Line);
    m_gl->glBindBuffer(GL_ARRAY_BUFFER, m_VBO_Line);
    m_gl->glDrawArrays(GL_LINES, 0, 2);
}

void PlanePicker::paintCurrentPlane(const glm::mat4& P, const glm::mat4& V) {
    //if (!m_firstClick) {
    //    return;
    //}
    
    // Setting the color of the symmetry axis.
    setColor(Color::red);

    glm::vec3 difference = m_currentMousePosition - glm::vec3(m_startPoint.x, m_startPoint.y, 0);
    if (m_secondClick) {
        difference = glm::vec3(m_endPoint.x, m_endPoint.y, m_endPoint.z) - glm::vec3(m_startPoint.x, m_startPoint.y, m_startPoint.z);
    }
    const double length = glm::length(difference);

    float angle = std::atan(difference.y / difference.x);
    if (difference.x < 0) {
        angle -= std::numbers::pi_v<float>;
    }

    glm::mat4 M = glm::mat4(1);
    M = glm::translate(M, glm::vec3(m_startPoint.x, m_voxelMesh.maxZ + 0.1, -m_startPoint.y));
    M = glm::rotate(M, angle, glm::vec3(0, 1, 0));
    M = glm::scale(M, glm::vec3(length, 1, 0.6));
    setPosition(P, V, M);
    m_gl->glBindVertexArray(m_VAO_Cube);
    m_gl->glBindBuffer(GL_ARRAY_BUFFER, m_VAO_Cube);
    m_gl->glDrawArrays(GL_TRIANGLES, 0, 36);
}



void PlanePicker::initializeGL() {
    // Loading OpenGL functions.
    m_gl = context()->extraFunctions();

    if (m_gl) {
        // Shader compliation.
        compileShaders();

        // Render parameters.
        m_gl->glEnable(GL_DEPTH_TEST);                            // Preventing triangle overlapping.
        m_gl->glDisable(GL_CULL_FACE);                            // Rendering both triangle faces.
        m_gl->glEnable(GL_MULTISAMPLE);                           // Enabling MSAA.
        m_gl->glEnable(GL_BLEND);                                 // Enabling blending (for transparency).
        m_gl->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  // Enabling alpha channel for transparency.

        // Creating a VAO for a cube.
        m_gl->glGenVertexArrays(1, &m_VAO_Cube);
        m_gl->glBindVertexArray(m_VAO_Cube);

        // Creating a VBO and a cube.
        glm::vec3 cube[] = {
            glm::vec3(0, 0, 1), glm::vec3(1, 0, 1), glm::vec3(1, 0, 0),   // Below.
            glm::vec3(0, 0, 0), glm::vec3(1, 0, 0), glm::vec3(0, 0, 1),   // Below.
            glm::vec3(1, 1, 1), glm::vec3(0, 1, 1), glm::vec3(1, 1, 0),   // Above.
            glm::vec3(1, 1, 0), glm::vec3(0, 1, 0), glm::vec3(0, 1, 1),   // Above.
            glm::vec3(0, 1, 1), glm::vec3(0, 0, 1), glm::vec3(1, 0, 1),   // Front.
            glm::vec3(1, 0, 1), glm::vec3(1, 1, 1), glm::vec3(0, 1, 1),   // Front.
            glm::vec3(1, 1, 0), glm::vec3(1, 0, 0), glm::vec3(0, 0, 0),   // Back.
            glm::vec3(0, 0, 0), glm::vec3(0, 1, 0), glm::vec3(1, 1, 0),   // Back.
            glm::vec3(0, 1, 1), glm::vec3(0, 1, 0), glm::vec3(0, 0, 0),   // Left.
            glm::vec3(0, 0, 1), glm::vec3(0, 0, 0), glm::vec3(0, 1, 1),   // Left.
            glm::vec3(1, 0, 1), glm::vec3(1, 0, 0), glm::vec3(1, 1, 1),   // Right.
            glm::vec3(1, 1, 0), glm::vec3(1, 1, 1), glm::vec3(1, 0, 0)    // Right.
        };
        m_gl->glGenBuffers(1, &m_VBO_Cube);
        m_gl->glBindBuffer(GL_ARRAY_BUFFER, m_VBO_Cube);
        m_gl->glBufferData(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
        m_gl->glEnableVertexAttribArray(0);
        m_gl->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

        m_gl->glGenBuffers(1, &m_VBO_Instance);
        m_gl->glBindBuffer(GL_ARRAY_BUFFER, m_VBO_Instance);

        m_gl->glGenVertexArrays(1, &m_VAO_Line);
        m_gl->glBindVertexArray(m_VAO_Line);

        // Creating a VBO.
        glm::vec3 line(1, 0, 0);
        m_gl->glGenBuffers(1, &m_VBO_Line);
        m_gl->glBindBuffer(GL_ARRAY_BUFFER, m_VBO_Line);
        m_gl->glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3), &(line)[0], GL_STATIC_DRAW);
        m_gl->glEnableVertexAttribArray(0);
        m_gl->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), 0);

        // Displaying a potential OpenGL error.
        displayError();

        m_gl->glUseProgram(m_program);
    }
    else {
        // If unable to load functions, an error is displayed.
        QMessageBox msgBox(
            QMessageBox::Icon::Critical,
            "Error",
            "Error loading OpenGL functions."
        );
        msgBox.exec();

        QApplication::exit(1);
    }
}

void PlanePicker::paintGL() {
    // Background and depth buffer cleaning.
    m_gl->glClearColor(1, 1, 1, 1);
    m_gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Calculating projection and view matrices.
    const glm::mat4 P = calculateProjectionMatrix(true);
    const glm::mat4 V = calculateViewMatrix();

    // Paint subprocedures.
    paintPoints(P, V);         // Painting points.
    paintBoundingBox(P, V);    // Painting the bounding box.
    paintCurrentPlane(P, V);

    // Displaying a potential OpenGL error.
    displayError();
    update();
}

void PlanePicker::resizeGL(int w, int h) {
    m_gl->glViewport(0, 0, w, h);  // Resizing an OpenGL viewport.
}

void PlanePicker::mousePressEvent(QMouseEvent* event) {
    // Getting mouse coordinates.
    const GLint x = event->pos().x();
    const GLint y = event->pos().y();

    const Point position = currentWorldPosition(x, y);


    //{
    //    std::stringstream ss;
    //    ss << position.x << ", " << position.y << ", " << position.z;
    //    QMessageBox(QMessageBox::Information, "TEST", QString::fromUtf8(ss.str())).exec();
    //}

    if (!m_firstClick) {
        m_startPoint = position;
        m_endPoint = position;
        m_currentMousePosition = glm::vec3(position.x, position.y, position.z);
        m_firstClick = true;
        m_secondClick = false;
    }
    else {
        m_endPoint = position;
        m_firstClick = false;
        m_secondClick = true;
        const Vector3d vector = (m_endPoint - m_startPoint).normalize();

        m_plane = Plane(m_startPoint, vector, Vector3d(0, 0, 1));

        //std::stringstream ss;
        //ss << plane.a << "x + " << plane.b << "y + " << plane.c << "z";
        //QMessageBox(QMessageBox::Information, "TEST", QString::fromUtf8(ss.str())).exec();
    }
}

void PlanePicker::mouseMoveEvent(QMouseEvent* event) {
    if (m_firstClick) {
        const GLint x = event->pos().x();
        const GLint y = event->pos().y();

        const Point position = currentWorldPosition(x, y);
        m_currentMousePosition = glm::vec3(position.x, position.y, 0);
    }
}




PlanePicker::PlanePicker(QWidget* parent) : QOpenGLWidget(parent) {
    setMouseTracking(true);  // Enabling mouse tracking (allowing camera rotation and move).
}

PlanePicker::~PlanePicker() {
    m_gl->glDeleteVertexArrays(1, &m_VAO_Cube);  // Deleting the cube VAO.
    m_gl->glDeleteProgram(m_program);            // Deleting the OpenGL instance.
}

void PlanePicker::setPoints(const std::vector<Point>& points) {
    // Clearing previous points.
    m_points = std::vector<Point>(points.size());

    // Since OpenGL has flipped Y and Z axes, Y and Z coordinates are changed.
    for (u128 i = 0; i < points.size(); i++) {
        const Point& point = points[i];
        m_points[i] = Point(point.x, point.z, -point.y);
    }

    // OpenGL rerendering.
    update();
}

void PlanePicker::setVoxelMesh(const VoxelMesh& voxelMesh) {
    m_voxelMesh = voxelMesh;
    update();
}

Plane PlanePicker::getPlane() const {
    return m_plane;
}