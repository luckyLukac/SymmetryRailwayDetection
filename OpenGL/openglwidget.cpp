#include <fstream>
#include <glm/gtx/normal.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <QApplication>
#include <QDir>
#include <QGuiApplication>
#include <QMessageBox>
#include <QMouseEvent>
#include <iostream>
#include <sstream>
#include <string>

#include "../Assets/ObjectLoader.hpp"
#include "LocalSymmetryDetector/Structs/Color.hpp"
#include "openglwidget.h"



// OBJECT METHODS
// Compiling and adding shaders to the program.
void OpenGLWidget::compileShaders() {
    // Creating an object for OpenGL rendering.
    if (gl == nullptr) {
        return;
    }

    // Creating a new OpenGL instance.
    id_Program = gl->glCreateProgram();

    // Vertex shader.
    {
        // Creating a vertex shader.
        GLint vShader = gl->glCreateShader(GL_VERTEX_SHADER);

        // Reading vertex shader.
        //std::ifstream in(PROJECT_PATH "OpenGL/vshader.vsh");
        std::ifstream in("./OpenGL/vshader.vsh");
        std::stringstream buffer;
        buffer << in.rdbuf();
        std::string shaderContent = buffer.str();
        const GLchar* shaderContentArr = shaderContent.c_str();

        // Compiling and adding to the program.
        gl->glShaderSource(vShader, 1, &shaderContentArr, nullptr);
        gl->glCompileShader(vShader);
        gl->glAttachShader(id_Program, vShader);
    }
    // Fragment shader.
    {
        // Creating a fragment shader.
        GLuint fShader = gl->glCreateShader(GL_FRAGMENT_SHADER);

        // Reading fragment shader.
        //std::ifstream in(PROJECT_PATH "OpenGL/fshader.fsh");
        std::ifstream in("./OpenGL/fshader.fsh");
        std::stringstream buffer;
        buffer << in.rdbuf();
        std::string shaderContent = buffer.str();
        const GLchar* shaderContentArr = shaderContent.c_str();

        // Compiling and adding to the program.
        gl->glShaderSource(fShader, 1, &shaderContentArr, nullptr);
        gl->glCompileShader(fShader);
        gl->glAttachShader(id_Program, fShader);
    }

    // Linking the object to the program.
    gl->glLinkProgram(id_Program);
}

// Calculating the projection matrix.
glm::mat4 OpenGLWidget::getProjectionMatrix() {
    // If the projection is set to perspective, the perspective matrix is returned.
    if (perspectiveProjection) {
        return glm::perspective(glm::radians(FOV), float(width()) / height(), 0.01f, 1000.0f);
    }

    // If the projection is set to orthogonal, the orthogonal matrix is returned.
    GLint viewport[4];
    gl->glGetIntegerv(GL_VIEWPORT, viewport);  // Getting width and height of the viewport.
    const float ratio = (float)viewport[2] / (float)viewport[3];  // Calculating the ratio between the width and the height.
    float width = 0.6 * (voxelMesh.deltaX > voxelMesh.deltaY ? voxelMesh.deltaX : voxelMesh.deltaY);  // Calculating the width of the projection.
    float height = width;  // Calculating the height of the projection.
    return glm::ortho(-width * ratio, width * ratio, -height, height, -100.0f, 1000.0f);
}

// Calculating the view matrix.
glm::mat4 OpenGLWidget::getViewMatrix() {
    glm::mat4 V = glm::mat4(1);
    V = glm::translate(V, glm::vec3(cameraPosition[0], cameraPosition[1], cameraPosition[2]));
    V = glm::inverse(V);
    V = glm::rotate_slow(V, glm::radians(cameraRotation[1]), glm::vec3(0, 1, 0));
    V = glm::rotate_slow(V, glm::radians(cameraRotation[0]), glm::vec3(1, 0, 0));
    V = glm::inverse(V);

    return V;
}

// Setting the position of an object.
void OpenGLWidget::setPosition(const glm::mat4& P, const glm::mat4& V, const glm::mat4& M) {
    // Setting projection-view-model (PVM) matrix.
    GLint P_matrix = gl->glGetUniformLocation(id_Program, "projection");  // Getting P matrix position in vertex shader.
    gl->glUniformMatrix4fv(P_matrix, 1, GL_FALSE, glm::value_ptr(P));     // Setting P matrix position in vertex shader.

    GLint V_matrix = gl->glGetUniformLocation(id_Program, "view");        // Getting P matrix position in vertex shader.
    gl->glUniformMatrix4fv(V_matrix, 1, GL_FALSE, glm::value_ptr(V));     // Setting P matrix position in vertex shader.

    GLint M_matrix = gl->glGetUniformLocation(id_Program, "model");       // Getting P matrix position in vertex shader.
    gl->glUniformMatrix4fv(M_matrix, 1, GL_FALSE, glm::value_ptr(M));     // Setting P matrix position in vertex shader.
}

// Setting a color of objects to be rendered.
void OpenGLWidget::setColor(const glm::vec4& color) {
    GLint matrix = gl->glGetUniformLocation(id_Program, "Color");  // Getting color matrix position in fragment shader.
    gl->glUniform4f(matrix, color.r, color.g, color.b, color.a);   // Setting color matrix in fragment shader.
}

// Setting a color of objects to be rendered.
void OpenGLWidget::setColor(const LAS::Data::Vector4d& color) {
    GLint matrix = gl->glGetUniformLocation(id_Program, "Color");  // Getting color matrix position in fragment shader.
    gl->glUniform4f(matrix, color.x, color.y, color.z, color.w);   // Setting color matrix in fragment shader.
}

// Bounding box render procedure.
void OpenGLWidget::paintBoundingBox(const glm::mat4& P, const glm::mat4& V) {
    // If painting bounding box is disabled, the function does nothing.
    if (!renderBoundingBox) {
        return;
    }

    // Creating the bounding box lines.
    glm::vec3 bbox[] = {
        glm::vec3(voxelMesh.minX, voxelMesh.minZ, voxelMesh.minY),
        glm::vec3(voxelMesh.maxX, voxelMesh.minZ, voxelMesh.minY),
        glm::vec3(voxelMesh.minX, voxelMesh.minZ, voxelMesh.maxY),
        glm::vec3(voxelMesh.minX, voxelMesh.minZ, voxelMesh.minY),
    };

    // Drawing 4 bounding box lines.
    for (int i = 0; i < 4; i++) {
        // Model matrix.
        glm::vec3 position = bbox[i];
        glm::mat4 M = glm::mat4(1);
        M = glm::translate(M, glm::vec3(position.x, position.y, -position.z));
        if (i % 2 == 0) {
            M = glm::scale(M, glm::vec3(voxelMesh.deltaX, 1, 1));
        }
        else {
            M = glm::rotate(M, glm::radians(90.0f), glm::vec3(0, 1, 0));
            M = glm::scale(M, glm::vec3(voxelMesh.deltaY, 1, 1));
        }

        // Setting the position of the line segment.
        setPosition(P, V, M);

        // Setting the color of the line segment.
        setColor(Color::black);

        // Drawing the line segment.
        gl->glBindVertexArray(id_VAO_cube);
        gl->glDrawArrays(GL_LINES, 3, 2);
    }
}

// Point render procedure.
void OpenGLWidget::paintPoints(const glm::mat4& P, const glm::mat4& V) {
    // If painting points is disabled, the function does nothing.
    if (!renderPoints) {
        return;
    }

    std::vector<glm::mat4> modelMatrices(points.size());
    std::vector<glm::vec4> colors(points.size());

    // Setting the position of the point.
    setPosition(P, V, glm::mat4(1));

    // Drawing points.
    const double side = 0.3;  // Point side length given in meters.
    for (u128 i = 0; i < points.size(); i++) {
        const auto& [point, pointPosition, isRendered] = points[i];

        // If the point should not be rendered, the
        // current iteration of the loop is skipped.
        if (!isRendered) {
            continue;
        }

        // If a certain layer is rendered and the point is not in
        // that layer, the current iteration of the loop is skipped.
        if (renderLayer != -1 && VoxelFunctions::calculateLayerInVoxelMesh(point.y, voxelMesh) != renderLayer) {
            continue;
        }

        // Model matrix.
        glm::mat4 M = glm::translate(glm::mat4(1), glm::vec3(point.x - side / 2, point.y - side / 2, point.z + side / 2));
        M = glm::scale_slow(M, glm::vec3(side, side, side));
        modelMatrices[i] = M;

        // Setting the color according to the position from the plane.
        if (pointPosition == PointPosition::center) {
            colors[i] = glm::vec4(Color::green.x, Color::green.y, Color::green.z, Color::lightGreen.w);  // Points in the center are light green.
        }
        else if (pointPosition == PointPosition::left) {
            colors[i] = glm::vec4(Color::red.x, Color::red.y, Color::red.z, Color::red.w);  // Points to the left are red.
        }
        else if (pointPosition == PointPosition::right) {
            colors[i] = glm::vec4(Color::blue.x, Color::blue.y, Color::blue.z, Color::blue.w);  // Points to the right are blue.
        }
        else if (pointPosition == PointPosition::rotational) {
            colors[i] = glm::vec4(Color::green.x, Color::green.y, Color::green.z, Color::magenta.w);  // Points in the rotational symmetry are magenta.
        }
        else if (pointPosition == PointPosition::railwayCandidate) {
            colors[i] = glm::vec4(1.0, 0.0, 0.0, 1.0);
        }
        else if (pointPosition == PointPosition::asymmetry) {
            colors[i] = glm::vec4(Color::red.x, Color::red.y, Color::red.z, 1.0);  // Points that are not in symmetry are gray.
        }
        else {
            colors[i] = glm::vec4(Color::gray.x, Color::gray.y, Color::gray.z, 1.0);  // Points that are not in symmetry are gray.
        }
    }

    if (!points.empty()) {
        // Passing color instances to the vertex shader.
        gl->glBindBuffer(GL_ARRAY_BUFFER, id_VBO_instance_color);
        gl->glBufferData(GL_ARRAY_BUFFER, static_cast<qopengl_GLsizeiptr>(points.size() * sizeof(glm::vec4)), &colors[0], GL_STATIC_DRAW);
        gl->glEnableVertexAttribArray(2);
        gl->glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
        gl->glVertexAttribDivisor(2, 1);

        // Passing model matrices to the vertex shader.
        gl->glBindBuffer(GL_ARRAY_BUFFER, id_VBO_instance);
        gl->glBufferData(GL_ARRAY_BUFFER, static_cast<qopengl_GLsizeiptr>(points.size() * sizeof(glm::mat4)), &modelMatrices[0], GL_STATIC_DRAW);
        gl->glEnableVertexAttribArray(3);
        gl->glEnableVertexAttribArray(4);
        gl->glEnableVertexAttribArray(5);
        gl->glEnableVertexAttribArray(6);
        gl->glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)0);
        gl->glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(1 * sizeof(glm::vec4)));
        gl->glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(2 * sizeof(glm::vec4)));
        gl->glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(3 * sizeof(glm::vec4)));
        gl->glVertexAttribDivisor(3, 1);
        gl->glVertexAttribDivisor(4, 1);
        gl->glVertexAttribDivisor(5, 1);
        gl->glVertexAttribDivisor(6, 1);

        // Instanced draw of the cubes.
        gl->glBindVertexArray(id_VAO_cube);
        gl->glUniform1i(gl->glGetUniformLocation(id_Program, "instanced"), true);
        gl->glDrawArraysInstanced(GL_TRIANGLES, 0, 36, static_cast<GLsizei>(points.size()));
        gl->glUniform1i(gl->glGetUniformLocation(id_Program, "instanced"), false);
    }
}

// Voxel render procedure.gl->glVertexAttribDivisor(3, 1);
void OpenGLWidget::paintVoxels(const glm::mat4& P, const glm::mat4& V) {
    if (renderVoxelsAsEdges) {
        paintVoxelsAsEdges(P, V);  // Painting voxels as edges.
    }
    else {
        paintVoxelsAsCubes(P, V);  // Painting voxels as cubes.
    }
}

// Voxel edge render procedure.
void OpenGLWidget::paintVoxelsAsEdges(const glm::mat4& P, const glm::mat4& V) {
    // If painting voxels is disabled, the function does nothing.
    if (!renderVoxels || voxelMesh.count <= 0) {
        return;
    }

    std::vector<glm::mat4> modelMatrices(voxelMesh.count);
    std::vector<glm::vec4> colors(voxelMesh.count);

    // Setting the position of the point.
    setPosition(P, V);

    const double edgeX = voxelMesh.sideX / 10.0;
    const double edgeY = voxelMesh.sideY / 10.0;
    const double edgeZ = voxelMesh.sideZ / 10.0;
    for (uint z = 0; z < voxelMesh.zSize; z++) {
        for (uint y = 0; y < voxelMesh.ySize; y++) {
            for (uint x = 0; x < voxelMesh.xSize; x++) {
                // Retrieving voxel parameters.
                const unsigned int i = z * voxelMesh.ySize * voxelMesh.xSize + y * voxelMesh.xSize + x;
                const Voxel& voxel = voxels[i];

                bool inSymmetry = voxel.inSymmetry;              // Checking whether the voxel is a part of currently displayed symmetry.
                bool superInteresting = voxel.interesting;  // Checking whether the voxel is super-interesting.
                bool interesting = voxel.material;            // Checking whether the voxel is interesting.

                // If the voxel is in symmetry and we hate
                // to paint those, we move on to the next one.
                if (inSymmetry && !renderSymmetryVoxels) {
                    continue;
                }

                // If the voxel is super-interesting and we hate
                // to paint those, we move on to the next one.
                if (!inSymmetry && superInteresting && !renderSuperInterestingVoxels) {
                    continue;
                }

                // If the voxel is interesting and we hate
                // to paint those, we move on to the next one.
                if (!inSymmetry && interesting && !superInteresting && !renderInterestingVoxels) {
                    continue;
                }

                // If the voxel is ordinary and we hate
                // to paint those, we move on to the next one.
                if (!inSymmetry && !interesting && !superInteresting && !renderOrdinaryVoxels) {
                    continue;
                }

                // If a certain layer is rendered and the voxel is not in
                // that layer, the current iteration of the loop is skipped.
                if (renderLayer != -1 && z != renderLayer) {
                    continue;
                }

                // Model matrix.
                glm::mat4 M = glm::mat4(1);
                M = glm::translate(M, glm::vec3(voxel.x / 10.0, voxel.y / 10.0, static_cast<int>(0 - voxelMesh.sideZ - voxel.z) / 10.0));
                M = glm::scale_slow(M, glm::vec3(edgeX, edgeY, edgeZ));
                modelMatrices[i] = M;

                colors[i] = glm::vec4(0.0, 0.0, 0.0, 1.0);
            }
        }
    }

    // Passing color instances to the vertex shader.
    gl->glBindBuffer(GL_ARRAY_BUFFER, id_VBO_instance_color);
    gl->glBufferData(GL_ARRAY_BUFFER, static_cast<qopengl_GLsizeiptr>(voxelMesh.count * sizeof(glm::vec4)), &colors[0], GL_STATIC_DRAW);
    gl->glEnableVertexAttribArray(2);
    gl->glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
    gl->glVertexAttribDivisor(2, 1);

    // Passing model matrices to the vertex shader.
    gl->glBindBuffer(GL_ARRAY_BUFFER, id_VBO_instance);
    gl->glBufferData(GL_ARRAY_BUFFER, static_cast<qopengl_GLsizeiptr>(voxelMesh.count * sizeof(glm::mat4)), &modelMatrices[0], GL_STATIC_DRAW);
    gl->glEnableVertexAttribArray(3);
    gl->glEnableVertexAttribArray(4);
    gl->glEnableVertexAttribArray(5);
    gl->glEnableVertexAttribArray(6);
    gl->glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)0);
    gl->glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(1 * sizeof(glm::vec4)));
    gl->glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(2 * sizeof(glm::vec4)));
    gl->glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(3 * sizeof(glm::vec4)));
    gl->glVertexAttribDivisor(3, 1);
    gl->glVertexAttribDivisor(4, 1);
    gl->glVertexAttribDivisor(5, 1);
    gl->glVertexAttribDivisor(6, 1);

    // Instanced draw of the cubes.
    gl->glBindVertexArray(id_VAO_cube);
    gl->glUniform1i(gl->glGetUniformLocation(id_Program, "instanced"), true);
    for (int i = 0; i < 12; i++) {
        gl->glDrawArraysInstanced(GL_LINES, 3 * i, 2, voxelMesh.count);
    }
    gl->glUniform1i(gl->glGetUniformLocation(id_Program, "instanced"), false);
}

// Voxel cube render procedure.
void OpenGLWidget::paintVoxelsAsCubes(const glm::mat4& P, const glm::mat4& V) {
    // If painting voxels is disabled, the function does nothing.
    if (!renderVoxels || voxelMesh.count == 0) {
        return;
    }

    std::vector<glm::mat4> modelMatrices(voxelMesh.count);
    std::vector<glm::vec4> colors(voxelMesh.count);

    // Setting the position of the point.
    setPosition(P, V, glm::mat4(1));

    const double edgeX = voxelMesh.sideX / 10.0;
    const double edgeY = voxelMesh.sideY / 10.0;
    const double edgeZ = voxelMesh.sideZ / 10.0;
    for (uint z = 0; z < voxelMesh.zSize; z++) {
        for (uint y = 0; y < voxelMesh.ySize; y++) {
            for (uint x = 0; x < voxelMesh.xSize; x++) {
                // Retrieving voxel parameters.
                const unsigned int i = z * voxelMesh.ySize * voxelMesh.xSize + y * voxelMesh.xSize + x;
                const Voxel& voxel = voxels[i];

                //bool rotational = voxelPositions[i] == PointPosition::rotational;  // Checking whether the voxel is in rotational symmetry.
                bool inAsymmetryAcross = voxel.inAsymmetryAcross;                                    // Checking whether the voxel is a part of currently displayed symmetry.
                bool inAsymmetry = voxel.inAsymmetry;                                    // Checking whether the voxel is a part of currently displayed symmetry.
                bool inSymmetry = voxel.inSymmetry;                                    // Checking whether the voxel is a part of currently displayed symmetry.
                bool bright = VoxelFunctions::isVoxelBright(x, y, z);               // Checking whether the voxel should be painted brightly.
                bool superInteresting = voxel.interesting;                        // Checking whether the voxel is super-interesting.
                bool interesting = voxel.material;                                  // Checking whether the voxel is interesting.
                //bool left = voxelPositions[i] == PointPosition::left;              // Checking whether the voxel is on the left side of the symmetry plane.
                //bool right = voxelPositions[i] == PointPosition::right;            // Checking whether the voxel is on the left side of the symmetry plane.

                // If the voxel is in symmetry and we hate
                // to paint those, we move on to the next one.
                if (inSymmetry && !renderSymmetryVoxels) {
                    continue;
                }

                // If the voxel is super-interesting and we hate
                // to paint those, we move on to the next one.
                if (!inSymmetry && superInteresting && !renderSuperInterestingVoxels) {
                    continue;
                }

                // If the voxel is interesting and we hate
                // to paint those, we move on to the next one.
                if (!inSymmetry && interesting && !superInteresting && !renderInterestingVoxels) {
                    continue;
                }

                // If the voxel is ordinary and we hate
                // to paint those, we move on to the next one.
                if (!inAsymmetryAcross && !inAsymmetry && !inSymmetry && !interesting && !superInteresting && !renderOrdinaryVoxels) {
                    continue;
                }

                // If a certain layer is rendered and the voxel is not in
                // that layer, the current iteration of the loop is skipped.
                if (renderLayer != -1 && VoxelFunctions::calculateLayerInVoxelMesh(voxel.y, voxelMesh) != renderLayer) {
                    continue;
                }

                // Model matrix.
                glm::mat4 M = glm::mat4(1);
                M = glm::translate(M, glm::vec3(voxel.x / 10.0, voxel.y / 10.0, static_cast<int>(-voxelMesh.sideZ /*- voxelMesh.sideZ*/ - voxel.z) / 10.0));
                M = glm::scale_slow(M, glm::vec3(edgeX, edgeZ, edgeY));
                modelMatrices[i] = M;

                // Setting the color according to voxel parameters.
                if (inAsymmetry) {
                    colors[i] = glm::vec4(Color::red.x, Color::red.y, Color::red.z, 0.6);
                }
                else if (inAsymmetryAcross) {
                    colors[i] = glm::vec4(Color::magenta.x, Color::magenta.y, Color::magenta.z, 0.6);
                }
                else if (inSymmetry && !interesting) {
                    colors[i] = glm::vec4(Color::blue.x, Color::blue.y, Color::blue.z, 0.6);
                }
                else if (inSymmetry) {
                    colors[i] = glm::vec4(Color::green.x, Color::green.y, Color::green.z, 0.6);
                }
                else if (superInteresting && bright) {
                    colors[i] = glm::vec4(Color::interestingVoxelLight.x, Color::interestingVoxelLight.y, Color::interestingVoxelLight.z, Color::interestingVoxelLight.w);
                }
                else if (superInteresting && !bright) {
                    colors[i] = glm::vec4(Color::interestingVoxelDark.x, Color::interestingVoxelDark.y, Color::interestingVoxelDark.z, Color::interestingVoxelDark.w);
                }
                else if (interesting && bright) {
                    colors[i] = glm::vec4(Color::interestingVoxelLight.x, Color::interestingVoxelLight.y, Color::interestingVoxelLight.z, Color::interestingVoxelLight.w);
                }
                else if (interesting && !bright) {
                    colors[i] = glm::vec4(Color::interestingVoxelDark.x, Color::interestingVoxelDark.y, Color::interestingVoxelDark.z, Color::interestingVoxelDark.w);
                }
                else if (bright) {
                    colors[i] = glm::vec4(Color::voxelLight.x, Color::voxelLight.y, Color::voxelLight.z, Color::voxelLight.w);
                }
                else if (!bright) {
                    colors[i] = glm::vec4(Color::voxelDark.x, Color::voxelDark.y, Color::voxelDark.z, Color::voxelDark.w);
                }
            }
        }
    }

    // Passing color instances to the vertex shader.
    gl->glBindBuffer(GL_ARRAY_BUFFER, id_VBO_instance_color);
    gl->glBufferData(GL_ARRAY_BUFFER, static_cast<qopengl_GLsizeiptr>(voxelMesh.count * sizeof(glm::vec4)), &colors[0], GL_STATIC_DRAW);
    gl->glEnableVertexAttribArray(2);
    gl->glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
    gl->glVertexAttribDivisor(2, 1);

    // Passing model matrices to the vertex shader.
    gl->glBindBuffer(GL_ARRAY_BUFFER, id_VBO_instance);
    gl->glBufferData(GL_ARRAY_BUFFER, static_cast<qopengl_GLsizeiptr>(voxelMesh.count * sizeof(glm::mat4)), &modelMatrices[0], GL_STATIC_DRAW);
    gl->glEnableVertexAttribArray(3);
    gl->glEnableVertexAttribArray(4);
    gl->glEnableVertexAttribArray(5);
    gl->glEnableVertexAttribArray(6);
    gl->glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)0);
    gl->glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(1 * sizeof(glm::vec4)));
    gl->glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(2 * sizeof(glm::vec4)));
    gl->glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(3 * sizeof(glm::vec4)));
    gl->glVertexAttribDivisor(3, 1);
    gl->glVertexAttribDivisor(4, 1);
    gl->glVertexAttribDivisor(5, 1);
    gl->glVertexAttribDivisor(6, 1);

    // Instanced draw of the cubes.
    gl->glBindVertexArray(id_VAO_cube);
    gl->glUniform1i(gl->glGetUniformLocation(id_Program, "instanced"), true);
    gl->glDrawArraysInstanced(GL_TRIANGLES, 0, 36, voxelMesh.count);
    gl->glUniform1i(gl->glGetUniformLocation(id_Program, "instanced"), false);
}

// Line segment render procedure.
void OpenGLWidget::paintLineSegments(const glm::mat4& P, const glm::mat4& V) {
    // If painting line segments is disabled, the function does nothing.
    if (!renderLineSegments) {
        return;
    }

    // Line segment rendering.
    for (unsigned long long i = 0; i < lineSegmentsPos.size(); i++) {
        // If a certain layer is rendered and the line segment is not in
        // that layer, the current iteration of the loop is skipped.
        if (renderLayer != -1 && VoxelFunctions::calculateLayerInVoxelMesh(lineSegmentsPos[i].y, voxelMesh) != renderLayer) {
            continue;
        }

        // Model matrix.
        glm::vec3 position = lineSegmentsPos[i];
        glm::mat4 M = glm::mat4(1);
        M = glm::translate(M, glm::vec3(position.x, position.y, -position.z));

        // Setting the position of the line segment.
        setPosition(P, V, M);

        // Setting the color of the line segment.
        setColor(Color::lineSegment);

        // Drawing the line segment.
        gl->glBindVertexArray(id_VAO_line);
        gl->glDrawArrays(GL_LINES, static_cast<GLint>(2 * i), 2);
    }
}

// Symmetry axis render procedure.
void OpenGLWidget::paintSymmetryAxis(const glm::mat4& P, const glm::mat4& V) {
    // If painting a symmetry axis is disabled, the function does nothing.
    if (!std::get<1>(symmetryAxis) || !renderSymmetryAxis) {
        return;
    }

    // Symmetry axis rendering.
    // Model matrix.
    glm::mat4 M = glm::mat4(1);
    glm::vec3 position = std::get<0>(symmetryAxis);
    M = glm::translate(M, glm::vec3(position.x, -100, -position.y));

    // Setting the position of the symmetry axis.
    setPosition(P, V, M);

    // Setting the color of the symmetry axis.
    setColor(Color::symmetryAxis);

    // Drawing the symmetry axis.
    gl->glBindVertexArray(id_VAO_axis);
    gl->glDrawArrays(GL_LINES, 0, 2);
}

// Symmetry plane render procedure.
void OpenGLWidget::paintSymmetryPlane(const glm::mat4& P, const glm::mat4& V) {
    // If painting a symmetry plane is disabled, the function does nothing.
    if (!renderSymmetryPlane) {
        return;
    }

    VoxelMesh voxelMesh = getVoxelMesh();
    const auto [position, distance, distanceByX, distanceByY] = plane.calculateStartPoint(voxelMesh);
    double planeY = -plane.parallelVector().y;

    // Model matrix.
    glm::mat4 M = glm::mat4(1);
    M = glm::translate(M, glm::vec3(position.x / 10.0, position.z / 10.0, -position.y / 10.0));
    M = glm::rotate(M, glm::orientedAngle(glm::vec3(0, 0, 1), glm::vec3(0, 0, planeY), glm::vec3(1, 0, 0)), glm::vec3(0, 1, 0));
    M = glm::scale_slow(M, glm::vec3(0.1 * cameraPosition[1] / 10.0, 1 * voxelMesh.deltaZ / 10.0, distance / 10.0));

    // Setting the position of the symmetry plane.
    setPosition(P, V, M);

    // Setting the color of the symmetry axis.
    setColor(LAS::Data::Vector4d(1.0, 0.4, 0.0, 1.0));

    // Nalaganje položaja kamere.
    glm::vec3 cameraPos = glm::vec3(cameraPosition[0], cameraPosition[1], cameraPosition[2]);
    gl->glUniform3fv(gl->glGetUniformLocation(id_Program, "cameraPosition"), 1, glm::value_ptr(cameraPos));
    glm::vec3 ambientInt = glm::vec3(0.0, 0.0, 0.0);
    gl->glUniform3fv(gl->glGetUniformLocation(id_Program, "ambientIntensity"), 1, glm::value_ptr(ambientInt));
    gl->glUniform1f(gl->glGetUniformLocation(id_Program, "ambientReflectivity"), 1.0);

    // Nalaganje položaja luči.
    glm::vec3 lightPos = glm::vec3(550, 0, -20);
    gl->glUniform3fv(gl->glGetUniformLocation(id_Program, "lightPosition"), 1, glm::value_ptr(lightPos));
    glm::vec3 lightInt = glm::vec3(0.5, 0.5, 0.5);
    gl->glUniform3fv(gl->glGetUniformLocation(id_Program, "lightIntensity"), 1, glm::value_ptr(lightInt));
    gl->glUniform1f(gl->glGetUniformLocation(id_Program, "diffuseReflectivity"), 0.8);
    gl->glUniform1f(gl->glGetUniformLocation(id_Program, "specularReflectivity"), 0.8);
    gl->glUniform1i(gl->glGetUniformLocation(id_Program, "exponent"), 32.0);
    gl->glUniform1f(gl->glGetUniformLocation(id_Program, "shading"), 1.0);

    // Drawing the symmetry plane.
    gl->glBindVertexArray(id_VAO_cube);
    gl->glDrawArrays(GL_TRIANGLES, 0, 36);

    gl->glUniform1f(gl->glGetUniformLocation(id_Program, "shading"), 0.0);
}

// Displaying a potential OpenGL error.
void OpenGLWidget::displayError() {
    // Potential error display.
    const GLuint err = gl->glGetError();
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



// SLOTS
// Camera translation or rotation.
void OpenGLWidget::mouseMoveEvent(QMouseEvent* event) {
    // COMMAND: Ctrl + Mid Mouse Button
    // Forward camera movement.
    if (QGuiApplication::queryKeyboardModifiers() == Qt::KeyboardModifier::ControlModifier &&
        event->buttons() == Qt::LeftButton)
    {
        const float SPEED = 40;
        const float z = event->pos().y();  // Z coordinate (flipped with Y).

        // Camera move, new Z is set.
        if (lastClick[2] == -1.f) {
            lastClick[2] = z;
            return;
        }

        // Forward camera movement.
        moveCameraForward(-SPEED * (z - lastClick[2]) / 100);

        // Previous coordinate update.
        lastClick[2] = -1.f;
    }
    // COMMAND: Shift + Mid Mouse Button
    // Side camera movement.
    else if (
        QGuiApplication::queryKeyboardModifiers() == Qt::KeyboardModifier::ShiftModifier &&
        event->buttons() == Qt::LeftButton)
    {
        const float SPEED = 20;
        const float x = event->pos().x();  // X coordinate.
        const float y = event->pos().y();  // Y coordinate.

        // Camera movement, new X and Y are set.
        if (lastClick[0] == -1.f) {
            lastClick[0] = x;
            lastClick[1] = y;
            return;
        }

        // Camera movement.
        moveCameraUp(SPEED * (y - lastClick[1]) / 100);      // Up/down camera movement.
        moveCameraToSide(SPEED * (x - lastClick[0]) / 100);  // Side camera movement.

        // Previous coordinates update.
        lastClick[0] = x;
        lastClick[1] = y;
    }
    // COMMAND: Mid Mouse Button
    // Camera rotation.
    else if (event->buttons() == Qt::LeftButton) {
        const float SPEED = 50;
        const float x = event->pos().x();  // X coordinate.
        const float y = event->pos().y();  // Y coordinate.

        // Camera rotation, new X and Y are set.
        if (lastClick[0] == -1.f) {
            lastClick[0] = x;
            lastClick[1] = y;
            return;
        }

        // Camera rotation.
        rotateCamera(SPEED * (lastClick[1] - y) / 100, SPEED * (lastClick[0] - x) / 100, 0);

        // Previous coordinates update.
        lastClick[0] = x;
        lastClick[1] = y;
    }

    // Coordinates reset.
    lastClick[0] = -1.f;
    lastClick[1] = -1.f;
    lastClick[2] = -1.f;
}

// Field of view change in the interval of [15°, 140°].
void OpenGLWidget::wheelEvent(QWheelEvent* event) {
    // If scrolling down, the Field of View (FOV) gets narrower.
    if (event->angleDelta().y() < 0) {
        // 140° is the maximum FOV.
        if (FOV > 140) {
            return;
        }

        // Increasement of the FOV by 2°.
        increaseFOV(2);
    }
    else {
        // 15° is the minimum FOV.
        if (FOV < 15) {
            return;
        }

        // Decreasement of the FOV by 2°.
        increaseFOV(-2);
    }
}



// OPENGL OBJECT METHODS
// OpenGL initialization procedure.
void OpenGLWidget::initializeGL() {
    // Loading OpenGL functions.
    gl = context()->extraFunctions();

    if (gl) {
        // Shader compliation.
        compileShaders();

        // Render parameters.
        gl->glEnable(GL_DEPTH_TEST);                            // Preventing triangle overlapping.
        gl->glDisable(GL_CULL_FACE);                            // Rendering both triangle faces.
        gl->glEnable(GL_MULTISAMPLE);                           // Enabling MSAA.
        gl->glEnable(GL_BLEND);                                 // Enabling blending (for transparency).
        gl->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  // Enabling alpha channel for transparency.

        // Creating a VAO for a cube.
        gl->glGenVertexArrays(1, &id_VAO_cube);
        gl->glBindVertexArray(id_VAO_cube);

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
        glm::vec3 normals[] = {
            glm::vec3(0, -1, 0), glm::vec3(0, -1, 0), glm::vec3(0, -1, 0),   // Below.
            glm::vec3(0, -1, 0), glm::vec3(0, -1, 0), glm::vec3(0, -1, 0),   // Below.
            glm::vec3(0,  1, 0), glm::vec3(0,  1, 0), glm::vec3(0,  1, 0),   // Above.
            glm::vec3(0,  1, 0), glm::vec3(0,  1, 0), glm::vec3(0,  1, 0),   // Above.
            glm::vec3(0,  0, 1), glm::vec3(0,  0, 1), glm::vec3(0,  0, 1),   // Front.
            glm::vec3(0,  0, 1), glm::vec3(0,  0, 1), glm::vec3(0,  0, 1),   // Front.
            glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1),   // Back.
            glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1),   // Back.
            glm::vec3(-1, 0, 0), glm::vec3(-1, 0, 0), glm::vec3(-1, 0, 0),   // Left.
            glm::vec3(-1, 0, 0), glm::vec3(-1, 0, 0), glm::vec3(-1, 0, 0),   // Left.
            glm::vec3(1, 0, 0), glm::vec3(1, 0, 0), glm::vec3(1, 0, 0),   // Right.
            glm::vec3(1, 0, 0), glm::vec3(1, 0, 0), glm::vec3(1, 0, 0),   // Right.
        };
        
        gl->glGenBuffers(1, &id_VBO_cube);
        gl->glBindBuffer(GL_ARRAY_BUFFER, id_VBO_cube);
        gl->glBufferData(GL_ARRAY_BUFFER, sizeof(cube), &cube, GL_STATIC_DRAW);
        gl->glEnableVertexAttribArray(0);
        gl->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

        gl->glGenBuffers(1, &id_cube_normals);
        gl->glBindBuffer(GL_ARRAY_BUFFER, id_cube_normals);
        gl->glBufferData(GL_ARRAY_BUFFER, sizeof(normals), &normals, GL_STATIC_DRAW);
        gl->glEnableVertexAttribArray(1);
        gl->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

        gl->glGenBuffers(1, &id_VBO_instance);
        gl->glBindBuffer(GL_ARRAY_BUFFER, id_VBO_instance);

        gl->glGenBuffers(1, &id_VBO_instance_color);
        gl->glBindBuffer(GL_ARRAY_BUFFER, id_VBO_instance_color);

        // Displaying a potential OpenGL error.
        displayError();

        gl->glUseProgram(id_Program);
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

// OpenGL render procedure.
void OpenGLWidget::paintGL() {
    // Background and depth buffer cleaning.
    gl->glClearColor(1, 1, 1, 1);
    gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Calculating projection and view matrices.
    glm::mat4 P = getProjectionMatrix();
    glm::mat4 V = getViewMatrix();

    // Paint subprocedures.
    paintPoints(P, V);         // Painting points.
    paintLineSegments(P, V);   // Painting line segments.
    paintSymmetryAxis(P, V);   // Painting a symmetry axis.
    paintSymmetryPlane(P, V);  // Painting a symmetry plane.
    paintBoundingBox(P, V);    // Painting the bounding box.
    paintVoxels(P, V);         // Painting voxels.

    // Displaying a potential OpenGL error.
    displayError();
}

// OpenGL resize procedure.
void OpenGLWidget::resizeGL(int w, int h) {
    gl->glViewport(0, 0, w, h);  // Resizing an OpenGL viewport.
}



// CONSTRUCTOR AND DESTRUCTOR
// Constructor of the widget.
OpenGLWidget::OpenGLWidget(QWidget* parent) : QOpenGLWidget(parent) {
    this->setMouseTracking(true);  // Enabling mouse tracking (allowing camera rotation and move).
}

// Destructor of the widget.
OpenGLWidget::~OpenGLWidget() {
    gl->glDeleteVertexArrays(1, &id_VAO_cube);  // Deleting the cube VAO.
    gl->glDeleteProgram(id_Program);            // Deleting the OpenGL instance.
}



// GETTERS AND SETTERS
// Setter for points to be rendered.
void OpenGLWidget::setPoints(const std::vector<Point>& pointVector, const std::vector<PointPosition>& positions) {
    // Clearing previous points.
    points.clear();

    // Since OpenGL has flipped Y and Z axes,
    // Y and Z coordinates are changed.
    for (u128 i = 0; i < pointVector.size(); i++) {
        // Adding the point to the list.
        points.push_back(
            std::make_tuple(
                Point(pointVector[i].x / 10, pointVector[i].z / 10, -pointVector[i].y / 10),
                positions[i],
                true
            )
        );
    }

    // OpenGL rerendering.
    update();
}

void OpenGLWidget::setPointPositions(const std::vector<Point>& points, const PointPosition& position) {
    for (const Point& point : points) {
        const Point openGLPoint(point.x / 10, point.z / 10, -point.y / 10);

        for (auto& pointTuple : this->points) {
            if (openGLPoint == std::get<0>(pointTuple)) {
                std::get<1>(pointTuple) = position;
            }
        }
    }

    // OpenGL rerendering.
    update();
}

// Voxel mesh setter.
void OpenGLWidget::setVoxelMesh(const VoxelMesh& vm) {
    this->voxelMesh = vm;

    update();
}

// Setter for voxels to be rendered.
void OpenGLWidget::setVoxels(
    const std::vector<std::vector<std::vector<Voxel>>>& voxels,
    const VoxelMesh& vm,
    const std::vector<PointPosition>& positions)
{
    this->voxels = std::vector<Voxel>();

    // Since OpenGL has flipped Y and Z axes,
    // Y and Z coordinates are changed.
    for (unsigned int i = 0; i < voxels.size(); i++) {
        for (unsigned int j = 0; j < voxels[i].size(); j++) {
            for (unsigned int k = 0; k < voxels[i][j].size(); k++) {
                Voxel v(
                    voxels[i][j][k].x,
                    voxels[i][j][k].z,
                    voxels[i][j][k].y,
                    voxels[i][j][k].interesting,
                    voxels[i][j][k].material
                );
                //if (v.inSymmetry) {
                v.inSymmetry = false;
                //}

                this->voxels.push_back(v);
                this->voxelPositions.push_back(positions[i]);
            }
        }
    }

    this->voxelMesh = vm;

    update();
}

// Setter for line segments to be rendered.
void OpenGLWidget::setLineSegments(const std::vector<LineSegment>& lineSegments, const Plane& Plane) {
    makeCurrent();

    this->lineSegmentsVectors.clear();
    this->lineSegmentsPos.clear();

    // Since OpenGL has flipped Y and Z axes,
    // Y and Z coordinates are changed.
    for (unsigned int i = 0; i < lineSegments.size(); i++) {
        Point p1 = lineSegments[i].p1;
        Point p2 = lineSegments[i].p2;

        this->lineSegmentsVectors.push_back(glm::vec3(0, 0, 0));
        this->lineSegmentsVectors.push_back(
            glm::vec3(
                (p2.x - p1.x) / 10.0,
                (p2.z - p1.z) / 10.0,
                (p1.y - p2.y) / 10.0
            )
        );
        this->lineSegmentsPos.push_back(glm::vec3(p1.x / 10.0, p1.z / 10.0, p1.y / 10.0));
    }

    gl->glGenVertexArrays(1,& id_VAO_line);
    gl->glBindVertexArray(id_VAO_line);

    // Creating a VBO.
    gl->glGenBuffers(1,& id_VBO_line);
    gl->glBindBuffer(GL_ARRAY_BUFFER, id_VBO_line);
    gl->glBufferData(GL_ARRAY_BUFFER, static_cast<qopengl_GLsizeiptr>(this->lineSegmentsVectors.size() * sizeof(glm::vec3)), &(this->lineSegmentsVectors)[0], GL_STATIC_DRAW);
    gl->glEnableVertexAttribArray(0);
    gl->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), 0);

    // Potential error display.
    const GLuint err = gl->glGetError();
    if (err) {
        QMessageBox msgBox(
            QMessageBox::Icon::Critical,
            "Error",
            "OpenGL error: " + QString::number(err)
        );
        msgBox.show();

        QApplication::exit(2);
    }

    this->plane = Plane;

    update();
}

// Setter for symmetry axis.
void OpenGLWidget::setSymmetryAxis(const Point& axis, bool active) {
    makeCurrent();

    // Creating a symmetry axis.
    std::vector<glm::vec3> symmetryAxisRender {glm::vec3(0, 0, 0), glm::vec3(0, 2000, 0)};
    symmetryAxis = std::make_pair(glm::vec3(axis.x, axis.y, axis.z), active);

    gl->glGenVertexArrays(1, &id_VAO_axis);
    gl->glBindVertexArray(id_VAO_axis);

    // Creating a VBO.
    gl->glGenBuffers(1, &id_VBO_axis);
    gl->glBindBuffer(GL_ARRAY_BUFFER, id_VBO_axis);
    gl->glBufferData(GL_ARRAY_BUFFER, static_cast<qopengl_GLsizeiptr>(symmetryAxisRender.size() * sizeof(glm::vec3)), &(symmetryAxisRender)[0], GL_STATIC_DRAW);
    gl->glEnableVertexAttribArray(0);
    gl->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), 0);

    // Potential error display.
    const GLuint err = gl->glGetError();
    if (err) {
        QMessageBox msgBox(
            QMessageBox::Icon::Critical,
            "Error",
            "OpenGL error: " + QString::number(err)
        );
        msgBox.show();

        QApplication::exit(2);
    }

    update();
}

// Voxels getter.
std::vector<Voxel>& OpenGLWidget::getVoxels() {
    return voxels;
}

// Voxel mesh getter.
VoxelMesh& OpenGLWidget::getVoxelMesh() {
    return this->voxelMesh;
}

// Perspective projection boolean setter.
void OpenGLWidget::setRenderPerspectiveProjection(const bool b) {
    perspectiveProjection = b;

    // OpenGL rerendering.
    update();
}

// Bounding box boolean setter.
void OpenGLWidget::setRenderBoundingBox(const bool b) {
    renderBoundingBox = b;

    // OpenGL rerendering.
    update();
}

// Point rendering setter.
void OpenGLWidget::setRenderPoints(const bool b) {
    makeCurrent();
    renderPoints = b;
    update();
}

// Voxel rendering setter.
void OpenGLWidget::setRenderVoxels(const bool b) {
    makeCurrent();
    renderVoxels = b;
    update();
}

// Voxel rendering as edges setter.
void OpenGLWidget::setRenderVoxelsAsEdges(const bool b) {
    renderVoxelsAsEdges = b;

    // OpenGL rerendering.
    update();
}

// Ordinary voxels rendering setter.
void OpenGLWidget::setRenderOrdinaryVoxels(const bool b) {
    makeCurrent();
    renderOrdinaryVoxels = b;
    update();
}

// Super-interesting voxels rendering setter.
void OpenGLWidget::setRenderSuperInterestingVoxels(const bool b) {
    makeCurrent();
    renderSuperInterestingVoxels = b;
    update();
}

// Interesting voxels rendering setter.
void OpenGLWidget::setRenderInterestingVoxels(const bool b) {
    makeCurrent();
    renderInterestingVoxels = b;
    update();
}

// Symmetry voxels rendering setter.
void OpenGLWidget::setRenderSymmetryVoxels(const bool b) {
    makeCurrent();
    renderSymmetryVoxels = b;
    update();
}

// Symmetry plane rendering setter.
void OpenGLWidget::setRenderSymmetryPlane(const bool b) {
    makeCurrent();
    renderSymmetryPlane = b;
    update();
}

// Line segment rendering setter.
void OpenGLWidget::setRenderLineSegments(const bool b) {
    makeCurrent();
    renderLineSegments = b;
    update();
}

// Symmetry axis render procedure.
void OpenGLWidget::setRenderSymmetryAxis(const bool b) {
    makeCurrent();
    renderSymmetryAxis = b;
    update();
}

// Layer rendering setter.
void OpenGLWidget::setRenderLayer(const int layer) {
    renderLayer = layer;
}

// Setting voxel positions according to a symmetry plane.
void OpenGLWidget::setVoxelPositions(const std::vector<PointPosition>& positions) {
    this->voxelPositions = positions;
}

// Getting voxel positions.
std::vector<PointPosition>& OpenGLWidget::getVoxelPositions() {
    return voxelPositions;
}



// OBJECT METHODS
// Getting camera parameters.
Camera OpenGLWidget::getCameraParameters() {
    // Preparing data structures.
    const std::vector<double> positions({ cameraPosition[0], cameraPosition[1], cameraPosition[2] });  // Getting camera positions.
    const std::vector<double> rotations({ cameraRotation[0], cameraRotation[1], cameraRotation[2] });  // Getting camera rotations.

    return Camera(FOV, positions, rotations);
}

// Setting camera Field of View (FOV).
void OpenGLWidget::setFOV(const float FOV) {
    this->FOV = FOV;

    // OpenGL rerendering.
    update();
}

// Camera Field of View (FOV) increasement.
void OpenGLWidget::increaseFOV(const float delta) {
    FOV += delta;

    // OpenGL rerendering.
    update();
}

// Forward camera movement.
void OpenGLWidget::moveCameraForward(const float forward) {
    // Getting yaw and pitch.
    const float yaw = cameraRotation[1];
    const float pitch = cameraRotation[0];

    // Calculating coordinates for the forward (backward) movement.
    cameraPosition[0] += forward * sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    cameraPosition[1] -= forward * sin(glm::radians(pitch));
    cameraPosition[2] += forward * cos(glm::radians(yaw)) * cos(glm::radians(pitch));

    // OpenGL rerendering.
    update();
}

// Upward camera movement.
void OpenGLWidget::moveCameraUp(const float up) {
    // Getting yaw and pitch.
    const float yaw = cameraRotation[1];
    const float pitch = cameraRotation[0];

    // Calculating coordinates for the upward (downward) movement.
    cameraPosition[0] -= up * sin(glm::radians(pitch)) * sin(glm::radians(yaw));
    cameraPosition[1] -= up * cos(glm::radians(pitch));
    cameraPosition[2] -= up * sin(glm::radians(pitch)) * cos(glm::radians(yaw));

    // OpenGL rerendering.
    update();
}

// Side camera movement.
void OpenGLWidget::moveCameraToSide(const float side) {
    // Getting yaw.
    const float yaw = cameraRotation[1];

    // Calculating coordinates for the side movement.
    cameraPosition[0] += side * cos(glm::radians(yaw));
    cameraPosition[2] -= side * sin(glm::radians(yaw));

    // OpenGL rerendering.
    update();
}

// Camera movement to a certain point.
void OpenGLWidget::moveCameraToPoint(const float x, const float y, const float z) {
    // Setting point coordinates.
    cameraPosition[0] = -x;
    cameraPosition[1] = z;
    cameraPosition[2] = y;

    // OpenGL rerendering.
    update();
}

// Camera rotation.
void OpenGLWidget::rotateCamera(const float x, const float y, const float z) {
    // Setting a new camera rotation.
    cameraRotation[0] += x;
    cameraRotation[1] += y;
    cameraRotation[2] += z;

    // OpenGL rerendering.
    update();
}

// Camera reset to initial values.
void OpenGLWidget::resetCamera() {
    // Resetting camera parameters to initial values.
    FOV = 60.0;
    cameraPosition[0] = cameraPosition[1] = 0.0;
    cameraPosition[2] = -15.0;
    cameraRotation[0] = cameraRotation[1] = cameraRotation[2] = 0.0;

    // OpenGL rerendering.
    update();
}