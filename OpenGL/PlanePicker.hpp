#pragma once

#include <glm/glm.hpp>
#include <LAS/Data/nVector.h>
#include <QOpenGLExtraFunctions>
#include <QtOpenGLWidgets/QOpenGLWidget>

#include "LocalSymmetryDetector/Structs/Plane.tpp"
#include "LocalSymmetryDetector/Structs/Voxel.hpp"


namespace Symmetry {
	/// <summary>
	/// Class for picking a symmetry plane to be displayed in OpenGLWidget.
	/// </summary>
	class PlanePicker : public QOpenGLWidget {
	private:
		QOpenGLExtraFunctions* m_gl = nullptr;  // OpenGL functions.
		GLuint m_program = 0;                   // OpenGL program handle.
		GLuint m_VAO_Cube = 0;                  // Vertex Array Object for a cube (additional vertex data).
		GLuint m_VBO_Cube = 0;                  // Vertex Buffer Object for a cube (basic vertex data).
		GLuint m_VBO_Instance = 0;              // Vertex Buffer Object for the instanced painting of points.
		GLuint m_VAO_Line = 0;					// Vertex Array Object for a line.
		GLuint m_VBO_Line = 0;					// Vertex Buffer Object for a line.
		std::vector<Point> m_points;			// Points to be drawn.
		VoxelMesh m_voxelMesh;                  // Voxel mesh for data about voxels.
		bool m_firstClick = false;				// Drawing line if true, false otherwise.
		bool m_secondClick = false;				// Drawing line if true, false otherwise.
		Point m_startPoint;						// Start point of the line.
		Point m_endPoint;						// End point of the line.
		glm::vec3 m_currentMousePosition = glm::vec3(0, 0, 0);
		Plane m_plane;
		float m_aspectRatio = 1.0f;


		/// <summary>
		/// Compilation of vertex and fragment shader.
		/// </summary>
		void compileShaders();

		/// <summary>
		/// Setting the position of an object (Projection-View-Model).
		/// </summary>
		/// <param name="P">: projection matrix</param>
		/// <param name="V">: view matrix</param>
		/// <param name="M">: model matrix</param>
		void setPosition(const glm::mat4& P, const glm::mat4& V, const glm::mat4& M = glm::mat4(1));

		/// <summary>
		/// Getting the current world position according to mouse coordinates.
		/// </summary>
		/// <param name="x">: X coordinate</param>
		/// <param name="y">: Y coordinate</param>
		/// <returns>current world position</returns>
		Point currentWorldPosition(const uint x, const uint y);

		/// <summary>
		/// Calculation of a projection matrix.
		/// </summary>
		/// <param name="ratio">: aspect ratio included if true</param>
		/// <returns>projection matrix</returns>
		glm::mat4 calculateProjectionMatrix(const bool ratio);

		/// <summary>
		/// Calculation of a view matrix.
		/// </summary>
		/// <returns>view matrix</returns>
		glm::mat4 calculateViewMatrix() const;

		/// <summary>
		/// Displaying a potential OpenGL error.
		/// </summary>
		void displayError() const;

		/// <summary>
		/// Setting a color of rendered objects.
		/// </summary>
		/// <param name="color">: color in a vector of [red, green, blue, alpha]</param>
		void setColor(const LAS::Data::Vector4d& color);

		/// <summary>
		/// Point render procedure.
		/// </summary>
		/// <param name="P">: projection matrix</param>
		/// <param name="V">: view matrix</param>
		void paintPoints(const glm::mat4& P, const glm::mat4& V);

		/// <summary>
		/// Bounding box render procedure.
		/// </summary>
		/// <param name="P">: projection matrix</param>
		/// <param name="V">: view matrix</param>
		void paintBoundingBox(const glm::mat4& P, const glm::mat4& V);

		/// <summary>
		/// Current user selected plane render procedure.
		/// </summary>
		/// <param name="P">: projection matrix</param>
		/// <param name="V">: model matrix</param>
		void paintCurrentPlane(const glm::mat4& P, const glm::mat4& V);

	protected:
		/// <summary>
		/// OpenGL initialization procedure.
		/// </summary>
		void initializeGL() override;

		/// <summary>
		/// OpenGL render procedure.
		/// </summary>
		void paintGL() override;

		/// <summary>
		/// OpenGL resize procedure.
		/// </summary>
		/// <param name="w">: width of the widget</param>
		/// <param name="h">: height of the widget</param>
		void resizeGL(int w, int h) override;

		/// <summary>
		/// Mouse click event.
		/// </summary>
		/// <param name="event">captured event</param>
		void mousePressEvent(QMouseEvent* event) override;

		/// <summary>
		/// Mouse move event
		/// </summary>
		/// <param name="event">captured event</param>
		void mouseMoveEvent(QMouseEvent* event) override;


	public:
		/// <summary>
		/// Constructor of the widget.
		/// </summary>
		/// <param name="parent">: parent widget</param>
		PlanePicker(QWidget* parent);

		/// <summary>
		/// Destructor of the widget.
		/// </summary>
		~PlanePicker();

		/// <summary>
		/// Point vector setter.
		/// </summary>
		/// <param name="points">: points to be rendered</param>
		void setPoints(const std::vector<Point>& points);

		/// <summary>
		/// Voxel mesh setter.
		/// </summary>
		void setVoxelMesh(const VoxelMesh& voxelMesh);

		/// <summary>
		/// Plane getter.
		/// </summary>
		/// <returns>selected symmetry plane</returns>
		Plane getPlane() const;
	};
}