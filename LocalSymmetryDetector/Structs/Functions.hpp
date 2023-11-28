#pragma once

#include <set>

#include "Point.tpp"
#include "Voxel.hpp"


namespace Symmetry {
    /// <summary>
    /// Namespace for Point related functions.
    /// </summary>
    namespace PointFunctions {
        /// <summary>
        /// Calculation of the Euclidean distance between two points.
        /// </summary>
        /// <param name="p1">: first point</param>
        /// <param name="p2">: second point</param>
        /// <returns>Euclidean distance</returns>
        double distance(const Point& p1, const Point& p2);

        /// <summary>
        /// Calculation of the midpoint between two points.
        /// </summary>
        /// <param name="p1">: first point</param>
        /// <param name="p2">: second point</param>
        /// <returns>Middle point</returns>
        Point calculateMidPoint(const Point& p1, const Point& p2);

        /// <summary>
        /// Calculation of a barycenter (centroid) point of a polygon, defined with the set of points.
        /// </summary>
        /// <param name="points"></param>
        /// <returns>Barycenter point</returns>
        Point calculateBarycenterPoint(const std::set<Point>& points);

        /// <summary>
        /// Calculation of the voxel position in the voxel mesh according to the given point.
        /// </summary>
        /// <param name="point">: the given point</param>
        /// <param name="voxelMesh">: voxel mesh</param>
        /// <returns>Voxel object with indices</returns>
        Voxel voxelPositionFromPoint(const Point& point, const VoxelMesh& voxelMesh);
    };


    /// <summary>
    /// Namespace for Vector related functions.
    /// </summary>
    namespace VectorFunctions {
        /// <summary>
        /// Cross product between two vectors.
        /// </summary>
        /// <param name="v1">: first vector</param>
        /// <param name="v2">: second vector</param>
        /// <returns>Vector, produced with cross product</returns>
        Vector3d crossProduct(const Vector3d& v1, const Vector3d& v2);

        /// <summary>
        /// Calculation of the angle between two vectors.
        /// </summary>
        /// <param name="v1">: first vector</param>
        /// <param name="v2">: second vector</param>
        /// <returns>Angle in radians</returns>
        double angle(const Vector3d& v1, const Vector3d& v2);
    };


    /// <summary>
    /// Namespace for Voxel related functions.
    /// </summary>
    namespace VoxelFunctions {
        /// <summary>
        /// Calculation of whether the voxel should be painted brightly according to its position.
        /// </summary>
        /// <param name="x">: X position in voxel mesh</param>
        /// <param name="y">: Y position in voxel mesh</param>
        /// <param name="z">: Z position in voxel mesh</param>
        /// <returns>True if the voxel should be bright, false otherwise</returns>
        bool isVoxelBright(const uint x, const uint y, const uint z);

        /// <summary>
        /// Calculation of the layer of the voxel mesh according to the Z coordinate.
        /// </summary>
        /// <param name="z">: Z coordinate</param>
        /// <param name="vm">: voxel mesh</param>
        /// <returns>Layer in voxel mesh</returns>
        uint calculateLayerInVoxelMesh(const double z, const VoxelMesh& vm);

        /// <summary>
        /// Function for clearing all set variables in a voxel vector.
        /// </summary>
        /// <param name="voxelVector">: input voxel vector</param>
        void clearVoxelVector(VoxelVector& voxelVector);

        /// <summary>
        /// Check for whether the voxel is in bounds of the voxel mesh.
        /// </summary>
        /// <param name="voxel">: input voxel</param>
        /// <param name="vm">: voxel mesh</param>
        /// <returns>True if voxel is inside of the voxel mesh, false otherwise</returns>
        bool isVoxelInsideVoxelMesh(const Voxel& voxel, const VoxelMesh& vm);
    }
}