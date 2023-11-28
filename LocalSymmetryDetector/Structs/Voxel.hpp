#pragma once
#include <array>
#include <limits>
#include <vector>

#include "Constants.hpp"
#include "Point.tpp"


namespace Symmetry {
    /// <summary>
    /// Structure for voxel mesh description.
    /// </summary>
    struct VoxelMesh {
        uint count = 0;                                    // Number of voxels.
        uint sideX = 0;                                    // Length of a voxel edge by X.
        uint sideY = 0;                                    // Length of a voxel edge by Y.
        uint sideZ = 0;                                    // Length of a voxel edge by Z.
        uint xSize = 0;                                    // Number of voxels by X.
        uint ySize = 0;                                    // Number of voxels by Y.
        uint zSize = 0;                                    // Number of voxels by Z.
        double minX = std::numeric_limits<double>::max();  // Minimum X.
        double maxX = std::numeric_limits<double>::min();  // Maximum X.
        double minY = std::numeric_limits<double>::max();  // Minimum Y.
        double maxY = std::numeric_limits<double>::min();  // Maximum Y.
        double minZ = std::numeric_limits<double>::max();  // Minimum Z.
        double maxZ = std::numeric_limits<double>::min();  // Maximum Z.
        double deltaX = 0.0;                               // Distance between the maximum and the minimum X coordinate.
        double deltaY = 0.0;                               // Distance between the maximum and the minimum Y coordinate.
        double deltaZ = 0.0;                               // Distance between the maximum and the minimum Z coordinate.
    };


    /// <summary>
    /// Voxel structure.
    /// </summary>
    struct Voxel {
        /******************************************************************************************************************************************************/
        /* STRUCT VARIABLES                                                                                                                                   */
        /******************************************************************************************************************************************************/

        uint x;                              // X coordinate.
        uint y;                              // Y coordinate.
        uint z;                              // Z coordinate.
        bool interesting;                    // True if the voxel contains at least one point and is not surrounded with other interesting voxels.
        bool material;                       // True if the voxel contains at least one point.
        bool checked;                        // True if the voxel has been checked during the geometry search.
        bool inSymmetry;                     // True if the voxel is a part of a local symmetry.
        bool inAsymmetry = false;                    // True if the voxel is not a part of a local symmetry.
        bool inAsymmetryAcross = false;                    // True if the voxel is not a part of a local symmetry.



        /******************************************************************************************************************************************************/
        /* CONSTRUCTORS                                                                                                                                       */
        /******************************************************************************************************************************************************/

        /// <summary>
        /// Default constructor.
        /// </summary>
        Voxel() :
            x(0),
            y(0),
            z(0),
            interesting(false),
            material(false),
            checked(false),
            inSymmetry(false) {}

        /// <summary>
        /// Constructor with X, Y and Z coordinates.
        /// </summary>
        /// <param name="x">: X coordinate</param>
        /// <param name="y">: Y coordinate</param>
        /// <param name="z">: Z coordinate</param>
        Voxel(uint x, uint y, uint z) :
            x(x),
            y(y),
            z(z),
            interesting(false),
            material(false),
            checked(false),
            inSymmetry(false) {}

        /// <summary>
        /// Constructor with X, Y and Z coordinates along with interesting and material properties.
        /// </summary>
        /// <param name="x">: X coordinate</param>
        /// <param name="y">: Y coordinate</param>
        /// <param name="z">: Z coordinate</param>
        /// <param name="interesting">: interesting property</param>
        /// <param name="material">: material property</param>
        Voxel(const uint x, const uint y, const uint z, const bool interesting, const bool material) :
            x(x),
            y(y),
            z(z),
            interesting(interesting),
            material(material),
            checked(false),
            inSymmetry(false) {}

        /// <summary>
        /// Constructor with X, Y and Z coordinates along with symmetry property.
        /// </summary>
        /// <param name="x">: X coordinate</param>
        /// <param name="y">: Y coordinate</param>
        /// <param name="z">: Z coordinate</param>
        /// <param name="inSymmetry">: part of a symmetry if true</param>
        Voxel(const uint x, const uint y, const uint z, const bool inSymmetry) :
            x(x),
            y(y),
            z(z),
            interesting(false),
            material(false),
            checked(false),
            inSymmetry(inSymmetry) {}

        /// <summary>
        /// Copy constructor from IntPoint and a boolean for an interesting voxel.
        /// </summary>
        /// <param name="intPoint">: integer point</param>
        /// <param name="interesting">: interesting property</param>
        Voxel(const IntPoint& intPoint, const bool interesting) :
            x(intPoint.x),
            y(intPoint.y),
            z(intPoint.z),
            interesting(interesting),
            material(interesting),
            checked(false),
            inSymmetry(false) {}

        /// <summary>
        /// Copy constructor from Point.
        /// </summary>
        /// <param name="point">: floating value point</param>
        Voxel(const Point& point) :
            x(static_cast<uint>(point.x)),
            y(static_cast<uint>(point.y)),
            z(static_cast<uint>(point.z)),
            interesting(false),
            material(false),
            checked(false),
            inSymmetry(false) {}



        /******************************************************************************************************************************************************/
        /* STRUCT METHODS                                                                                                                                     */
        /******************************************************************************************************************************************************/

        /// <summary>
        /// Operator for voxel equality check (current voxel with the given one).
        /// </summary>
        /// <param name="voxel">: the given voxel</param>
        /// <returns>True if voxels are equal, false otherwise</returns>
        bool operator == (const Voxel& voxel) const {
            if (x == voxel.x && y == voxel.y && z == voxel.z) {
                return true;
            }

            return false;
        }

        /// <summary>
        /// Operator for the sort method as the comparator.
        /// </summary>
        /// <param name="voxel">: second voxel</param>
        /// <returns>True if current voxel has smaller coordinates than the given one, false otherwise</returns>
        bool operator < (const Voxel& voxel) const {
            return (
                (x < voxel.x) ||
                (x == voxel.x && y < voxel.y) ||
                (x == voxel.x && y == voxel.y && z < voxel.z)
                );
        }

        /// <summary>
        /// Multiplication operator with the given scalar.
        /// </summary>
        /// <param name="factor">: the given scalar</param>
        /// <returns>Multiplied voxel</returns>
        Voxel operator * (const uint factor) const {
            return Voxel(x * factor, y * factor, z * factor);
        }

        /// <summary>
        /// Calculation of a center coordinate of the voxel.
        /// </summary>
        /// <param name="voxelMesh">: voxel mesh</param>
        /// <returns>Center point</returns>
        Point centerPoint(const VoxelMesh& voxelMesh) const {
            const double pointX = x + (0.5 * voxelMesh.sideX);
            const double pointY = y + (0.5 * voxelMesh.sideY);
            const double pointZ = z + (0.5 * voxelMesh.sideZ);

            return Point(pointX, pointY, pointZ);
        }

        /// <summary>
        /// Method for retrieving normalized coordinates of the voxel.
        /// </summary>
        /// <param name="voxelMesh">: voxel mesh</param>
        /// <returns>Array of coordinates</returns>
        std::array<uint, 3> getNormalizedCoordinates(const VoxelMesh& voxelMesh) const {
            const uint xCoordinate = static_cast<uint>(x / voxelMesh.sideX);
            const uint yCoordinate = static_cast<uint>(y / voxelMesh.sideY);
            const uint zCoordinate = static_cast<uint>(z / voxelMesh.sideZ);

            return std::array<uint, 3>({ xCoordinate, yCoordinate, zCoordinate });
        }

        /// <summary>
        /// Method for retrieving coordinates of the voxel.
        /// </summary>
        /// <returns>Array of coordinates</returns>
        std::array<uint, 3> getCoordinates() const {
            return std::array<uint, 3>({ x, y, z });
        }

        /// <summary>
        /// Method for dividing voxel coordinates with a scalar.
        /// </summary>
        /// <param name="factor">: scalar value</param>
        /// <returns>New voxel</returns>
        Voxel divideVoxelWithScalar(const uint factor) const {
            const uint xCoordinate = static_cast<uint>(x / factor);
            const uint yCoordinate = static_cast<uint>(y / factor);
            const uint zCoordinate = static_cast<uint>(z / factor);

            return Voxel(xCoordinate, yCoordinate, zCoordinate);
        }

        /// <summary>
        /// Method for dividing voxel coordinates with a vector.
        /// </summary>
        /// <param name="v">: 3D vector</param>
        /// <returns>New voxel</returns>
        Voxel divideVoxelWithVector(const Vector3d& v) const {
            const uint xCoordinate = static_cast<uint>(x / v.x);
            const uint yCoordinate = static_cast<uint>(y / v.y);
            const uint zCoordinate = static_cast<uint>(z / v.z);

            return Voxel(xCoordinate, yCoordinate, zCoordinate);
        }

        /// <summary>
        /// Reset of all voxel properties to their default values.
        /// </summary>
        void resetAllProperties() {
            material = interesting = inSymmetry = checked = false;
        }
    };


    struct VoxelWithPoints : public Voxel {
        double minZ;
        double maxZ;

        VoxelWithPoints() : minZ(-1), maxZ(-1) {};
    };



    /******************************************************************************************************************************************************/
    /* ALIASES                                                                                                                                            */
    /******************************************************************************************************************************************************/

    using VoxelVector = std::vector<std::vector<std::vector<Voxel>>>;
    using VoxelGrid = std::vector<std::vector<VoxelWithPoints>>;
};