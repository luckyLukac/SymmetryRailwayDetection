#include <algorithm>

#include "Functions.hpp"

using namespace Symmetry;


/******************************************************************************************************************************************************/
/* POINT FUNCTIONS                                                                                                                                    */
/******************************************************************************************************************************************************/

double PointFunctions::distance(const Point& p1, const Point& p2) {
    return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2) + std::pow(p1.z - p2.z, 2));
}

Point PointFunctions::calculateMidPoint(const Point& p1, const Point& p2) {
    // Midpoint can be calculated as the "average point" between two points: (p1 + p2) / 2.
    Point midpoint((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0, (p1.z + p2.z) / 2.0);

    return midpoint;
}

Point PointFunctions::calculateBarycenterPoint(const std::set<Point>& points) {
    Point point(0, 0, 0);

    // Barycenter is calculated as sum of all points' coordinates divided by a number of points.
    std::for_each(points.begin(), points.end(), [&point](const Point& currentPoint) { point += currentPoint; });
    point /= static_cast<double>(points.size());

    return point;
}

Voxel PointFunctions::voxelPositionFromPoint(const Point& point, const VoxelMesh& voxelMesh) {
    const uint x = static_cast<uint>(point.x / voxelMesh.sideX) * voxelMesh.sideX;
    const uint y = static_cast<uint>(point.y / voxelMesh.sideY) * voxelMesh.sideY;
    const uint z = static_cast<uint>(point.z / voxelMesh.sideZ) * voxelMesh.sideZ;

    return Voxel(x, y, z);
}



/******************************************************************************************************************************************************/
/* VECTOR FUNCTIONS                                                                                                                                   */
/******************************************************************************************************************************************************/

Vector3d VectorFunctions::crossProduct(const Vector3d& v1, const Vector3d& v2) {
    return {
        v1.y * v2.z - v1.z * v2.y,
        v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x
    };
}

double VectorFunctions::angle(const Vector3d& v1, const Vector3d& v2) {
    const Vector3d vn1 = v1.normalize();
    const Vector3d vn2 = v2.normalize();

    const double dot = vn1.x * vn2.x + vn1.y * vn2.y;
    const double det = vn1.x * vn2.y - vn1.y * vn2.x;
    const double angle = std::atan2(det, dot);

    // Angle lies inside the interval [0, 2 * PI].
    return angle >= 0 ? angle : angle + 2 * PI;
}



/******************************************************************************************************************************************************/
/* VOXEL FUNCTIONS                                                                                                                                    */
/******************************************************************************************************************************************************/

bool VoxelFunctions::isVoxelBright(const uint x, const uint y, const uint z) {
    // First bright option.
    if (z % 2 == 0 && y % 2 == 0 && x % 2 == 0) {
        return true;
    }

    // Second bright option.
    if (z % 2 == 0 && y % 2 == 1 && x % 2 == 1) {
        return true;
    }

    // Third bright option.
    if (z % 2 == 1 && y % 2 == 1 && x % 2 == 0) {
        return true;
    }

    // Fourth bright option.
    if (z % 2 == 1 && y % 2 == 0 && x % 2 == 1) {
        return true;
    }

    return false;
}

uint VoxelFunctions::calculateLayerInVoxelMesh(const double z, const VoxelMesh& vm) {
    const uint layer = static_cast<uint>((z - vm.minZ) / vm.sideZ);
    return layer;
}

void VoxelFunctions::clearVoxelVector(VoxelVector& voxelVector) {
    for (auto& vector2D : voxelVector) {
        for (auto& vector1D : vector2D) {
            for (Voxel& voxel : vector1D) {
                voxel.resetAllProperties();
            }
        }
    }
}

bool VoxelFunctions::isVoxelInsideVoxelMesh(const Voxel& voxel, const VoxelMesh& vm) {
    if (voxel.x < 0 || voxel.x > vm.maxX || voxel.y < 0 || voxel.y > vm.maxY || voxel.z < 0 || voxel.z > vm.maxZ) {
        return false;
    }

    return true;
}