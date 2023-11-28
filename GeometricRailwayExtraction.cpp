#include <LocalSymmetryDetector/LocalSymmetry.hpp>

#include "GeometricRailwayExtraction.hpp"



VoxelMesh GeometricRailwayExtraction::calculate2DVoxelMeshByVoxelSideLength(const uint userInputX, const uint userInputY, const double deltaX, const double deltaY) {
    // Calculating X, Y, Z and total voxel count.
    const uint x = static_cast<uint>(floor(deltaX / userInputX)) + 1;  // X coordinate voxel count.
    const uint y = static_cast<uint>(floor(deltaY / userInputY)) + 1;  // Y coordinate voxel count.
    const uint count = x * y;                                          // Total voxel count.

    VoxelMesh voxelMesh;
    voxelMesh.count = count;
    voxelMesh.sideX = userInputX;
    voxelMesh.sideY = userInputY;
    voxelMesh.sideZ = 1;
    voxelMesh.xSize = x;
    voxelMesh.ySize = y;
    voxelMesh.zSize = 1;

    return voxelMesh;
}



GeometricRailwayExtraction::GeometricRailwayExtraction(const std::vector<Point>& points, const VoxelMesh& voxelMesh) : m_Points(points) {
    m_VoxelMesh = calculate2DVoxelMeshByVoxelSideLength(20, 20, voxelMesh.deltaX, voxelMesh.deltaY);
    m_VoxelGrid = std::vector<std::vector<VoxelWithPoints>>(m_VoxelMesh.ySize, std::vector<VoxelWithPoints>(m_VoxelMesh.xSize));
}

std::vector<Point> GeometricRailwayExtraction::extractCandidateRailwayPoints(const double minHeightAboveGround, const double maxHeightAboveGround) {
    // Iterating through all input point of the loaded point cloud.
    for (const Point& point : m_Points) {
        // Getting voxel grid X and Y coordinates.
        const uint voxelX = static_cast<uint>(point.x / m_VoxelMesh.sideX);
        const uint voxelY = static_cast<uint>(point.y / m_VoxelMesh.sideY);

        // Detection of a new minimum Z at the current grid position.
        if (Tolerance::isInTolerance(m_VoxelGrid[voxelY][voxelX].minZ, -1, 0.00001) || point.z < m_VoxelGrid[voxelY][voxelX].minZ) {
            m_VoxelGrid[voxelY][voxelX].minZ = point.z;
        }
        // Detection of a new maximum Z at the current grid position.
        if (Tolerance::isInTolerance(m_VoxelGrid[voxelY][voxelX].minZ, -1, 0.00001) || point.z > m_VoxelGrid[voxelY][voxelX].maxZ) {
            m_VoxelGrid[voxelY][voxelX].maxZ = point.z;
        }
    }

    // Second iteration through all points and extraction of candidates.
    std::vector<Point> candidatePoints;
    for (const Point& point : m_Points) {
        // Getting voxel grid X and Y coordinates.
        const uint voxelX = static_cast<uint>(point.x / m_VoxelMesh.sideX);
        const uint voxelY = static_cast<uint>(point.y / m_VoxelMesh.sideY);

        // If the point is inside the allowed tolerance, it becomes a candidate for a railway. Yipeeeeeee!
        if (Difference::difference(point.z, m_VoxelGrid[voxelY][voxelX].minZ) > minHeightAboveGround &&
            Difference::difference(point.z, m_VoxelGrid[voxelY][voxelX].minZ) < maxHeightAboveGround
        ) 
        {
            candidatePoints.push_back(point);
        }
    }

    return candidatePoints;
}
