#pragma once

#include <LocalSymmetryDetector/Structs/Point.tpp>
#include <LocalSymmetryDetector/Structs/Voxel.hpp>

using namespace Symmetry;


/// <summary>
/// Class for geometric extraction of candidate railway points from the point cloud.
/// </summary>
class GeometricRailwayExtraction {
private:
    std::vector<Point> m_Points;
    VoxelMesh m_VoxelMesh;
    VoxelGrid m_VoxelGrid;


    /// <summary>
    /// 2D voxel mesh calculation according to voxel side lengths.
    /// </summary>
    /// <param name="userInputX">: desired length of X voxel side</param>
    /// <param name="userInputY">: desired length of Y voxel side</param>
    /// <param name="deltaX">: distance between minimum and maximum points by X</param>
    /// <param name="deltaY">: distance between minimum and maximum points by Y</param>
    VoxelMesh calculate2DVoxelMeshByVoxelSideLength(const uint userInputX, const uint userInputY, const double deltaX, const double deltaY);

public:
    /// <summary>
    /// Main constructor of the class.
    /// </summary>
    /// <param name="points">: point cloud with railways</param>
    /// <param name="voxelMesh">: scene voxel mesh</param>
    GeometricRailwayExtraction(const std::vector<Point>& points, const VoxelMesh& voxelMesh);

    /// <summary>
    /// Extraction of points that represent candidates for railways.
    /// </summary>
    /// <param name="minHeightAboveGround">: minimum height above the ground</param>
    /// <param name="maxHeightAboveGround">: maximum height above the ground</param>
    /// <returns>Candidate railway points</returns>
    std::vector<Point> extractCandidateRailwayPoints(const double minHeightAboveGround = 10, const double maxHeightAboveGround = 20);
};