#pragma once
#include <vector>

#include <LocalSymmetryDetector/ReflectionSymmetry.hpp>


using namespace Symmetry;


/// <summary>
/// Class for detection of rails and obstacles on a railway.
/// </summary>
class RailObstacleDetector {
private:
    std::vector<ReflectionSymmetry> m_Symmetries;
    std::vector<Point> m_Points;
    std::vector<Point> m_FilteredPoints;



    /// <summary>
    /// Merging of all symmetries.
    /// </summary>
    void mergeAllSymmetries();

    /// <summary>
    /// Basic rail detection using all symmetry line segments.
    /// </summary>
    /// <param name="symmetry">: detected local reflection symmetry</param>
    /// <returns>Vector of voxels that lie on rails</returns>
    std::vector<Voxel> basicRailDetection(ReflectionSymmetry& symmetry);

    /// <summary>
    /// Obtaining points upon which accurate plane is detected.
    /// </summary>
    /// <param name="railVoxels">: list of voxels that represent rails</param>
    /// <param name="voxelMesh">: voxel mesh</param>
    /// <returns>Vector of points for MSE</returns>
    std::vector<Point> obtainPointsForAccuratePlaneDetection(std::vector<Voxel> railVoxels, const VoxelMesh& voxelMesh) const;

    /// <summary>
    /// Calculation of the actual railway symmetry plane according to rail voxels using Minimal Square Error method.
    /// </summary>
    /// <param name="railVoxels">: list of voxels that represent rails</param>
    /// <param name="voxelMesh">: voxel mesh</param>
    /// <returns>Accurate symmetry plane</returns>
    Plane accuratePlane(const std::vector<Voxel>& railVoxels, const VoxelMesh& voxelMesh) const;

    /// <summary>
    /// Obstacle detection on the railway track.
    /// </summary>
    /// <param name="symmetry">: symmetry that represents a railway track</param>
    /// <returns>(rail symmetry voxels, voxels under the obstacle, voxels across the obstacle on another rail)</returns>
    std::tuple<std::vector<Voxel>, std::vector<Voxel>, std::vector<Voxel>> detectObstaclesOnRailway(ReflectionSymmetry& symmetry);

    /// <summary>
    /// Converting non-material voxels, which are material in the phase before the geometric point extraction, to material.
    /// </summary>
    /// <param name="symmetry">Railway scene</param>
    void addOriginalSceneMaterialVoxels(ReflectionSymmetry& symmetry);

public:
    /// <summary>
    /// Main constructor.
    /// </summary>
    /// <param name="symmetry"></param>
    RailObstacleDetector(const std::vector<ReflectionSymmetry>& symmetries, const std::vector<Point>& points, const std::vector<Point>& filteredPoints);

    /// <summary>
    /// Rail detection according to the given reflection symmetries.
    /// </summary>
    /// <param name="mergeAll">: merging all symmetries if true</param>
    /// <returns>Vector of symmetries with detected rails</returns>
    std::vector<ReflectionSymmetry> detectRails(const bool mergeAll = false);
};