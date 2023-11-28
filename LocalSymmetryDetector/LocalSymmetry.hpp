#pragma once
#include <set>
#include <stack>
#include <string>
#include <vector>

#include "Structs/LineSegment.tpp"
#include "Structs/Point.tpp"
#include "Structs/Voxel.hpp"


namespace Symmetry {
    /// <summary>
    /// Base class for local symmetries.
    /// </summary>
    class LocalSymmetry {
    protected:
        /******************************************************************************************************************************************************/
        /* OBJECT VARIABLES                                                                                                                                   */
        /******************************************************************************************************************************************************/
        
        std::vector<Point> m_Points;             // Vector of LAS points.
        VoxelMesh m_VoxelMesh;                   // Object with the data about the voxel mesh.
        VoxelVector m_VoxelVector;               // 3D vector of voxels.
        std::vector<Voxel> m_InterestingVoxels;  // Vector of interesting voxels.
        std::vector<Voxel> m_MaterialVoxels;     // Vector of material voxels that are not interesting.



    public:
        /******************************************************************************************************************************************************/
        /* OBJECT METHODS                                                                                                                                     */
        /******************************************************************************************************************************************************/

        /// <summary>
        /// Points setter (creates voxel mesh in the background).
        /// </summary>
        /// <param name="points">: vector of points</param>
        /// <param name="moveToCoordinateCenter">: bounding box move to (0, 0, 0)</param>
        void setPoints(std::vector<Point>& points, bool moveToCoordinateCenter = true);

        /// <summary>
        /// Points getter.
        /// </summary>
        /// <returns>Vector of points</returns>
        const std::vector<Point>& getPoints() const;

        /// <summary>
        /// Voxel mesh getter.
        /// </summary>
        /// <returns>Voxel mesh</returns>
        VoxelMesh& getVoxelMesh();

        /// <summary>
        /// Voxel mesh getter.
        /// </summary>
        /// <returns>Voxel mesh</returns>
        VoxelMesh getVoxelMesh() const;

        /// <summary>
        /// Voxel vector getter.
        /// </summary>
        /// <returns>3D voxel vector</returns>
        VoxelVector& getVoxelVector();

        /// <summary>
        /// Voxel vector const getter.
        /// </summary>
        /// <returns>3D voxel vector</returns>
        VoxelVector getVoxelVector() const;

        /// <summary>
        /// Interesting voxels getter.
        /// </summary>
        /// <returns>Vector of interesting voxels</returns>
        std::vector<Voxel>& getInterestingVoxels();

        /// <summary>
        /// Material voxels getter.
        /// </summary>
        /// <returns>Vector of material voxels</returns>
        std::vector<Voxel>& getMaterialVoxels();



        /// <summary>
        /// Voxel mesh calculation according to voxel side length.
        /// </summary>
        /// <param name="userInputX">: desired length of X voxel side</param>
        /// <param name="userInputY">: desired length of Y voxel side</param>
        /// <param name="userInputZ">: desired length of Z voxel side</param>
        /// <returns>Voxel mesh</returns>
        VoxelMesh calculateVoxelMeshByVoxelSideLength(const uint userInputX, const uint userInputY, const uint userInputZ);

        /// <summary>
        /// Voxel mesh calculation according to maximal voxel count.
        /// </summary>
        /// <param name="voxels"></param>
        /// <param name="voxelMesh"></param>
        /// <param name="userInput"></param>
        void calculateVoxelMeshByMaximalVoxelCount(const uint userInput);

        /// <summary>
        /// Building a 3D voxel vector from a voxel mesh.
        /// </summary>
        /// <returns>3D voxel vector</returns>
        VoxelVector buildVoxelVector();

        /// <summary>
        /// Method for searching interesting voxels according to points in the 3D voxel vector.
        /// </summary>
        /// <param name="minClusterSize">: minimal size of interesting voxels cluster</param>
        /// <returns>(3D voxel vector, interesting voxels, material voxels)</returns>
        std::tuple<VoxelVector, std::vector<Voxel>, std::vector<Voxel>> findInterestingVoxels(const uint minClusterSize);

        /// <summary>
        /// Calculation of line segments between all pairs of points (fully connected graph).
        /// </summary>
        /// <param name="points">: vector of points</param>
        /// <param name="tolerance">: distance tolerance, within two floating point values are equal</param>
        /// <param name="minDistance">: minimum distance of a line segment</param>
        /// <param name="minMaterialRatio">: minimum value line segment material ratio</param>
        /// <param name="forbiddenAngles">: forbidden angle ranges where no line segments are formed</param>
        /// <returns>Vector of line segments</returns>
        std::vector<LineSegment> calculateLineSegmentsBetweenPoints(const std::vector<Point>& points, const double tolerance, const double minDistance, const double minMaterialRatio = 0.8, const std::vector<std::pair<double, double>>& forbiddenAngles = std::vector<std::pair<double, double>>());

        /// <summary>
        /// Splitting a line segment vector into parts, where the lenghts of line segments are below the tolerance.
        /// </summary>
        /// <param name="lineSegments">: vector of line segments</param>
        /// <param name="tolerance">: the desired tolerance</param>
        /// <returns>Vector of vectors of line segments with similar lengths</returns>
        SplitLineSegments splitLineSegmentVectorIntoParts(const std::vector<LineSegment>& lineSegments, const double tolerance);

        /// <summary>
        /// Extraction of points from material voxels.
        /// </summary>
        /// <returns>Vector of points</returns>
        std::vector<Point> pointsFromMaterialVoxels() const;

        /// <summary>
        /// Extraction of points from interesting voxels.
        /// </summary>
        /// <returns>Vector of points</returns>
        std::vector<Point> pointsFromInterestingVoxels() const;
    };


    /// <summary>
    /// Namespace for Symmetry related functions.
    /// </summary>
    namespace SymmetryFunctions {
        /// <summary>
        /// Getting a normalized 3D voxel vector (voxel edge size equals 1).
        /// </summary>
        /// <param name="voxelVector">: 3D voxel vector</param>
        /// <param name="voxelMesh">: voxel mesh</param>
        /// <returns>Normalized voxel vector</returns>
        VoxelVector getNormalizedVoxelVector(const VoxelVector& voxelVector, const VoxelMesh& voxelMesh);

        /// <summary>
        /// Clustering step.
        /// </summary>
        /// <param name="cluster">: current cluster</param>
        /// <param name="voxelVector">: 3D voxel vector</param>
        /// <param name="voxelMesh">: voxel mesh</param>
        /// <param name="stack">: stack of voxels in the vicinity</param>
        /// <param name="x">: current X position</param>
        /// <param name="y">: current Y position</param>
        /// <param name="z">: current Z position</param>
        void findMaterialClusterItem(std::set<Voxel>& cluster, VoxelVector& voxelVector, const VoxelMesh& voxelMesh, std::stack<Voxel>& stack, const uint x, const uint y, const uint z);

        /// <summary>
        /// Method for extraction of clusters of material voxels in a voxel mesh.
        /// </summary>
        /// <param name="normalizedVoxels">: normalized 3D voxel vector</param>
        /// <param name="voxelMesh">: voxel mesh</param>
        /// <returns>Vector of voxel clusters</returns>
        std::vector<std::set<Voxel>> findClustersOfMaterialVoxels(VoxelVector& normalizedVoxels, const VoxelMesh& voxelMesh);

        /// <summary>
        /// Method for removing small clusters of material voxels.
        /// </summary>
        /// <param name="materialVoxels">: vector of material voxels</param>
        /// <param name="voxels">: 3D voxel vector</param>
        /// <param name="voxelMesh">: voxel mesh</param>
        /// <param name="minimalClusterSize">: minimal cluster size</param>
        void removeSmallMaterialClusters(std::vector<Voxel>& materialVoxels, VoxelVector& voxels, const VoxelMesh& voxelMesh, const uint minimalClusterSize);


        /// <summary>
        /// Calculation of the ratio of material line segment sections.
        /// </summary>
        /// <param name="lineSegment">: input line segment</param>
        /// <param name="voxels">3D voxel vector</param>
        /// <param name="voxelMesh">: voxel mesh</param>
        /// <returns>Ratio of material sections</returns>
        double calculateMaterialSectionRatioInLineSegment(const LineSegment& lineSegment, const VoxelVector& voxels, const VoxelMesh& voxelMesh);
    }
};