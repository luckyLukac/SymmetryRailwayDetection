#pragma once
#include <set>
#include <tuple>
#include <vector>

#include "LocalSymmetry.hpp"
#include "ReflectionSymmetry.hpp"
#include "Structs/LineSegment.tpp"
#include "Structs/Plane.tpp"


namespace Symmetry {
    /// <summary>
    /// Class for storing and detecting reflection symmetries.
    /// </summary>
    class ReflectionSymmetry : public LocalSymmetry {
    private:
        /******************************************************************************************************************************************************/
        /* OBJECT VARIABLES                                                                                                                                   */
        /******************************************************************************************************************************************************/

        uint m_Index;
        Plane m_Plane;
        std::vector<LineSegment> m_LineSegments;
        std::vector<Voxel> m_SymmetryVoxels;
        std::vector<Voxel> m_AsymmetryVoxels;
        std::vector<Voxel> m_AsymmetryVoxelsAcross;



        /******************************************************************************************************************************************************/
        /* PRIVATE OBJECT METHODS                                                                                                                             */
        /******************************************************************************************************************************************************/

        /// <summary>
        /// Detection of a trivial reflection symmetry consisting of two line segments.
        /// </summary>
        /// <param name="l1">: first line segment</param>
        /// <param name="l2">: second line segment</param>
        /// <param name="tolerance">: the allowed distance tolerance</param>
        /// <param name="angleTolerance">: the allowed angle tolerance in degrees</param>
        /// <param name="forbiddenDistances">: forbidden distance ranges (between the two line segments)</param>
        /// <returns>(Reflection symmetry, validity)</returns>
        std::pair<ReflectionSymmetry, bool> detectTrivialReflectionSymmetry(const LineSegment& l1, const LineSegment& l2, const double tolerance, const double angleTolerance, const std::vector<std::pair<double, double>>& forbiddenDistances = std::vector<std::pair<double, double>>());                                                                                                                                                  // Checking whether two line segments represent a reflection symmetry.

        /// <summary>
        /// Detection of trivial symmetries based on the given line segments.
        /// </summary>
        /// <param name="lineSegments">: vector of line segments</param>
        /// <param name="tolerance">: distance tolerance</param>
        /// <param name="angleTolerance">: angle tolerance in degrees</param>
        /// <param name="forbiddenDistances">: forbidden distance ranges (between pairs of line segments)</param>
        /// <returns>Vector of trivial reflection symmetries</returns>
        std::vector<ReflectionSymmetry> detectTrivialSymmetries(const std::vector<LineSegment>& lineSegments, const double tolerance, const double angleTolerance, const std::vector<std::pair<double, double>>& forbiddenDistances = std::vector<std::pair<double, double>>());                                                                                                                                                                  // Finding simple symmetries.

        /// <summary>
        /// Merge of reflection symmetries that have an equal symmetry plane.
        /// </summary>
        /// <param name="reflectionSymmetries">: vector of reflection symmetries</param>
        /// <param name="tolerance">: distance tolerance</param>
        /// <param name="angleTolerance">: angle tolerance in degrees</param>
        /// <returns>Vector of merged reflection symmetries</returns>
        std::vector<ReflectionSymmetry> mergeTrivialSymmetries(std::vector<ReflectionSymmetry>& reflectionSymmetries, const double tolerance, const double angleTolerance);

        /// <summary>
        /// Adding voxels that form the reflection symmetry to the object.
        /// </summary>
        void addSymmetryVoxels();

        /// <summary>
        /// Adding material voxels to the object.
        /// <param name="maxDistanceFromPlane">: maximum distance from a plane to the voxel center</param>
        /// </summary>
        void addMaterialVoxels(const double maxDistanceFromPlane = MAX);

        /// <summary>
        /// Removal of too small clusters of symmetry voxels.
        /// </summary>
        /// <param name="minimumClusterSize">: minimum cluster size in voxels</param>
        void removeSmallClusters(const uint minimumClusterSize);

        /// <summary>
        /// Adding the interesting cluster voxel to the cluster.
        /// </summary>
        /// <param name="cluster">: current voxel cluster</param>
        /// <param name="voxelVector">: 3D voxel vector</param>
        /// <param name="voxelMesh">: voxel mesh</param>
        /// <param name="stack">: current stack of voxels</param>
        /// <param name="x">: X coordinate</param>
        /// <param name="y">: Y coordinate</param>
        /// <param name="z">: Z coordinate</param>
        void addInterestingClusterItem(std::set<Voxel>& cluster, VoxelVector& voxelVector, const VoxelMesh& voxelMesh, std::stack<Voxel>& stack, const uint x, const uint y, const uint z);                                                                                                         // Clustering step.

        /// <summary>
        /// Recursive cluster search step.
        /// </summary>
        /// <param name="cluster">: current voxel cluster</param>
        /// <param name="voxelVector">: 3D voxel vector</param>
        /// <param name="x">: X coordinate</param>
        /// <param name="y">: Y coordinate</param>
        /// <param name="z">: Z coordinate</param>
        void clusterRecursiveStep(std::set<Voxel>& cluster, VoxelVector& voxelVector, const uint x, const uint y, const uint z) const;

        /// <summary>
        /// Function for the detection of clusters of symmetry voxels.
        /// </summary>
        /// <returns>Vector of clusters</returns>
        std::vector<std::set<Voxel>> findClusters();

        /// <summary>
        /// Processing of possible points that lie on the voxel edge and the symmetry plane simultaneously.
        /// </summary>
        void processPointsOnThePlaneAndVoxelEdge();

        /// <summary>
        /// Postprocessing of the obtained symmetries (removing empty symmetries, indexing ...).
        /// </summary>
        /// <param name="symmetries">: vector of symmetries</param>
        void postprocessSymmetries(std::vector<ReflectionSymmetry>& symmetries);

    public:
        /******************************************************************************************************************************************************/
        /* CONSTRUCTORS                                                                                                                                       */
        /******************************************************************************************************************************************************/
        /// <summary>
        /// Default constructor.
        /// </summary>
        ReflectionSymmetry();

        /// <summary>
        /// Constructor with a symmetry plane and a vector of line segments.
        /// </summary>
        /// <param name="plane">: symmetry plane</param>
        /// <param name="lineSegments">: line segments</param>
        ReflectionSymmetry(const Plane& plane, const std::vector<LineSegment>& lineSegments);

        /// <summary>
        /// Constructor with a symmetry plane and a vector of voxels.
        /// </summary>
        /// <param name="plane">: symmetry plane</param>
        /// <param name="voxels">: symmetry voxels</param>
        ReflectionSymmetry(const Plane& plane, const std::vector<Voxel>& voxels);

        /// <summary>
        /// Constructor for merging two reflection symmetries into one.
        /// </summary>
        /// <param name="symmetry1"></param>
        /// <param name="symmetry2"></param>
        ReflectionSymmetry(const ReflectionSymmetry& symmetry1, const ReflectionSymmetry& symmetry2);

        /// <summary>
        /// Operator == is used for checking whether the two symmetries are equal.
        /// </summary>
        /// <param name="symmetry">: compared symmetry</param>
        /// <returns>True if both symmetries are equal, false otherwise</returns>
        bool operator == (const ReflectionSymmetry& symmetry) const;

        /// <summary>
        /// Operator &= is used for merging the current symmetry with another.
        /// </summary>
        /// <param name="symmetry">: symmetry to be merged</param>
        /// <returns>Merged reflection symmetry</returns>
        ReflectionSymmetry operator &= (const ReflectionSymmetry& symmetry) const;



        /******************************************************************************************************************************************************/
        /* GETTERS AND SETTERS                                                                                                                                */
        /******************************************************************************************************************************************************/

        /// <summary>
        /// Getter of the symmetry index.
        /// </summary>
        /// <returns>Symmetry index</returns>
        uint getIndex() const;

        /// <summary>
        /// Setter of the symmetry index.
        /// </summary>
        /// <param name="index">: symmetry index</param>
        void setIndex(const uint index);

        /// <summary>
        /// Getter of the symmetry plane.
        /// </summary>
        /// <returns>Symmetry plane</returns>
        Plane getPlane() const;

        /// <summary>
        /// Setter of the symmetry plane.
        /// </summary>
        /// <param name="plane"></param>
        void setPlane(const Plane& plane);

        /// <summary>
        /// Getter of symmetry line segments.
        /// </summary>
        /// <returns>Line segments</returns>
        std::vector<LineSegment> getLineSegments() const;

        /// <summary>
        /// Getter of line segment count.
        /// </summary>
        /// <returns>Number of line segments</returns>
        u128 getLineSegmentCount() const;

        /// <summary>
        /// Adder of a new line segment to the symmetry.
        /// </summary>
        /// <param name="lineSegment">: line segment to be added</param>
        void addLineSegment(const LineSegment& lineSegment);

        /// <summary>
        /// Getter of symmetry voxels.
        /// </summary>
        /// <returns>Vector of symmetry voxels</returns>
        std::vector<Voxel> getSymmetryVoxels() const;

        /// <summary>
        /// Getter of obstacle voxels.
        /// </summary>
        /// <returns>Vector of obstacle voxels</returns>
        std::vector<Voxel> getObstacleVoxels() const;

        /// <summary>
        /// Getter of voxels across the obstacle.
        /// </summary>
        /// <returns>Vector of voxels across the obstacle</returns>
        std::vector<Voxel> getVoxelsAcrossObstacle() const;

        /// <summary>
        /// Setter of symmetry voxels.
        /// </summary>
        /// <param name="voxels">: new symmetry voxel vector</param>
        void setSymmetryVoxels(const std::vector<Voxel>& voxels, const std::vector<Voxel>& aVoxels, const std::vector<Voxel>& aVoxelsAcross);

        /// <summary>
        /// Getter of symmetry voxel count.
        /// </summary>
        /// <returns>Number of symmetry voxels</returns>
        u128 getVoxelCount() const;

        /// <summary>
        /// Setter of voxel vector.
        /// </summary>
        /// <param name="voxelVector">: new voxel vector</param>
        void setVoxelVector(const VoxelVector& voxelVector);

        /// <summary>
        /// Adder of a new symmetry voxel to the symmetry.
        /// </summary>
        /// <param name="voxel">: voxel to be added</param>
        void addSymmetryVoxel(const Voxel& voxel);

        /// <summary>
        /// Check whether a voxel contains at least one point.
        /// </summary>
        /// <param name="voxel">: arbitrary voxel</param>
        /// <param name="vm">: voxel mesh</param>
        /// <param name="points">: list of points</param>
        /// <returns>True if voxel is material, false otherwise</returns>
        bool doesVoxelContainPoint(const Voxel& voxel, const VoxelMesh& vm, const std::vector<Point>& points);

        /// <summary>
        /// Retrieving number of points in a voxel.
        /// </summary>
        /// <param name="voxel">: arbitrary voxel</param>
        /// <param name="vm">: voxel mesh</param>
        /// <param name="points">: list of points</param>
        /// <returns>True if voxel is material, false otherwise</returns>
        uint numberOfPointsInVoxel(const Voxel& voxel, const VoxelMesh& vm, const std::vector<Point>& points);

        /******************************************************************************************************************************************************/
        /* PUBLIC OBJECT METHODS                                                                                                                              */
        /******************************************************************************************************************************************************/

        /// <summary>
        /// Detection of partial reflection symmetries in a point cloud.
        /// </summary>
        /// <param name="voxelSideX">: desired X side of a voxel</param>
        /// <param name="voxelSideY">: desired Y side of a voxel</param>
        /// <param name="voxelSideZ">: desired Z side of a voxel</param>
        /// <param name="tolerance">: distance tolerance</param>
        /// <param name="angleTolerance">: angle tolerance in degrees</param>
        /// <param name="minimumClusterSize">: minimum size of a voxel cluster</param>
        /// <param name="minimumLineSegmentDistance">: minimum length of a line segment</param>
        /// <param name="postprocessing">: postprocessing of symmetries if true</param>
        /// <returns>Vector of reflection symmetries</returns>
        std::vector<ReflectionSymmetry> detectReflectionSymmetries(const uint voxelSideX, const uint voxelSideY, const uint voxelSideZ, const double tolerance, const double angleTolerance, const uint minimumClusterSize, const double minimumLineSegmentDistance, const bool postprocessing = true);

        /// <summary>
        /// Detection of local symmetries from a partial reflection symmetries.
        /// </summary>
        /// <param name="reflectionSymmetries">: vector of detected reflection symmetries</param>
        /// <param name="minSymmetrySize">: minimum number of voxels in a local symmetry</param>
        /// <returns>Vector of local reflection symmetries</returns>
        std::vector<ReflectionSymmetry> detectLocalSymmetries(std::vector<ReflectionSymmetry>& reflectionSymmetries, const uint minSymmetrySize);

        /// <summary>
        /// Calculation of point positions according to the symmetry plane.
        /// </summary>
        /// <param name="points">: vector of points</param>
        /// <param name="voxelVector">: 3D voxel vector</param>
        /// <param name="voxelMesh">: voxel mesh</param>
        /// <param name="plane">: symmetry plane</param>
        /// <param name="validPlane">: validity of symmetry plane</param>
        /// <param name="normalize">: normalization if true</param>
        /// <returns>Point position vector</returns>
        std::vector<PointPosition> calculatePositionsFromPlane(const std::vector<Point>& points, const VoxelVector& voxelVector, const VoxelMesh& voxelMesh, const Plane& plane, const bool validPlane, const bool normalize) const;

        /// <summary>
        /// Calculation of nearest symmetries according to the given plane.
        /// </summary>
        /// <param name="symmetries">: vector of symmetries</param>
        /// <param name="desiredPlane">: the given plane</param>
        /// <param name="voxelMesh">: voxel mesh</param>
        /// <returns>Sorted symmetries</returns>
        std::vector<ReflectionSymmetry> calculateNearestSymmetriesByPlane(const std::vector<ReflectionSymmetry>& symmetries, const Plane& desiredPlane, const VoxelMesh& voxelMesh);
    };
};