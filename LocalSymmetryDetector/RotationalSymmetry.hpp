#ifndef ROTATIONALSYMMETRY_H
#define ROTATIONALSYMMETRY_H

#include <vector>

#include "Structs/LineSegment.tpp"
#include "Structs/Plane.tpp"
#include "Structs/Voxel.hpp"
#include "LocalSymmetry.hpp"


namespace Symmetry {
    // CLASS DECLARATIONS
    class Circumcircle;
    class Edge;
    class LineSegmentGraph;
    class LineSegmentPair;
    class LineSegmentPairGroup;
    class LineSegmentSubgroup;
    class RotationalSymmetry;
    class Vertex;

    // ALIASES
    using SplitLineSegmentPairs = std::vector<std::vector<std::vector<std::vector<LineSegmentPair>>>>;

    // CLASSES
    // Class for rotational symmetries.
    class RotationalSymmetry : public LocalSymmetry {
    private:
        // OBJECT VARIABLES
        unsigned int rotation;                  // Rotation level (1, 2, ..., n).
        Point symmetryAxis;                     // Symmetry axis point.
        std::vector<LineSegment> lineSegments;  // Line segments of the rotational symmetry.
        std::vector<Voxel> voxels;              // Vector of voxels, included in the symmetry.

        // CLASS METHODS
        SplitLineSegmentPairs splitLineSegmentsInPairs(const std::vector<LineSegment>& lineSegments, const VoxelMesh& voxelMesh, const double lengthTolerance, const double angleTolerance);  // Finding simple symmetries.
        void buildCombinationsVectorRecursively(std::vector<std::vector<int> >& combinations, std::vector<int>& tempVector, int N, int left, int K);  // Getting one recursive vector combination.
        std::vector<std::vector<int>> allVectorCombinations(const int N, const int K);  // Getting all vector combinations with size K from a vector with the size of N.
        std::vector<LineSegmentPairGroup> splitPairsIntoGroups(const SplitLineSegmentPairs& pairs);  // Splitting line segment pairs by length between pair center points.
        std::vector<LineSegmentSubgroup> splitLineSegmentGroupsIntoSubgroups(const std::vector<LineSegmentPairGroup>& groups);  // Splitting line segment groups into subgroups.
        std::vector<Circumcircle> buildLineSegmentMidpointsCircumcircles(const SplitLineSegmentPairs& pairs);  // Building line segment midpoints circumcircles.
        std::vector<RotationalSymmetry> findSimpleSymmetries(const std::vector<Circumcircle>& circumcircles, const double angleTolerance);  // Finding simple rotational symmetries from circumcircles.
        std::vector<RotationalSymmetry> mergeSymmetries(std::vector<RotationalSymmetry>& symmetries);  // Merging symmetries by layers, where symmetries with the same rotation axis and rotation level are merged.
        void getVoxelsInSymmetries(std::vector<RotationalSymmetry>& rotationalSymmetries, const VoxelVector& voxelVector, const VoxelMesh& voxelMesh);  // Getting voxels for each rotational symmetry.
    public:
        // OBJECT METHODS
        RotationalSymmetry(const Point& symmetryAxis, const unsigned int rotation, const std::vector<LineSegment>& lineSegments);  // Main constructor.
        std::vector<PointPosition> calculatePositionsFromAxis(const std::vector<Point>& points, const VoxelMesh& voxelMesh, bool validAxis);  // Calculating the point positions according to the axis.

        // GETTERS
        Point getSymmetryAxis() const;                            // Symmetry axis getter.
        unsigned int getRotation() const;                         // Rotation getter.
        const std::vector<LineSegment>& getLineSegments() const;  // Line segments getter.
        const std::vector<Voxel>& getVoxels() const;              // Voxel getter.

        // OVERLOADED OPERATORS
        RotationalSymmetry operator &= (RotationalSymmetry& s);   // Operator &= is used for merging the current symmetry with another.

        // CLASS METHODS
        std::vector<RotationalSymmetry> calculateRotationalSymmetries(const std::vector<Point>& points, const VoxelVector& voxelVector, const VoxelMesh& voxelMesh, const double lengthTolerance, const double angleTolerance, const double minimumLineSegmentDistance = 1);  // Calculating rotational symmetries.
    };


    // Line segment graph class.
    class LineSegmentGraph {
    private:
        int rotation;                          // Rotation level (1, 2, ..., n).
        std::vector<Vertex> vertices;          // Graph vertices (line segments).
        std::vector<std::vector<bool>> edges;  // 2D matrix of edges.
    public:
        LineSegmentGraph(const std::vector<LineSegmentPair>& lineSegmentPairs, const unsigned int rotation);  // Main constructor.
        void buildGraphFromLineSegmentsPairs(const std::vector<LineSegmentPair>& lineSegmentPairs);           // Building a graph from line segment pairs.

        std::vector<Vertex> getVertices() const;          // Vertices getter.
        std::vector<std::vector<bool>> getEdges() const;  // Edges getter.
        int getRotation() const;                          // Rotation getter.
    };


    // Line segment graph vertex class.
    class Vertex {
        unsigned int index;
        bool checked;
        LineSegment lineSegment;
    public:
        Vertex(const unsigned int index, const LineSegment& lineSegment);

        unsigned int getIndex() const;               // Index getter.
        bool getChecked() const;                     // Checked getter.
        void setChecked();                           // Checked setter (true).
        LineSegment getLineSegment() const;          // Line segment getter.
        LineSegment getInvertedLineSegment() const;  // Inverted line segment getter.
    };


    // Circumcircle class.
    class Circumcircle {
    public:
        unsigned int rotation;
        double radius;
        Point center;
        std::vector<std::pair<LineSegment, double>> lineSegmentsWithAngles;
        std::vector<bool> lineSegmentsInCircumcircle;

        Circumcircle(const unsigned int rotation, const Point& center, double radius, const std::vector<bool>& lineSegmentsInCircumcircle, const std::vector<LineSegment>& lineSegments);  // Main constructor.
        void addLineSegment(const LineSegment& ls);  // Adding a new line segment.
    };


    // Line segment pair class.
    class LineSegmentPair {
    public:
        LineSegment ls1;
        LineSegment ls2;

        LineSegmentPair(const LineSegment& ls1, const LineSegment& ls2);  // Main constructor.
    };

    // Line segment pair group class.
    class LineSegmentPairGroup {
    public:
        double distanceBetweenPairCenterPoints = -1.0;
        unsigned int rotation;
        unsigned int level;
        unsigned int layer;
        std::vector<LineSegmentPair> pairs;

        LineSegmentPairGroup(const std::vector<LineSegmentPair>& pairs, const unsigned int rotation, unsigned int level, unsigned int layer);  // Main constructor.
        void addPair(const LineSegmentPair& pair);                                                                                             // Adding a new pair to the group.
    };

    // Line segment pair subgroup class.
    class LineSegmentSubgroup {
    public:
        unsigned int rotation;
        std::vector<LineSegment> lineSegments;

        LineSegmentSubgroup(const std::vector<LineSegment>& lineSegments, const unsigned int rotation);  // Main constructor.
    };
};

#endif // ROTATIONALSYMMETRY_H
