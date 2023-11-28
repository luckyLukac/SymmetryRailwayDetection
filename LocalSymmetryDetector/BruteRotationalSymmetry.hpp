#ifndef BRUTEROTATIONALSYMMETRY_H
#define BRUTEROTATIONALSYMMETRY_H

#include <vector>

#include "Structs/Plane.tpp"
#include "Structs/Point.tpp"
#include "Structs/Voxel.hpp"
#include "LocalSymmetry.hpp"


namespace Symmetry {
    // Class for angle ranges (while searching rotational symmetries).
    struct AngleRangeEvent {
        // CLASS VARIABLES
        double angle;  // Angle in radians.
        int level;     // Level the event is part of (see documentation for more details).
        bool start;    // True if this is a starting event, false otherwise.
        Voxel voxel;   // Voxel of the angle event.

        // CONSTRUCTORS
        // Default constructor without arguments.
        AngleRangeEvent() :
            angle(0.0),
            level(0),
            start(false)
        {}

        // Constructor with all necessary arguments.
        AngleRangeEvent(const double angle, const int level, const bool start, const Voxel& v) :
            angle(angle),
            level(level),
            start(start),
            voxel(v)
        {}
    };


    // Class for rotational symmetries.
    class BruteRotationalSymmetry : public LocalSymmetry {
    private:
        // OBJECT VARIABLES
        int rotation;                 // Rotation level (1, 2, ..., n).
        int upperLayer;               // Index of the starting (upper) layer.
        int lowerLayer;               // Index of the ending (lower) layer.
        std::vector<Voxel> voxels;    // Vector of pointers to voxels, included in the symmetry.
        Point symmetryAxis;   // Voxel, where the symmetry axis is.

        // CLASS VARIABLES
        static inline std::vector<float> symmetryAngles;              // Angles of the rotational symmetry.
        static inline std::vector<std::vector<bool>> geometrySearch;  // Vector for voxel search in 2D (layer).
        static inline int32_t geometrySearchSizeX;                    // Dimension X size in the voxel field.
        static inline int32_t geometrySearchSizeY;                    // Dimension Y size in the voxel field.

        // PRIVATE OBJECT METHODS
        std::vector<BruteRotationalSymmetry> getRotationalSymmetriesAroundCenterPoint(std::vector<std::vector<Voxel>>& field, Point centerPoint, const int maxRadius, int pixelX, int pixelY);  // Function that increments the radius of the circumference around the center point with the step of 1. In every step it selects all interesting voxels the circumference overlaps with and calls the function FindRotationalSymmetries. In the end a union of all interesting voxels from obtained symmetries with different radii is made. The function returns a list of rotational symmetries (2 through 16), where each rotational symmetry is present once at most.
        std::vector<BruteRotationalSymmetry> mergeSymmetries(std::vector<BruteRotationalSymmetry>& symmetries);  // Merging symmetries by layers. Symmetries with the same search axis and rotation level are merged.
        std::vector<BruteRotationalSymmetry> getRotationalSymmetriesOnRadius(Point centerPoint, const std::vector<Voxel>& interestingVoxels);  // A main function for searching symmetries gets a center point and a list of interesting voxels. Symmetries are searched for only in the set of interesting voxels that are on the same radius from the center point (or inside the same interval).
    public:
        // CONSTRUCTOR AND DESTRUCTOR
        BruteRotationalSymmetry();   // Default constructor.
        BruteRotationalSymmetry(const Point& symmetryAxis, const int& rotation, const int& upperLayer, const int& lowerLayer);  // Constructor with all parameters.

        // OVERLOADED OPERATORS
        BruteRotationalSymmetry operator & (const BruteRotationalSymmetry& s) const;  // Operator & is used for merging two symmetries.
        BruteRotationalSymmetry operator &= (BruteRotationalSymmetry& s);             // Operator &= is used for merging the current symmetry with another.

        // GETTERS AND SETTERS
        int getVoxelCount() const;                         // Getting voxel count.
        std::vector<Voxel>& getVoxels();                 // Getting all voxels.
        Voxel getVoxel(unsigned int index) const;         // Getting a voxel at a certain index.
        void addVoxel(Voxel v);                           // Adding a voxel to the vector.
        Point getSymmetryAxis() const;             // Getting a symmetry axis.
        void setSymmetryAxis(Point symmetryAxis);  // Setting a symmetry axis.
        int getRotation() const;                       // Getting a rotation.
        void setRotation(int rotation);                // Setting a rotation.
        int getUpperLayer() const;                    // Getting an upper layer.
        void setUpperLayer(int upperLayer);           // Setting an upper layer.
        int getLowerLayer() const;                    // Getting a lower layer.
        void setLowerLayer(int lowerLayer);           // Setting a lower layer.

        // CLASS METHODS
        static void setGeometrySearch(int32_t x, int32_t y);    // Symmetry search class initialization method.
        static void clearGeometrySearch(int32_t x, int32_t y);  // Cleaning the history of geometry search.
        static void setSymmetryAngles();                        // Calculation of angles of rotational symmetries (interval [1, 16]).
        static void getDistanceRanges(Point centerPoint, Voxel v, double& minD, double& maxD);  // Calculation of the minimum and the maximum distance from the center point to the voxel (in 2D).
        static float getAngle(const Point a, const Point b);  // Getting the angle between two vectors (represented as Points).
        static void getAngleRanges(const Point& centerPoint, Voxel v, float& minAngle, float& maxAngle);  // Calculation of the minimum and the maximum angle between horizontal vector and the vectors from the center point to the vertices of the voxel (in 2D).

        // OBJECT METHODS
        std::vector<BruteRotationalSymmetry> getRotationalSymmetries(std::vector<std::vector<std::vector<Voxel>>> voxels, VoxelMesh& voxelMesh);  // Finding rotational symmetries.
        std::vector<PointPosition> calculatePositionsFromAxis(const Point* axis);  // Calculating the point positions according to the axis.
    };
};

#endif // BRUTEROTATIONALSYMMETRY_H
