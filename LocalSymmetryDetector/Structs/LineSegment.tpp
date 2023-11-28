#pragma once

#include <sstream>
#include <string>
#include <vector>

#include "Functions.hpp"
#include "Point.tpp"


namespace Symmetry {
    /// <summary>
    /// Line segment, bounded by two endpoints.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    template <typename T>
    struct LineSegmentT {
        /******************************************************************************************************************************************************/
        /* STRUCT VARIABLES                                                                                                                                   */
        /******************************************************************************************************************************************************/
        Point p1;      // First endpoint that limits the line segment.
        Point p2;      // Second endpoint that limits the line segment.
       


        /******************************************************************************************************************************************************/
        /* CONSTRUCTORS                                                                                                                                       */
        /******************************************************************************************************************************************************/
        
        /// <summary>
        /// Main constructor of the line segment.
        /// </summary>
        LineSegmentT() :
            p1(Point()),
            p2(Point())
        {}

        /// <summary>
        /// Constructor with two endpoints.
        /// </summary>
        /// <param name="p1">: first endpoint</param>
        /// <param name="p2">: second endpoint</param>
        LineSegmentT(const Point& p1, const Point& p2) :
            p1(p1),
            p2(p2)
        {}



        /******************************************************************************************************************************************************/
        /* STRUCT METHODS                                                                                                                                     */
        /******************************************************************************************************************************************************/
        
        /// <summary>
        /// Calculation of line segment length.
        /// </summary>
        /// <returns>Length of the line segment</returns>
        T length() const {
            return static_cast<T>(std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2) + std::pow(p1.z - p2.z, 2)));
        }

        /// <summary>
        /// Transforming line segment to a string.
        /// </summary>
        /// <returns>Line segment basic description as string</returns>
        std::string toString() const {
            std::stringstream ss;
            ss << "LineSegment(" << p1.toString() << ", " << p2.toString() << ")";
            return ss.str();
        }

        /// <summary>
        /// Operator for equality check between current and given line segments.
        /// </summary>
        /// <param name="ls">: the given line segment</param>
        /// <returns>True if both line segments are equal, false otherwise</returns>
        bool operator == (const LineSegmentT& ls) const {
            return (
                (p1 == ls.p1 && p2 == ls.p2) ||
                (p1 == ls.p2 && p2 == ls.p1)
            );
        }

        /// <summary>
        /// Calculation of the angle between current and given line segments.
        /// </summary>
        /// <param name="ls">: the given line segment</param>
        /// <returns>Angle in radians</returns>
        double angle(const LineSegmentT& ls) const {
            const Vector3d v1 = p2 - p1;
            const Vector3d v2 = ls.p2 - ls.p1;

            return VectorFunctions::angle(v1, v2);
        }

        /// <summary>
        /// Calculation of the center point of the line segment.
        /// </summary>
        /// <returns>Midpoint</returns>
        Point midpoint() const {
            return Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2);
        }

        /// <summary>
        /// Calculation of the perpendicular line segment to the current one that passes through the midpoint of the current line segment.
        /// </summary>
        /// <returns>Long perpendicular line segment</returns>
        LineSegmentT<T> perpendicularLineSegment() const {
            const Vector3d lsVector = p2 - p1;                                                              // Calculating a line segment vector.
            const Vector3d perpendicular(-lsVector.y, lsVector.x, 0);                                       // Calculating a perpendicular vector to the line segment.
            const Point center = midpoint();                                                                // Calculating the center point of the line segment.
            const LineSegmentT<T> bisection(center - perpendicular * 1000, center + perpendicular * 1000);  // Calculating the bisection.

            return bisection;
        }
    };



    /******************************************************************************************************************************************************/
    /* ALIASES                                                                                                                                            */
    /******************************************************************************************************************************************************/

    using LineSegment = LineSegmentT<double>;
    using SplitLineSegments = std::vector<std::vector<LineSegment>>;



    /******************************************************************************************************************************************************/
    /* LINE SEGMENT FUNCTIONS                                                                                                                             */
    /******************************************************************************************************************************************************/

    /// <summary>
    /// Namespace for LineSegment related functions.
    /// </summary>
    namespace LineSegmentFunctions {
        /// <summary>
        /// Function for determining whether two line segments have an intersection point.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="ls1">: first line segment</param>
        /// <param name="ls2">: second line segment</param>
        /// <param name="tolerance">: tolerance (to compensate round-off errors)</param>
        /// <param name="edgesIncluded">: line segment pairs' common point (if existing) is considered an intersection if true</param>
        /// <returns>True if line segments intersect with each other, false otherwise</returns>
        template <typename T>
        bool doLineSegmentsIntersect(const LineSegmentT<T>& ls1, const LineSegmentT<T>& ls2, const double tolerance, const bool edgesIncluded = false) {
            // If the two line segments are parallel and at the same location,
            // they have infinite number of intersections. Yaaaaaayyyy!
            if (ls1 == ls2) {
                return true;
            }

            // Getting points from line segments.
            const Point P1 = ls1.p1;
            const Point P2 = ls1.p2;
            const Point P3 = ls2.p1;
            const Point P4 = ls2.p2;

            // Calculating vectors between points.
            const Vector3d V12 = P2 - P1;
            const Vector3d V34 = P4 - P3;
            const Vector3d V31 = P1 - P3;

            // Coefficient calculation.
            const T D = (V12.x * V34.y) - (V12.y * V34.x);
            const T A = (V34.x * V31.y) - (V34.y * V31.x);
            const T B = (V12.x * V31.y) - (V12.y * V31.x);

            // If D == 0, the line segments are parallel, therefore,
            // they do not have an intersection point.
            if (Tolerance::isInTolerance(D, 0.0, tolerance)) {
                return false;
            }

            // Normalizing coefficients.
            const T Ua = A / D;
            const T Ub = B / D;


            // If edges of line segments are not considered intersections and the two
            // line segments have a point in common, it is not considered an intersection.
            if (!edgesIncluded &&
                (Tolerance::isInTolerance(Ua, 0, 0.0001) ||
                 Tolerance::isInTolerance(Ua, 1, 0.0001) ||
                 Tolerance::isInTolerance(Ub, 0, 0.0001) ||
                 Tolerance::isInTolerance(Ub, 1, 0.0001)
                )
            )
            {
                return false;
            }

            // If the two line segment intersect,
            // true is returned as a result.
            if (Ua > -0.0001 && Ua < 1.0001 && Ub > -0.0001 && Ub < 1.0001) {
                return true;
            }

            return false;
        }

        /// <summary>
        /// Calculation of an intersection point of two line segments.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="ls1">: first line segment</param>
        /// <param name="ls2">: second line segment</param>
        /// <param name="tolerance">: tolerance (to compensate round-off errors)</param>
        /// <returns></returns>
        template <typename T>
        std::pair<Point, bool> calculateIntersection(const LineSegmentT<T>& ls1, const LineSegmentT<T>& ls2, const double tolerance = 0.0001) {
            // Getting points from line segments.
            const Point P1 = ls1.p1;
            const Point P2 = ls1.p2;
            const Point P3 = ls2.p1;
            const Point P4 = ls2.p2;

            // Calculating vectors between points.
            const Vector3d V12 = P2 - P1;
            const Vector3d V34 = P4 - P3;
            const Vector3d V31 = P1 - P3;

            // Coefficient calculation.
            const T D = (V12.x * V34.y) - (V12.y * V34.x);
            const T A = (V34.x * V31.y) - (V34.y * V31.x);
            const T B = (V12.x * V31.y) - (V12.y * V31.x);

            // If D == 0, the line segments are parallel, therefore,
            // they do not have an intersection point.
            if (Tolerance::isInTolerance(D, 0.0, tolerance)) {
                return { Point(), false };
            }

            // Normalizing coefficients.
            const T Ua = A / D;
            const T Ub = B / D;

            // If the two line segment intersect,
            // true is returned as a result.
            if (Ua >= 0 && Ua <= 1 && Ub >= 0 && Ub <= 1) {
                const T x = P1.x + Ua * (P2.x - P1.x);
                const T y = P1.y + Ua * (P2.y - P1.y);
                const T z = P1.z + Ua * (P2.z - P1.z);

                return { Point(x, y, z), true };
            }

            return { Point(), false };
        }
    };
};