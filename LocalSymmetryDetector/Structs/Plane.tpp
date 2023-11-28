#pragma once

#include <iomanip>
#include <string>

#include "Functions.hpp"
#include "Point.tpp"
#include "Voxel.hpp"


namespace Symmetry {
    /// <summary>
    /// Algebraic plane (ax + by + cz = d).
    /// </summary>
    /// <typeparam name="T"></typeparam>
    template <typename T>
    struct PlaneT {
        /******************************************************************************************************************************************************/
        /* STRUCT VARIABLES                                                                                                                                   */
        /******************************************************************************************************************************************************/

        T a;
        T b;
        T c;
        T d;



        /******************************************************************************************************************************************************/
        /* CONSTRUCTORS                                                                                                                                       */
        /******************************************************************************************************************************************************/
    private:
        PlaneT(const T a, const T b, const T c, const T d) : a(a), b(b), c(c), d(d) {}

    public:
        /// <summary>
        /// Default constructor.
        /// </summary>
        PlaneT() {
            a = b = c = d = std::numeric_limits<T>::max();
        }

        /// <summary>
        /// Constructor that builds a plane with two vectors and a point on it.
        /// </summary>
        /// <param name="p">: point that lies on the plane</param>
        /// <param name="v1">: first vector</param>
        /// <param name="v2">: second vector</param>
        PlaneT(const Point& p, const Vector3d& v1, const Vector3d& v2) {
            // Calculating a, b, c and d coefficients.
            // According to Linear Algebra (TM), (a, b, c) represents the plane normal vector.
            const Vector3d normal = VectorFunctions::crossProduct(v1, v2);
            a = static_cast<T>(normal.x);
            b = static_cast<T>(normal.y);
            c = static_cast<T>(normal.z);
            d = static_cast<T>(p.dot(normal));

            // To make our lives easier, A coefficient should always be >=0.
            if (a < 0) {
                a = static_cast<T>(-a);
                b = static_cast<T>(-b);
                c = static_cast<T>(-c);
                d = static_cast<T>(-d);
            }
        }

        /// <summary>
        /// Average plane constructor.
        /// </summary>
        /// <param name="plane1">: first plane</param>
        /// <param name="plane2">: second plane</param>
        PlaneT(const PlaneT<T>& plane1, const PlaneT<T>& plane2) {
            //if (VectorFunctions::angle(Vector3d(plane1.a, plane1.b, plane1.c), Vector3d(plane2.a, plane2.b, plane2.c)) > PI / 2) {
            //    a = (plane1.a - plane2.a) / 2.0;
            //    b = (plane1.b - plane2.b) / 2.0;
            //    c = (plane1.c - plane2.c) / 2.0;
            //    d = (plane1.d - plane2.d) / 2.0;
            //}
            //else {
                a = (plane1.a + plane2.a) / 2.0;
                b = (plane1.b + plane2.b) / 2.0;
                c = (plane1.c + plane2.c) / 2.0;
                d = (plane1.d + plane2.d) / 2.0;
            //}

            Vector3d v(a, b, c);
            //Point p = v.normalize();
            //a = p.x;
            //b = p.y;
            //c = p.z;

            // To make our lives easier, A coefficient should always be >=0.
            if (a < 0) {
                a = static_cast<T>(-a);
                b = static_cast<T>(-b);
                c = static_cast<T>(-c);
                d = static_cast<T>(-d);
            }
        }



        /******************************************************************************************************************************************************/
        /* STRUCT METHODS                                                                                                                                     */
        /******************************************************************************************************************************************************/

        /// <summary>
        /// Operator for equality check between the current and a given plane.
        /// </summary>
        /// <param name="p">: the given plane</param>
        /// <returns>True if the planes are equal, false otherwise</returns>
        bool operator == (const PlaneT<T>& p) const {
            return (
                (a == p.a || a == -p.a) &&
                (b == p.b || b == -p.b) &&

                (c == p.c || c == -p.c) &&
                (d == p.d || d == -p.d)
                );
        }

        /// <summary>
        /// Transforming plane to a string.
        /// </summary>
        /// <returns>Plane basic description as string</returns>
        std::string toString() const {
            std::stringstream ss;
            ss << std::setprecision(2);  // Setting the maximum of 2 decimal places.

            // A coefficient output.
            ss << a << "x";

            // B coefficient output.
            if (b >= 0) {
                ss << "+" << std::abs(b) << "y";
            }
            else {
                ss << b << "y";
            }

            // C coefficient output.
            if (c >= 0) {
                ss << "+" << std::abs(c) << "z=";
            }
            else {
                ss << c << "z=";
            }

            // D coefficient output.
            if (d >= 0) {
                ss << std::abs(d);
            }
            else {
                ss << d;
            }

            return ss.str();
        }

        /// <summary>
        /// Calculation of a vector that represents a normal to the plane.
        /// </summary>
        /// <returns>Normal vector</returns>
        Vector3d normalVector() const {
            return Vector3d(a, b, c);
        }

        /// <summary>
        /// Calculation of a vector that represents a parallel vector to the plane (along x and y axes).
        /// </summary>
        /// <returns>Parallel vector</returns>
        Vector3d parallelVector() const {
            Vector3d v(-b, a, 0);  // Calculating a parallel vector.

            // To make our lives easier again, the X component
            // of the vector should be >=0. Always. Immer und ewig.
            if (v.x < 0) {
                v.x = -v.x;
                v.y = -v.y;
            }

            return v;
        }

        /// <summary>
        /// Calculation of a vector that represents a parallel vector to the plane (along z axis).
        /// </summary>
        /// <returns>Parallel vector</returns>
        Vector3d zVector() const {
            return Vector3d(0, 0, 1);
        }

        /// <summary>
        /// Calculation of the plane intersection with the bounding box, where Y==max (if possible).
        /// </summary>
        /// <param name="vm">: voxel mesh</param>
        /// <returns>Pair of intersection point and its validity</returns>
        std::pair<Point, bool> calculateTopBoundingBoxIntersection(const VoxelMesh& vm) const {
            const Vector3d v = parallelVector();  // Calculating a parallel vector.

            // If the plane is constant along the Y axis, it
            // will never reach the top side of the bounding box.
            if (v.y == 0 && d != 0) {
                return { Point(), false };
            }

            // Calculating the coordinates of the intersection point.
            const double x = (d - b * vm.maxY) / a;
            const double y = vm.maxY;
            const double z = 0;

            // If X lies outside the bounding box, there is no intersection.
            if (x < vm.minX || x > vm.maxX) {
                return { Point(), false };
            }

            return { Point(x, y, z), true };
        }

        /// <summary>
        /// Calculation of the plane intersection with the bounding box, where Y==min (if possible).
        /// </summary>
        /// <param name="vm">: voxel mesh</param>
        /// <returns>Pair of intersection point and its validity</returns>
        std::pair<Point, bool> calculateBottomBoundingBoxIntersection(const VoxelMesh& vm) const {
            const Vector3d v = parallelVector();  // Calculating a parallel vector.

            // If the plane is constant along the Y axis, it
            // will never reach the bottom side of the bounding box.
            if (v.y == 0 && d != 0) {
                return { Point(), false };
            }

            // Calculating the coordinates of the intersection point.
            const double x = (d - b * vm.minY) / a;
            const double y = vm.minY;
            const double z = 0;

            // If X lies outside the bounding box, there is no intersection.
            if (x < vm.minX || x > vm.maxX) {
                return { Point(), false };
            }

            return { Point(x, y, z), true };
        }

        /// <summary>
        /// Calculation of the plane intersection with the bounding box, where X==min (if possible).
        /// </summary>
        /// <param name="vm">: voxel mesh</param>
        /// <returns>Pair of intersection point and its validity</returns>
        std::pair<Point, bool> calculateLeftBoundingBoxIntersection(const VoxelMesh& vm) const {
            const Vector3d v = parallelVector();  // Calculating a parallel vector.

            // If the plane is constant along the X axis, it
            // will never reach the left side of the bounding box.
            if (v.x == 0 && d != 0) {
                return { Point(), false };
            }

            // Calculating the coordinates of the intersection point.
            const double x = vm.minX;
            const double y = (d - a * vm.minX) / b;
            const double z = 0;

            // If Y lies outside the bounding box, there is no intersection.
            if (y < vm.minY || y > vm.maxY) {
                return { Point(), false };
            }

            return { Point(x, y, z), true };
        }

        /// <summary>
        /// Calculation of the plane intersection with the bounding box, where X==max (if possible).
        /// </summary>
        /// <param name="vm">: voxel mesh</param>
        /// <returns>Pair of intersection point and its validity</returns>
        std::pair<Point, bool> calculateRightBoundingBoxIntersection(const VoxelMesh& vm) const {
            // The Y component of the point is calculated as d/b.
            const Vector3d v = parallelVector();  // Calculating a parallel vector.

            // If the plane is constant along the X axis, it
            // will never reach the left side of the bounding box.
            if (v.x == 0 && d != 0) {
                return { Point(), false };
            }

            // Calculating the coordinates of the intersection point.
            const double x = vm.maxX;
            const double y = (d - a * vm.maxX) / b;
            const double z = 0;

            // If Y lies outside the bounding box, there is no intersection.
            if (y < vm.minY || y > vm.maxY) {
                return { Point(), false };
            }

            return { Point(x, y, z), true };
        }

        /// <summary>
        /// Determining whether a point lies on the left side of the plane.
        /// </summary>
        /// <param name="p">: the given point</param>
        /// <param name="vm">: voxel mesh</param>
        /// <returns>True if point is on the left side, false otherwise</returns>
        bool isPointOnTheLeftSide(const Point& p, const VoxelMesh& vm) const {
            // Calculating top, bottom and left bounding box intersections.
            // Since all planes intersect with the bounding box at least
            // twice, 3 intersection points are enough for sure.
            const auto [pTop, topValid] = calculateTopBoundingBoxIntersection(vm);           // Calculating top intersection.
            const auto [pBottom, bottomValid] = calculateBottomBoundingBoxIntersection(vm);  // Calculating bottom intersection.
            const auto [pLeft, leftValid] = calculateLeftBoundingBoxIntersection(vm);        // Calculating left intersection.
            double cross = 0.0;                                                              // Cross product between two vectors.

            // Calculation of the side with a help of the top intersection point.
            const Vector3d parallel = parallelVector();
            if (topValid) {
                const Point p1 = pTop - parallel;
                const Vector3d v1 = (Vector3d(parallelVector().x, parallelVector().y, 0)).normalize();
                const Vector3d v2 = (Vector3d(p.x - p1.x, p.y - p1.y, 0)).normalize();
                cross = v1.x * v2.y - v2.x * v1.y;
            }
            // Calculation of the side with a help of the bottom intersection point.
            else if (bottomValid) {
                const Point p1 = pBottom - parallel;
                const Vector3d v1 = (Vector3d(parallelVector().x, parallelVector().y, 0)).normalize();
                const Vector3d v2 = (Vector3d(p.x - p1.x, p.y - p1.y, 0)).normalize();
                cross = v1.x * v2.y - v2.x * v1.y;
            }
            // Calculation of the side with a help of the left intersection point.
            else if (leftValid) {
                const Point p1 = pLeft - parallel;
                const Vector3d v1 = (Vector3d(parallelVector().x, parallelVector().y, 0)).normalize();
                const Vector3d v2 = (Vector3d(p.x - p1.x, p.y - p1.y, 0)).normalize();
                cross = v1.x * v2.y - v2.x * v1.y;
            }

            // If cross product is <0, the point lies on the left side of the plane.
            return cross < 0;
        }

        /// <summary>
        /// Calculation of a projection point on the plane of the given point.
        /// </summary>
        /// <param name="p">: the given point</param>
        /// <param name="vm">: voxel mesh</param>
        /// <returns>Projection point</returns>
        Point calculateProjectionPoint(const Point& p, const VoxelMesh& vm) const {
            Point projection;

            // Calculating top, bottom and leftbounding box intersections.
            const auto [pTop, topValid] = calculateTopBoundingBoxIntersection(vm);           // Calculating the top intersection.
            const auto [pBottom, bottomValid] = calculateBottomBoundingBoxIntersection(vm);  // Calculating the bottom intersection.
            const auto [pLeft, leftValid] = calculateLeftBoundingBoxIntersection(vm);        // Calculating the left intersection.
            const auto [pRight, rightValid] = calculateRightBoundingBoxIntersection(vm);     // Calculating the right intersection.

            // Calculation of the projection point with a help of the top intersection point.
            if (topValid) {
                // Base vector calculation.
                const Vector3d v1(parallelVector());
                const Vector3d v2(pTop.x - p.x + vm.minX, pTop.y - p.y + vm.minY, 0);
                const Vector3d vn = Vector3d(pTop.x - p.x + vm.minX, pTop.y - p.y + vm.minY, 0).normalize();
                const double length = v2.length();

                // Projection calculation.
                const double dot = -v1.dot(vn);
                projection.x = pTop.x + dot * v1.x * length + vm.minX;
                projection.y = pTop.y + dot * v1.y * length + vm.minX;
                projection.z = p.z;
            }
            // Calculation of the projection point with a help of the bottom intersection point.
            else if (bottomValid) {
                // Base vector calculation.
                const Vector3d v1(parallelVector());
                const Vector3d v2(pBottom.x - p.x + vm.minX, pBottom.y - p.y + vm.minY, 0);
                const Vector3d vn = Vector3d(pBottom.x - p.x + vm.minX, pBottom.y - p.y + vm.minY, 0).normalize();
                const double length = v2.length();

                // Projection calculation.
                const double dot = -v1.dot(vn);
                projection.x = pBottom.x + dot * v1.x * length + vm.minX;
                projection.y = pBottom.y + dot * v1.y * length + vm.minY;
                projection.z = p.z;
            }
            // Calculation of the projection point with a help of the left intersection point.
            else if (leftValid) {
                // Base vector calculation.
                const Vector3d v1(parallelVector());
                const Vector3d v2(pLeft.x - p.x + vm.minX, pLeft.y - p.y + vm.minY, 0);
                const Vector3d vn = Vector3d(pLeft.x - p.x + vm.minX, pLeft.y - p.y + vm.minY, 0).normalize();
                const double length = v2.length();

                // Projection calculation.
                const double dot = std::abs(v1.dot(vn));
                projection.x = pLeft.x + dot * v1.x * length + vm.minX;
                projection.y = pLeft.y + dot * v1.y * length + vm.minY;
                projection.z = p.z;
            }
            // Calculation of the projection point with a help of the right intersection point.
            else if (rightValid) {
                // Base vector calculation.
                const Vector3d v1(parallelVector());
                const Vector3d v2(pRight.x - p.x + vm.minX, pRight.y - p.y + vm.minY, 0);
                const Vector3d vn = Vector3d(pRight.x - p.x + vm.minX, pRight.y - p.y + vm.minY, 0).normalize();
                const double length = v2.length();

                // Projection calculation.
                const double dot = std::abs(v1.dot(vn));
                projection.x = pRight.x - dot * v1.x * length + vm.minX;
                projection.y = pRight.y - dot * v1.y * length + vm.minY;
                projection.z = p.z;
            }

            return projection;
        }

        /// <summary>
        /// Calculation of Euclidean distance between the plane and the given point.
        /// </summary>
        /// <param name="p">: given point</param>
        /// <returns>Distance</returns>
        double distanceToPoint(const Point& p) const {
            return std::abs(a * p.x + b * p.y + c * p.z - d) / std::sqrt(std::pow(a, 2) + std::pow(b, 2) + std::pow(c, 2));
        }

        /// <summary>
        /// Calculation of a point that lies across the plane (according to the given point).
        /// </summary>
        /// <param name="p">: the given point</param>
        /// <param name="vm">: voxel mesh</param>
        /// <returns>Point across the plane</returns>
        Point calculatePointAcrossPlane(const Point& p, const VoxelMesh& vm) const {
            const Point projectionPoint = calculateProjectionPoint(p, vm);  // Calculation of the projection point.
            const Vector3d vector = projectionPoint - p;                    // Calculation of a vector from the point to the projection point.
            const Point acrossPoint = projectionPoint + vector;             // Calculation of the opposite point.

            return acrossPoint;
        }

        /// <summary>
        /// Calculation of the Y coordinate on the plane according to the given X coordinate.
        /// </summary>
        /// <param name="x">: X coordinate</param>
        /// <returns>Y coordinate</returns>
        double calculateYCoordinateAtX(const double x) const {
            // Y can be calculated pretty easily from the plane equation (ax+by+cz=d).
            const double y = (d - a * x) / b;
            return y;
        }

        /// <summary>
        /// Calculation of the X coordinate on the plane according to the given Y coordinate.
        /// </summary>
        /// <param name="y">: Y coordinate</param>
        /// <returns>X coordinate</returns>
        double calculateXCoordinateAtY(const double y) const {
            // X can be calculated pretty easily from the plane equation (ax+by+cz=d).
            const double x = (d - b * y) / a;
            return x;
        }

        /// <summary>
        /// Detection of indices of the points that lie on the symmetry plane and voxel edge simultaneously.
        /// </summary>
        /// <param name="points">: vector of points (point cloud)</param>
        /// <param name="vm">: voxel mesh</param>
        /// <returns>Vector of indices</returns>
        std::vector<uint> getPointsIndicesOnPlaneAndVoxelEdge(const std::vector<Point>& points, const VoxelMesh& vm) const {
            std::vector<uint> indices;

            for (uint i = 0; i < points.size(); i++) {
                // If the projection point and the point do not have
                // the same coordinates, the point does not lie on the plane.
                const Point projection = calculateProjectionPoint(points[i], vm);
                if (points[i] != projection) {
                    continue;
                }

                // If the quotient between the point X coordinate
                if (!Tolerance::isInTolerance(std::fmod(points[i].x / vm.sideX, 1), 0.0, 0.0001) &&
                    !Tolerance::isInTolerance(std::fmod(points[i].y / vm.sideY, 1), 0.0, 0.0001)
                )
                {
                    continue;
                }

                indices.push_back(i);
            }

            return indices;
        }

        /// <summary>
        /// Calculation of the desired start point of the plane in a symmetry.
        /// </summary>
        /// <param name="vm">: voxel mesh</param>
        /// <returns>Tuple(start point, plane distance, plane distance by X, plane distance by Y)</returns>
        std::tuple<Point, double, double, double> calculateStartPoint(const VoxelMesh& vm) const {
            // Calculating the intersections betweeen the symmetry plane and the bounding box.
            const auto [pTop, topValid] = calculateTopBoundingBoxIntersection(vm);
            const auto [pBottom, bottomValid] = calculateBottomBoundingBoxIntersection(vm);
            const auto [pLeft, leftValid] = calculateLeftBoundingBoxIntersection(vm);
            const auto [pRight, rightValid] = calculateRightBoundingBoxIntersection(vm);

            // Setting the position of a plane.
            Point position;
            const double RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE = 0.5;
            const double planeY = parallelVector().y;
            double distance = 0.0;
            double distanceX = 0.0;
            double distanceY = 0.0;

            // Positioning the plane according to all the posibilities (8).
            if (planeY < 0) {
                // Option 1: left and bottom bounding box intersections.
                if (leftValid && bottomValid) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pBottom.x - pLeft.x;
                    distanceY = pLeft.y - pBottom.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pLeft.x - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX;
                    position.y = pLeft.y + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 2: top and bottom bounding box intersections.
                else if (topValid && bottomValid) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pBottom.x - pTop.x;
                    distanceY = pTop.y - pBottom.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pTop.x - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX;
                    position.y = pTop.y + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 3: top and right bounding box intersections.
                else if (topValid && rightValid) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pRight.x - pTop.x;
                    distanceY = pTop.y - pRight.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pTop.x - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX;
                    position.y = pTop.y + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 4: top and right bounding box intersections.
                else if (leftValid && rightValid) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pRight.x - pLeft.x;
                    distanceY = pLeft.y - pRight.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pLeft.x - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX;
                    position.y = pLeft.y + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
            }
            else {
                // Option 5: bottom and right bounding box intersections.
                if (bottomValid && rightValid) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pRight.x - pBottom.x;
                    distanceY = pRight.y - pBottom.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pRight.x + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX - 2 * distanceX;
                    position.y = pRight.y - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY - distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 6: bottom and top bounding box intersections.
                else if (bottomValid && topValid) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pTop.x - pBottom.x;
                    distanceY = pTop.y - pBottom.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pBottom.x + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX - distanceX;
                    position.y = pBottom.y - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 7: left and top bounding box intersections.
                else if (leftValid && topValid) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pTop.x - pLeft.x;
                    distanceY = pTop.y - pLeft.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pLeft.x + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX - distanceX;
                    position.y = pLeft.y - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 8: left and right bounding box intersections.
                else if (leftValid && rightValid) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pRight.x - pLeft.x;
                    distanceY = pRight.y - pLeft.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pLeft.x + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX - distanceX;
                    position.y = pLeft.y - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
            }

            return std::make_tuple(position, distance, distanceX, distanceY);
        }

        /// <summary>
        /// Calculation of a parallel plane to the current one.
        /// </summary>
        /// <param name="distance"></param>
        /// <returns></returns>
        PlaneT<T> parallelPlane(const double distance) const {
            return PlaneT<T>(a, b, c, d + distance);
        }
    };



    /******************************************************************************************************************************************************/
    /* ALIASES                                                                                                                                            */
    /******************************************************************************************************************************************************/

    using Plane = PlaneT<double>;
};