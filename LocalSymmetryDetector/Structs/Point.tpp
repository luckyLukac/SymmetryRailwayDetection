#pragma once
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "Constants.hpp"
#include "Tolerance.hpp"


namespace Symmetry {
    /// <summary>
    /// Templatized point struct with X, Y and Z coordinates.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    template <class T>
    struct PointT {
        /******************************************************************************************************************************************************/
        /* STRUCT VARIABLES                                                                                                                                   */
        /******************************************************************************************************************************************************/

        T x;
        T y;
        T z;



        /******************************************************************************************************************************************************/
        /* CONSTRUCTORS                                                                                                                                       */
        /******************************************************************************************************************************************************/

        /// <summary>
        /// Default constructor.
        /// </summary>
        PointT() :
            x(static_cast<T>(0)),
            y(static_cast<T>(0)),
            z(static_cast<T>(0)) {}

        /// <summary>
        /// Constructor with X, Y and Z coordinates.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        PointT(const T x, const T y, const T z) :
            x(x),
            y(y),
            z(z) {}

        /// <summary>
        /// Copy constructor.
        /// </summary>
        /// <param name="p"></param>
        PointT(const PointT<T>& p) :
            x(p.x),
            y(p.y),
            z(p.z) {}



        /******************************************************************************************************************************************************/
        /* STRUCT METHODS                                                                                                                                     */
        /******************************************************************************************************************************************************/

        /// <summary>
        /// Operator for adding a scalar to the point.
        /// </summary>
        /// <param name="number">: scalar</param>
        /// <returns>Updated point</returns>
        PointT<T> operator + (const double number) const {
            return PointT<T>(x + number, y + number, z + number);
        }

        /// <summary>
        /// Operator for adding a vector to the point.
        /// </summary>
        /// <param name="p">: vector</param>
        /// <returns>Updated point</returns>
        PointT<T> operator + (const PointT<T>& p) const {
            return PointT<T>(x + p.x, y + p.y, z + p.z);
        }

        /// <summary>
        /// Operator for adding a vector to the current point.
        /// </summary>
        /// <param name="p">: vector</param>
        /// <returns>Updated point</returns>
        PointT<T> operator += (const PointT<T>& p) const {
            return PointT<T>(x + p.x, y + p.y, z + p.z);
        }

        /// <summary>
        /// Operator for substracting a scalar from the current point.
        /// </summary>
        /// <param name="number">: scalar</param>
        /// <returns>Updated point</returns>
        PointT<T> operator - (const double number) const {
            return PointT<T>(x - number, y - number, z - number);
        }

        /// <summary>
        /// Operator for substracting a vector from the current point.
        /// </summary>
        /// <param name="p">: vector</param>
        /// <returns>Updated point</returns>
        PointT<T> operator - (const PointT<T>& p) const {
            return PointT<T>(x - p.x, y - p.y, z - p.z);
        }

        /// <summary>
        /// Operator for multiplying a scalar with the current point.
        /// </summary>
        /// <param name="factor">: scalar</param>
        /// <returns>Updated point</returns>
        PointT<T> operator * (const double factor) const {
            return PointT<T>(factor * x, factor * y, factor * z);
        }

        /// <summary>
        /// Operator for multiplying a vector with the current point.
        /// </summary>
        /// <param name="p">: vector</param>
        /// <returns>Updated point</returns>
        PointT<T> operator * (const PointT<T>& p) const {
            return PointT<T>(p.x * x, p.y * y, p.z * z);
        }

        /// <summary>
        /// Operator for dividing the current point with a scalar.
        /// </summary>
        /// <param name="factor">: scalar</param>
        /// <returns>Updated point</returns>
        PointT<T> operator / (const double factor) const {
            return PointT<T>(x / factor, y / factor, z / factor);
        }

        /// <summary>
        /// Operator for dividing the current point with a scalar.
        /// </summary>
        /// <param name="factor">: scalar</param>
        /// <returns>Updated point</returns>
        PointT<T> operator /= (const double factor) {
            return PointT<T>(x / factor, y / factor, z / factor);
        }

        /// <summary>
        /// Equality check of the current and the given point.
        /// </summary>
        /// <param name="p">: the given point</param>
        /// <returns>True if both points are equal, false otherwise</returns>
        bool operator == (const PointT<T>& p) const {
            return Tolerance::isInTolerance(x, p.x, 0.01) && Tolerance::isInTolerance(y, p.y, 0.01) && Tolerance::isInTolerance(z, p.z, 0.01);
        }

        /// <summary>
        /// Inquality check of the current and the given point.
        /// </summary>
        /// <param name="p">: the given point</param>
        /// <returns>True if both points are not equal, false otherwise</returns>
        bool operator != (const PointT<T>& p) const {
            return !(*this == p);
        }

        /// <summary>
        /// Transforming a point to a string.
        /// </summary>
        /// <returns>Point basic description as string</returns>
        std::string toString() const {
            std::stringstream ss;
            ss << "Point(" << x << "," << y << "," << z << ")";
            return ss.str();
        }

        /// <summary>
        /// Dot product of the current vector (point) and the given one.
        /// </summary>
        /// <param name="p">: the given vector</param>
        /// <returns>Dot product (scalar value)</returns>
        double dot(const PointT<T>& p) const {
            const T xDot = p.x * x;
            const T yDot = p.y * y;
            const T zDot = p.z * z;

            return xDot + yDot + zDot;
        }

        /// <summary>
        /// Calculation of a vector (point) length.
        /// </summary>
        /// <returns>Length of the vector</returns>
        double length() const {
            const T xDot = x * x;
            const T yDot = y * y;
            const T zDot = z * z;
            const T dot = xDot + yDot + zDot;

            return std::sqrt(dot);
        }

        /// <summary>
        /// Normalization of a vector.
        /// </summary>
        /// <returns>Normalized vector</returns>
        PointT<T> normalize() const {
            const double len = length();
            return PointT<T>(x / len, y / len, z / len);
        }

        /// <summary>
        /// Convertion of the point to a certain C++ type.
        /// </summary>
        /// <typeparam name="U"></typeparam>
        /// <returns>Converted point</returns>
        template <typename U>
        PointT<U> convert() const {
            return PointT<U>(static_cast<U>(x), static_cast<U>(y), static_cast<U>(z));
        }
    };



    /******************************************************************************************************************************************************/
    /* ALIASES                                                                                                                                            */
    /******************************************************************************************************************************************************/

    using IntPoint = PointT<int>;
    using Point = PointT<double>;
    using Vector3f = PointT<float>;
    using Vector3d = PointT<double>;



    /******************************************************************************************************************************************************/
    /* ENUM CLASSES                                                                                                                                       */
    /******************************************************************************************************************************************************/

    /// <summary>
    /// Enumeration class for point position according to the symmetry plane or the symmetry axis.
    /// </summary>
    enum class PointPosition {
        left,
        right,
        center,
        rotational,
        notInSymmetry,
        railwayCandidate,
        asymmetry
    };
};