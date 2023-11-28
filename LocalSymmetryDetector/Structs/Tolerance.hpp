#pragma once


namespace Symmetry {
    /// <summary>
    /// Namespace for difference related functions.
    /// </summary>
    namespace Difference {
        /// <summary>
        /// Calculation of the difference between two values.
        /// </summary>
        /// <param name="v1">: first value</param>
        /// <param name="v2">: second value</param>
        /// <returns>Calculated difference</returns>
        inline double difference(const double v1, const double v2) {
            return std::abs(v1 - v2);
        }
    };

    /// <summary>
    /// Namespace for tolerance related functions.
    /// </summary>
    namespace Tolerance {
        /// <summary>
        /// Calculation whether two floating point values are equal (inside of the given tolerance).
        /// </summary>
        /// <param name="v1">: first value</param>
        /// <param name="v2">: second value</param>
        /// <param name="tolerance">: tolerance</param>
        /// <returns>True if the two values differ for a value, smaller than allowed tolerance; false otherwise</returns>
        inline bool isInTolerance(const double v1, const double v2, const double tolerance) {
            if (std::abs(v1 - v2) <= tolerance) {
                return true;
            }

            return false;
        }
    };
};