#ifndef POINTREADER_H
#define POINTREADER_H

#include <vector>

#include "LocalSymmetryDetector/Structs/Point.tpp"


namespace Symmetry {
    namespace PointReader {
        /// <summary>
        /// Reading points from a point cloud file.
        /// </summary>
        /// <param name="path">: path to the data with the point cloud </param>
        /// <returns> vector of points </returns>
        std::vector<Point> readPointsFromPointCloud(const std::string path);
    };
};

#endif // POINTREADER_H
