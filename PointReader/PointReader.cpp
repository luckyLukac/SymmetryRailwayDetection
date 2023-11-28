#include <fstream>

#include "PointReader.hpp"


using namespace Symmetry;


/// <summary>
/// Point classification type enum.
/// </summary>
enum PointType {
    unclassified = 1,
    ground = 2,
    vegetation = 5,
    roofs = 6
};


std::vector<Point> PointReader::readPointsFromPointCloud(const std::string path) {
    // Opening the file.
    std::ifstream file(path);
    std::vector<Point> points;
    std::string str;

    const std::string fileType = path.substr(path.rfind('.') + 1);

    if (fileType == "pc") {
        // Reading the number of points.
        std::getline(file, str);
        const uint count = std::stoi(str);

        // Reading each point.
        for (uint i = 0; i < count; i++) {
            // Reading point coordinates.
            std::getline(file, str, ' ');
            const double x = std::stod(str);
            std::getline(file, str, ' ');
            const double y = (std::stod(str));
            std::getline(file, str);
            const double z = std::stod(str);

            // Translating the point to reference point local coordinate system.
            Point point = Point(x, y, z);
            Point referencePoint(0, 0, 0);
            const double angle = 0.00;
            const double xTranslated = point.x - referencePoint.x;
            const double yTranslated = point.y - referencePoint.y;

            // Rotating the point.
            const double rotatedX = xTranslated * std::cos(angle) - yTranslated * std::sin(angle);
            const double rotatedY = xTranslated * std::sin(angle) + yTranslated * std::cos(angle);

            // Translating the point back to the global coordinate system.
            const double newX = rotatedX + referencePoint.x;
            const double newY = rotatedY + referencePoint.y;

            points.push_back(Point(newX, newY, z) * 100);
        }
    }
    else if (fileType == "csv") {
        // Reading point coordinates.
        while (std::getline(file, str, ',')) {
            const double x = std::stod(str);
            std::getline(file, str, ',');
            const double z = -(std::stod(str));
            std::getline(file, str);
            const double y = std::stod(str);

            // Translating the point to reference point local coordinate system.
            Point point = Point(x, y, z);
            Point referencePoint(0, 0, 0);
            const double angle = 0 * PI / 180;
            const double xTranslated = point.x - referencePoint.x;
            const double yTranslated = point.y - referencePoint.y;

            // Rotating the point.
            const double rotatedX = xTranslated * std::cos(angle) - yTranslated * std::sin(angle);
            const double rotatedY = xTranslated * std::sin(angle) + yTranslated * std::cos(angle);

            // Translating the point back to the global coordinate system.
            const double newX = rotatedX + referencePoint.x;
            const double newY = rotatedY + referencePoint.y;

            points.push_back(Point(newX, newY, z) * 100);
        }

        points.push_back(Point(-1, -1, 0));
        points.push_back(Point(0, 0, 5));
    }

    return points;
}
