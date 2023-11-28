#include <algorithm>
#include <iostream>
#include <set>
#include <stack>

#include "Structs/Tolerance.hpp"
#include "LocalSymmetry.hpp"

using namespace Symmetry;


const std::vector<Point>& LocalSymmetry::getPoints() const {
    return m_Points;
}

void LocalSymmetry::setPoints(std::vector<Point>& points, bool moveToCoordinateCenter) {
    m_Points = points;
    m_VoxelMesh = {};

    for (Point& lasPoint : LocalSymmetry::m_Points) {
        // Getting point coordinates.
        const double x = lasPoint.x;
        const double y = lasPoint.y;
        const double z = lasPoint.z;

        // Finding minimum and maximum coordinates
        // of points to create the bounding box.
        if (x < m_VoxelMesh.minX)
            m_VoxelMesh.minX = x;
        if (x > m_VoxelMesh.maxX)
            m_VoxelMesh.maxX = x;
        if (y < m_VoxelMesh.minY)
            m_VoxelMesh.minY = y;
        if (y > m_VoxelMesh.maxY)
            m_VoxelMesh.maxY = y;
        if (z < m_VoxelMesh.minZ)
            m_VoxelMesh.minZ = z;
        if (z > m_VoxelMesh.maxZ)
            m_VoxelMesh.maxZ = z;
    }

    std::cout << m_VoxelMesh.minX << " " << m_VoxelMesh.minY << " " << m_VoxelMesh.minZ << std::endl;

    if (moveToCoordinateCenter) {
        // Calculation of the difference between minimum
        // and maximum point coordinates.
        m_VoxelMesh.deltaX = m_VoxelMesh.maxX - m_VoxelMesh.minX;
        m_VoxelMesh.deltaY = m_VoxelMesh.maxY - m_VoxelMesh.minY;
        m_VoxelMesh.deltaZ = m_VoxelMesh.maxZ - m_VoxelMesh.minZ;

        for (u128 i = 0; i < points.size(); i++) {
            m_Points[i].x = m_Points[i].x - m_VoxelMesh.minX;
            m_Points[i].y = m_Points[i].y - m_VoxelMesh.minY;
            m_Points[i].z = m_Points[i].z - m_VoxelMesh.minZ;
        }
        m_VoxelMesh.minX = m_VoxelMesh.minY = m_VoxelMesh.minZ = 0;
        m_VoxelMesh.maxX = m_VoxelMesh.deltaX;
        m_VoxelMesh.maxY = m_VoxelMesh.deltaY;
        m_VoxelMesh.maxZ = m_VoxelMesh.deltaZ;
    }
    else {
        m_VoxelMesh.minX = m_VoxelMesh.minY = m_VoxelMesh.minZ = 0;
        // Calculation of the difference between minimum
        // and maximum point coordinates.
        m_VoxelMesh.deltaX = m_VoxelMesh.maxX - 0;
        m_VoxelMesh.deltaY = m_VoxelMesh.maxY - 0;
        m_VoxelMesh.deltaZ = m_VoxelMesh.maxZ - 0;
    }


}

VoxelMesh& LocalSymmetry::getVoxelMesh() {
    return m_VoxelMesh;
}

VoxelMesh LocalSymmetry::getVoxelMesh() const {
    return m_VoxelMesh;
}

VoxelVector& LocalSymmetry::getVoxelVector() {
    return m_VoxelVector;
}

VoxelVector LocalSymmetry::getVoxelVector() const {
    return m_VoxelVector;
}

std::vector<Voxel>& LocalSymmetry::getInterestingVoxels() {
    return m_InterestingVoxels;
}

std::vector<Voxel>& LocalSymmetry::getMaterialVoxels() {
    return m_MaterialVoxels;
}




VoxelMesh LocalSymmetry::calculateVoxelMeshByVoxelSideLength(const uint userInputX, const uint userInputY, const uint userInputZ) {
    // Calculating X, Y, Z and total voxel count.
    const uint x = static_cast<uint>(floor(m_VoxelMesh.deltaX / userInputX)) + 1;
    const uint y = static_cast<uint>(floor(m_VoxelMesh.deltaY / userInputY)) + 1;
    const uint z = static_cast<uint>(floor(m_VoxelMesh.deltaZ / userInputZ)) + 1;
    const uint count = x * y * z;

    // Saving the optimal values to the voxel mesh.
    m_VoxelMesh.count = count;
    m_VoxelMesh.sideX = userInputX;
    m_VoxelMesh.sideY = userInputY;
    m_VoxelMesh.sideZ = userInputZ;
    m_VoxelMesh.xSize = x;
    m_VoxelMesh.ySize = y;
    m_VoxelMesh.zSize = z;

    return m_VoxelMesh;
}

void LocalSymmetry::calculateVoxelMeshByMaximalVoxelCount(const uint userInput) {
    uint bestVoxelCount = 0;
    uint bestVoxelSideSize = 0;
    uint bestVoxelX = 0;
    uint bestVoxelY = 0;
    uint bestVoxelZ = 0;

    // Fitting the voxel mesh by X.
    for (int i = 1; i < static_cast<int>(m_VoxelMesh.deltaX); i++) {
        m_VoxelMesh.sideX = static_cast<int>(ceil(m_VoxelMesh.deltaX / i));
        m_VoxelMesh.sideY = static_cast<int>(ceil(m_VoxelMesh.deltaY / i));
        m_VoxelMesh.sideZ = static_cast<int>(ceil(m_VoxelMesh.deltaZ / i));
        m_VoxelMesh.xSize = std::min(i, static_cast<int>(ceil(m_VoxelMesh.deltaX / m_VoxelMesh.sideX))) + 1;
        m_VoxelMesh.ySize = static_cast<int>(ceil(m_VoxelMesh.deltaY / m_VoxelMesh.sideY)) + 1;
        m_VoxelMesh.zSize = static_cast<int>(ceil(m_VoxelMesh.deltaZ / m_VoxelMesh.sideZ)) + 1;
        m_VoxelMesh.count = m_VoxelMesh.xSize * m_VoxelMesh.ySize * m_VoxelMesh.zSize;

        // If the voxel count is larger than the
        // maximum input value, the search is over.
        if (m_VoxelMesh.count > userInput) {
            break;
        }

        // If a better option is found, it is stored in the voxel mesh.
        if (m_VoxelMesh.count > bestVoxelCount) {
            bestVoxelCount = m_VoxelMesh.count;
            bestVoxelSideSize = m_VoxelMesh.sideX;
            bestVoxelX = m_VoxelMesh.xSize;
            bestVoxelY = m_VoxelMesh.ySize;
            bestVoxelZ = m_VoxelMesh.zSize;

            // Voxel edge size is minimum 1, therefore the search is over.
            if (m_VoxelMesh.sideX == 1) {
                break;
            }
        }
    }

    // Fitting the voxel mesh by Y.
    for (int i = 1; i < static_cast<int>(m_VoxelMesh.deltaY); i++) {
        m_VoxelMesh.sideX = static_cast<int>(ceil(m_VoxelMesh.deltaY / i));
        m_VoxelMesh.xSize = static_cast<int>(ceil(m_VoxelMesh.deltaX / m_VoxelMesh.sideX)) + 1;
        m_VoxelMesh.ySize = std::min(i, static_cast<int>(ceil(m_VoxelMesh.deltaY / m_VoxelMesh.sideY))) + 1;
        m_VoxelMesh.zSize = static_cast<int>(ceil(m_VoxelMesh.deltaZ / m_VoxelMesh.sideZ)) + 1;
        m_VoxelMesh.count = m_VoxelMesh.xSize * m_VoxelMesh.ySize * m_VoxelMesh.zSize;

        // If the voxel count is larger than the
        // maximum input value, the search is over.
        if (m_VoxelMesh.count > userInput) {
            break;
        }

        // If a better option is found, it is stored in the voxel mesh.
        if (m_VoxelMesh.count > bestVoxelCount) {
            bestVoxelCount = m_VoxelMesh.count;
            bestVoxelSideSize = m_VoxelMesh.sideX;
            bestVoxelX = m_VoxelMesh.xSize;
            bestVoxelY = m_VoxelMesh.ySize;
            bestVoxelZ = m_VoxelMesh.zSize;

            // Voxel edge size is minimum 1, therefore the search is over.
            if (m_VoxelMesh.sideX == 1) {
                break;
            }
        }
    }

    // Fitting the voxel mesh by Z.
    for (int i = 1; i < static_cast<int>(m_VoxelMesh.deltaZ); i++) {
        m_VoxelMesh.sideX = static_cast<int>(ceil(m_VoxelMesh.deltaZ / i));
        m_VoxelMesh.xSize = static_cast<int>(ceil(m_VoxelMesh.deltaX / m_VoxelMesh.sideX)) + 1;
        m_VoxelMesh.ySize = static_cast<int>(ceil(m_VoxelMesh.deltaY / m_VoxelMesh.sideY)) + 1;
        m_VoxelMesh.zSize = std::min(i, static_cast<int>(ceil(m_VoxelMesh.deltaZ / m_VoxelMesh.sideZ))) + 1;
        m_VoxelMesh.count = m_VoxelMesh.xSize * m_VoxelMesh.ySize * m_VoxelMesh.zSize;

        // Voxel edge size is minimum 1, therefore the search is over.
        if (m_VoxelMesh.count > userInput) {
            break;
        }

        // If a better option is found, it is stored in the voxel mesh.
        if (m_VoxelMesh.count > bestVoxelCount) {
            bestVoxelCount = m_VoxelMesh.count;
            bestVoxelSideSize = m_VoxelMesh.sideX;
            bestVoxelX = m_VoxelMesh.xSize;
            bestVoxelY = m_VoxelMesh.ySize;
            bestVoxelZ = m_VoxelMesh.zSize;

            // Voxel edge size is minimum 1, therefore the search is over.
            if (m_VoxelMesh.sideX == 1) {
                break;
            }
        }
    }

    // If no solution is found, the whole scene is represented as 1 voxel.
    // Not nice for a user but nice for the algorithm :)
    if (bestVoxelCount == 0) {
        m_VoxelMesh.count = 1;
        m_VoxelMesh.sideX = static_cast<int>(std::max({ m_VoxelMesh.deltaX, m_VoxelMesh.deltaY, m_VoxelMesh.deltaZ }));
        m_VoxelMesh.xSize = 1;
        m_VoxelMesh.ySize = 1;
        m_VoxelMesh.zSize = 1;
    }
    else {
        // Saving the optimal values to the voxel mesh.
        m_VoxelMesh.count = bestVoxelCount;
        m_VoxelMesh.sideX = bestVoxelSideSize;
        m_VoxelMesh.xSize = bestVoxelX;
        m_VoxelMesh.ySize = bestVoxelY;
        m_VoxelMesh.zSize = bestVoxelZ;
    }
}

VoxelVector LocalSymmetry::buildVoxelVector() {
    m_VoxelVector = VoxelVector(m_VoxelMesh.zSize, std::vector<std::vector<Voxel>>(m_VoxelMesh.ySize, std::vector<Voxel>(m_VoxelMesh.xSize, Voxel())));

    // Setting voxel coordinates.
    for (uint z = 0; z < m_VoxelMesh.zSize; z++) {
        for (uint y = 0; y < m_VoxelMesh.ySize; y++) {
            for (uint x = 0; x < m_VoxelMesh.xSize; x++) {
                // Calculation of voxel coordinates.
                const uint voxelX = static_cast<uint>(floor(x * m_VoxelMesh.sideX));  // X coordinate.
                const uint voxelY = static_cast<uint>(floor(y * m_VoxelMesh.sideY));  // Y coordinate.
                const uint voxelZ = static_cast<uint>(floor(z * m_VoxelMesh.sideZ));  // Z coordinate.

                // Setting voxel coordinates.
                m_VoxelVector[z][y][x].x = voxelX;
                m_VoxelVector[z][y][x].y = voxelY;
                m_VoxelVector[z][y][x].z = voxelZ;
            }
        }
    }

    return m_VoxelVector;
}

std::tuple<VoxelVector, std::vector<Voxel>, std::vector<Voxel>> LocalSymmetry::findInterestingVoxels(const uint minClusterSize) {
    m_InterestingVoxels.clear();
    m_MaterialVoxels.clear();
    VoxelFunctions::clearVoxelVector(m_VoxelVector);

    // Adding all interesting voxels.
    for (const Point& point : m_Points) {
        const uint x = static_cast<uint>(floor((point.x - m_VoxelMesh.minX) / m_VoxelMesh.sideX));
        const uint y = static_cast<uint>(floor((point.y - m_VoxelMesh.minY) / m_VoxelMesh.sideY));
        const uint z = static_cast<uint>(floor((point.z - m_VoxelMesh.minZ) / m_VoxelMesh.sideZ));

        // If the current point lies on the voxel edge, it is ignored.
        if (
            Tolerance::isInTolerance(x, (point.x - m_VoxelMesh.minX) / m_VoxelMesh.sideX, 0.0001) ||
            Tolerance::isInTolerance(y, (point.y - m_VoxelMesh.minY) / m_VoxelMesh.sideY, 0.0001)
        )
        {
            continue;
        }

        // Searching the voxel in the vector and setting
        // super-interesting (temporary) and interesting.
        if (!m_VoxelVector[z][y][x].interesting) {
            m_VoxelVector[z][y][x].interesting = true;
            m_VoxelVector[z][y][x].material = true;

            m_InterestingVoxels.push_back(Voxel(x * m_VoxelMesh.sideX, y * m_VoxelMesh.sideY, z * m_VoxelMesh.sideZ));
        }
    }

    // Removing too small clusters.
    SymmetryFunctions::removeSmallMaterialClusters(m_InterestingVoxels, m_VoxelVector, m_VoxelMesh, minClusterSize);

    m_MaterialVoxels = m_InterestingVoxels;

    // Removing interesting voxels in the middle,
    // so that the whole planes are not searched.
    VoxelVector tempVoxels(m_VoxelVector);
    for (uint z = 0; z < m_VoxelVector.size(); z++) {
        for (uint y = 1; y < m_VoxelVector[z].size() - 1; y++) {
            for (uint x = 1; x < m_VoxelVector[z][y].size() - 1; x++) {
                if (m_VoxelVector[z][y].size() >= 3 &&
                    m_VoxelVector[z].size() >= 3 &&
                    tempVoxels[z][y][x].material &&
                    tempVoxels[z][y][x - 1].material &&
                    tempVoxels[z][y + 1][x].material &&
                    tempVoxels[z][y + 1][x + 1].material &&
                    tempVoxels[z][y - 1][x].material &&
                    m_VoxelVector[z][y][x].interesting
                )
                {
                    m_VoxelVector[z][y][x].interesting = false;

                    // Moving the voxel from interesting to material voxels.
                    m_InterestingVoxels.erase(std::remove(m_InterestingVoxels.begin(), m_InterestingVoxels.end(), m_VoxelVector[z][y][x]), m_InterestingVoxels.end());
                }
            }
        }
    }

    //// Removing super-interesting voxels in the middle,
    //// so that the whole planes are not searched.
    //for (uint z = 1; z < voxels.size() - 1; z++) {
    //    for (uint y = 0; y < voxels[z].size(); y++) {
    //        for (uint x = 1; x < voxels[z][y].size() - 1; x++) {
    //            if (
    //                tempVoxels[z][y][x].material &&
    //                tempVoxels[z][y][x - 1].material &&
    //                tempVoxels[z + 1][y][x].material &&
    //                tempVoxels[z][y][x + 1].material &&
    //                tempVoxels[z - 1][y][x].material &&
    //                voxels[z][y][x].interesting
    //            )
    //            {
    //                voxels[z][y][x].interesting = false;

    //                // Moving the voxel from interesting to material voxels.
    //                interesting.erase(std::remove(interesting.begin(), interesting.end(), voxels[z][y][x]), interesting.end());
    //            }
    //        }
    //    }
    //}

    // Removing super-interesting voxels in the middle,
    // so that the whole planes are not searched.
    //for (uint z = 1; z < voxels.size() - 1; z++) {
    //    for (uint y = 1; y < voxels[z].size() - 1; y++) {
    //        for (uint x = 0; x < voxels[z][y].size(); x++) {
    //            if (
    //                tempVoxels[z][y][x].material &&
    //                tempVoxels[z][y - 1][x].material &&
    //                tempVoxels[z + 1][y][x].material &&
    //                tempVoxels[z][y + 1][x].material &&
    //                tempVoxels[z - 1][y][x].material &&
    //                voxels[z][y][x].interesting

    //            )
    //            {
    //                voxels[z][y][x].interesting = false;

    //                // Moving the voxel from interesting to material voxels.
    //                interesting.erase(std::remove(interesting.begin(), interesting.end(), voxels[z][y][x]), interesting.end());
    //            }
    //        }
    //    }
    //}

    // Removing super-interesting voxels in the middle,
    // so that the whole planes are not searched.
    //for (uint z = 1; z < voxels.size() - 1; z++) {
    //    for (uint y = 1; y < voxels[z].size() - 1; y++) {
    //        for (uint x = 1; x < voxels[z][y].size() - 1; x++) {
    //            if (
    //                (tempVoxels[z][y][x].material &&
    //                tempVoxels[z][y][x - 1].material &&
    //                tempVoxels[z][y][x + 1].material &&
    //                tempVoxels[z + 1][y + 1][x].material &&
    //                tempVoxels[z - 1][y - 1][x].material)
    //                ||
    //                (tempVoxels[z][y][x].material &&
    //                tempVoxels[z][y][x - 1].material &&
    //                tempVoxels[z][y][x + 1].material &&
    //                tempVoxels[z - 1][y + 1][x].material &&
    //                tempVoxels[z + 1][y - 1][x].material)
    //                ||
    //                (tempVoxels[z][y][x].material &&
    //                tempVoxels[z][y - 1][x].material &&
    //                tempVoxels[z - 1][y][x + 1].material &&
    //                tempVoxels[z][y - 1][x].material &&
    //                tempVoxels[z + 1][y][x - 1].material)
    //                ||
    //                (tempVoxels[z][y][x].material &&
    //                tempVoxels[z][y - 1][x].material &&
    //                tempVoxels[z + 1][y][x + 1].material &&
    //                tempVoxels[z][y - 1][x].material &&
    //                tempVoxels[z - 1][y][x - 1].material)
    //                ||
    //                (tempVoxels[z][y][x].material &&
    //                tempVoxels[z + 1][y][x].material &&
    //                tempVoxels[z - 1][y][x].material &&
    //                tempVoxels[z][y + 1][x - 1].material &&
    //                tempVoxels[z][y - 1][x + 1].material)
    //                ||
    //                (tempVoxels[z][y][x].material &&
    //                tempVoxels[z + 1][y][x].material &&
    //                tempVoxels[z - 1][y][x].material &&
    //                tempVoxels[z][y - 1][x - 1].material &&
    //                tempVoxels[z][y + 1][x + 1].material)
    //                &&
    //                voxels[z][y][x].interesting
    //            )
    //            {
    //                voxels[z][y][x].interesting = false;

    //                // Moving the voxel from interesting to material voxels.
    //                interesting.erase(std::remove(interesting.begin(), interesting.end(), voxels[z][y][x]), interesting.end());
    //            }
    //        }
    //    }
    //}

    return { m_VoxelVector, m_InterestingVoxels, m_MaterialVoxels };
}

std::vector<LineSegment> LocalSymmetry::calculateLineSegmentsBetweenPoints(const std::vector<Point>& points, const double tolerance, const double minDistance, const double minMaterialRatio, const std::vector<std::pair<double, double>>& forbiddenAngles) {
    std::vector<LineSegment> lineSegments;

    // Calculation of every pair of points.
    for (uint i = 0; i < points.size() - 1; i++) {
        for (uint j = i + 1; j < points.size(); j++) {
            // Line segments between points with a greater Z coordinate difference than tolerance are ignored.
            if (Tolerance::isInTolerance(points[i].z, points[j].z, tolerance)) {
                // If the line segment length is smaller than the minimum distance, the line segment is ignored.
                const LineSegment ls(points[i], points[j]);
                if (ls.length() < minDistance) {
                    continue;
                }

                // If the line segment angle is a forbidden one, the line segment is ignored.
                const double angle = 180 * ls.angle(LineSegment(Point(0, 0, 0), Point(1, 0, 0))) / PI;
                bool forbiddenAngle = false;
                for (const auto& forbiddenAngleRange : forbiddenAngles) {
                    if (angle > forbiddenAngleRange.first && angle < forbiddenAngleRange.second) {
                        forbiddenAngle = true;
                        break;
                    }
                }
                if (forbiddenAngle) {
                    continue;
                }

                // If the material voxel ratio is smaller than a minimum material voxel ratio, the line segment is ignored.
                const double percent = SymmetryFunctions::calculateMaterialSectionRatioInLineSegment(ls, m_VoxelVector, m_VoxelMesh);
                if (percent < minMaterialRatio) {
                    continue;
                }
                
                lineSegments.push_back(ls);
            }
        }
    }

    // Sorting all line segments by length.
    std::sort(lineSegments.begin(), lineSegments.end(), [](const LineSegment& ls1, const LineSegment& ls2) { return ls1.length() < ls2.length(); });

    return lineSegments;
}


std::vector<std::vector<LineSegment>> LocalSymmetry::splitLineSegmentVectorIntoParts(const std::vector<LineSegment>& lineSegments, const double tolerance) {
    // Creating a vector and pushing a first part.
    std::vector<std::vector<LineSegment>> splitLineSegments;
    splitLineSegments.push_back(std::vector<LineSegment>());

    // Iterating through all line segments.
    for (uint i = 0; i < lineSegments.size(); i++) {
        // If no elements in vector, a first line segment is added.
        if (splitLineSegments.back().size() == 0) {
            splitLineSegments.back().push_back(lineSegments[i]);
        }
        // If a line segment length is inside of the allowed tolerance,
        // it is added to the same part as the previous one.
        else if (
            Tolerance::isInTolerance(
                splitLineSegments.back()[0].length(),
                lineSegments[i].length(),
                100 * tolerance
            )
            ) {
            splitLineSegments.back().push_back(lineSegments[i]);
        }
        // If a line segment length is outside of the allowed tolerance,
        // a new part is created, line segment is added to the new part.
        else {
            // If the last part contains only one line
            // segment, there will be no symmetry.
            if (splitLineSegments.back().size() == 1) {
                splitLineSegments.pop_back();
            }

            splitLineSegments.push_back(std::vector<LineSegment>());
            splitLineSegments.back().push_back(lineSegments[i]);
        }
    }

    // If the last part contains only one line
    // segment, there will be no symmetry.
    if (splitLineSegments.back().size() == 1) {
        splitLineSegments.pop_back();
    }

    return splitLineSegments;
}

std::vector<Point> LocalSymmetry::pointsFromMaterialVoxels() const {
    std::vector<Point> points;

    // Creating a new point from a voxel.
    for (uint z = 0; z < m_VoxelVector.size(); z++) {
        for (uint y = 0; y < m_VoxelVector[z].size(); y++) {
            for (uint x = 0; x < m_VoxelVector[z][y].size(); x++) {
                if (m_VoxelVector[z][y][x].material) {
                    points.push_back(
                        Point(
                            m_VoxelVector[z][y][x].x + m_VoxelMesh.sideX / 2,
                            m_VoxelVector[z][y][x].y + m_VoxelMesh.sideY / 2,
                            m_VoxelVector[z][y][x].z + m_VoxelMesh.sideZ / 2
                        )
                    );
                }
            }
        }
    }

    return points;
}

std::vector<Point> LocalSymmetry::pointsFromInterestingVoxels() const {
    std::vector<Point> points;

    // Creating a new point from a voxel.
    for (uint z = 0; z < m_VoxelVector.size(); z++) {
        for (uint y = 0; y < m_VoxelVector[z].size(); y++) {
            for (uint x = 0; x < m_VoxelVector[z][y].size(); x++) {
                if (m_VoxelVector[z][y][x].interesting) {
                    points.push_back(
                        Point(
                            m_VoxelVector[z][y][x].x + m_VoxelMesh.sideX / 2,
                            m_VoxelVector[z][y][x].y + m_VoxelMesh.sideY / 2,
                            m_VoxelVector[z][y][x].z + m_VoxelMesh.sideZ / 2
                        )
                    );
                }
            }
        }
    }

    return points;
}


VoxelVector SymmetryFunctions::getNormalizedVoxelVector(const VoxelVector& voxelVector, const VoxelMesh& voxelMesh) {
    // Vector for storing normalized voxels.
    VoxelVector normalisedVoxelVector(voxelMesh.zSize, std::vector<std::vector<Voxel>>(voxelMesh.ySize, std::vector<Voxel>(voxelMesh.xSize, Voxel())));

    for (uint z = 0; z < voxelMesh.zSize; z++) {
        for (uint y = 0; y < voxelMesh.ySize; y++) {
            for (uint x = 0; x < voxelMesh.xSize; x++) {
                normalisedVoxelVector[z][y][x] = Voxel(x, y, z, voxelVector[z][y][x].interesting, voxelVector[z][y][x].material);
            }
        }
    }

    return normalisedVoxelVector;
}

void SymmetryFunctions::findMaterialClusterItem(std::set<Voxel>& cluster, VoxelVector& voxelVector, const VoxelMesh& voxelMesh, std::stack<Voxel>& stack, const uint x, const uint y, const uint z) {
    // Checking whether the point is in voxel mesh bounds.
    if (x >= 0 && x < voxelMesh.xSize &&
        y >= 0 && y < voxelMesh.ySize &&
        z >= 0 && z < voxelMesh.zSize &&
        !voxelVector[z][y][x].checked && voxelVector[z][y][x].material
    )
    {
        voxelVector[z][y][x].checked = true;    // Setting the current voxel to checked. 
        stack.push(voxelVector[z][y][x]);       // Pushing a voxel to the stack.
        cluster.emplace(Voxel(x, y, z, true));  // Adding a voxel to the current cluster.
    }
}

std::vector<std::set<Voxel>> SymmetryFunctions::findClustersOfMaterialVoxels(VoxelVector& normalizedVoxels, const VoxelMesh& voxelMesh) {
    std::vector<std::set<Voxel>> clusters;

    // Searching for clusters in each voxel in symmetry (depth-first-search).
    for (uint z = 0; z < normalizedVoxels.size(); z++) {
        for (uint y = 0; y < normalizedVoxels[z].size(); y++) {
            for (uint x = 0; x < normalizedVoxels[z][y].size(); x++) {
                // If the current voxel has not been checked and is in symmetry, a new cluster has been found.
                if (!normalizedVoxels[z][y][x].checked && normalizedVoxels[z][y][x].material) {
                    std::set<Voxel> cluster;                                 // Creating a new cluster.
                    std::stack<Voxel> stack({ normalizedVoxels[z][y][x] });  // Creating a stack for flood fill.

                    // While stack is not empty, we push the new items on it.
                    while (!stack.empty()) {
                        // Getting the top item from the stack.
                        const Voxel top = stack.top();
                        stack.pop();

                        // Casting the coordinates of the top voxel.
                        const uint topX = static_cast<uint>(top.x / voxelMesh.sideX);
                        const uint topY = static_cast<uint>(top.y / voxelMesh.sideY);
                        const uint topZ = static_cast<uint>(top.z / voxelMesh.sideZ);

                        // 26 steps to for depth-first search of the neighborhood.
                        // 9 in the upper layer, 8 (all but current) in the current layer, 9 in the lower layer.
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY - 1, topZ + 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY - 1, topZ + 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY - 1, topZ + 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 0, topZ + 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 0, topZ + 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 0, topZ + 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 1, topZ + 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 1, topZ + 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 1, topZ + 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY - 1, topZ + 0);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY - 1, topZ + 0);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY - 1, topZ + 0);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 0, topZ + 0);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 0, topZ + 0);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 0, topZ + 0);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 1, topZ + 0);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 1, topZ + 0);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 1, topZ + 0);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY - 1, topZ - 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY - 1, topZ - 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY - 1, topZ - 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 0, topZ - 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 0, topZ - 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 0, topZ - 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 1, topZ - 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 1, topZ - 1);
                        findMaterialClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 1, topZ - 1);
                    }


                    clusters.push_back(cluster);  // Adding a new cluster to the list.
                }
            }
        }
    }

    return clusters;
}

void SymmetryFunctions::removeSmallMaterialClusters(std::vector<Voxel>& materialVoxels, VoxelVector& voxels, const VoxelMesh& voxelMesh, const uint minimalClusterSize) {
    // If minimum cluster size equals 1, no clusters have to be removed. Yaaayyy!
    if (minimalClusterSize == 1) {
        return;
    }

    // Searching for clusters of interesting voxels.
    const std::vector<std::set<Voxel>> clusters = SymmetryFunctions::findClustersOfMaterialVoxels(voxels, voxelMesh);

    // Removing clusters that contain too few voxels.
    for (uint i = 0; i < clusters.size(); i++) {
        if (clusters[i].size() < minimalClusterSize) {
            for (const Voxel& v : clusters[i]) {
                voxels[v.z][v.y][v.x].material = false;
                voxels[v.z][v.y][v.x].interesting = false;

                // Removing the obsolete interesting voxel from the vector.
                materialVoxels.erase(std::remove(materialVoxels.begin(), materialVoxels.end(), voxels[v.z][v.y][v.x]), materialVoxels.end());
            }
        }
    }
}

double SymmetryFunctions::calculateMaterialSectionRatioInLineSegment(const LineSegment& lineSegment, const VoxelVector& voxels, const VoxelMesh& voxelMesh) {
    // Calculating the number of segments.
    uint interestingSegments = 0;
    const uint segments = static_cast<uint>(10 * lineSegment.length());

    // Calculating the vector from P1 to P2.
    const Point p1 = lineSegment.p1;
    const Point p2 = lineSegment.p2;
    const Vector3d difference = p2 - p1;

    // Marching on the line segment and calculating
    // interesting line segment sections.
    for (uint i = 0; i < segments; i++) {
        double factor = static_cast<double>(i) / segments;
        const Vector3d p = Vector3d(
            p1.x + factor * difference.x,
            p1.y + factor * difference.y,
            p1.z + factor * difference.z
        );

        // Calculating a voxel position.
        const uint x = static_cast<uint>(floor((p.x - voxelMesh.minX) / voxelMesh.sideX));
        const uint y = static_cast<uint>(floor((p.y - voxelMesh.minY) / voxelMesh.sideY));
        const uint z = static_cast<uint>(floor((p.z - voxelMesh.minZ) / voxelMesh.sideZ));

        // If a voxel is interesting, the number of
        // interesting sections is incremented.
        if (voxels[z][y][x].material) {
            interestingSegments++;
        }
    }

    return static_cast<double>(interestingSegments) / segments;
}