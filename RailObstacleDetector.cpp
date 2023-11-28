#include <algorithm>
#include <iostream>

#include "LinearRegression/LinearRegression.hpp"
#include "RailObstacleDetector.hpp"


void RailObstacleDetector::mergeAllSymmetries() {
    std::vector<ReflectionSymmetry> symmetries{ m_Symmetries[0] };

    for (uint i = 1; i < m_Symmetries.size(); i++) {
        symmetries[0] = symmetries[0] &= m_Symmetries[i];
    }

    m_Symmetries = symmetries;
}

std::vector<Voxel> RailObstacleDetector::basicRailDetection(ReflectionSymmetry& symmetry) {
    const std::vector<LineSegment> lineSegments = symmetry.getLineSegments();
    const VoxelMesh voxelMesh = symmetry.getVoxelMesh();
    const Plane plane = symmetry.getPlane();

    // Iterating through each symmetry line segment.
    for (LineSegment lineSegment : lineSegments) {
        if (lineSegment.p1.x > lineSegment.p2.x) {
            auto temp = lineSegment.p1.x;
            lineSegment.p1.x = lineSegment.p2.x;
            lineSegment.p2.x = temp;
        }

        // Calculation of line segment vector and other parameters.
        const Vector3d directionalVector = lineSegment.p2 - lineSegment.p1;
        Point start = lineSegment.p1;
        Point end = lineSegment.p2;
        const uint numberOfSteps = 1000;

        // Detection of rail voxels.
        for (uint i = 0; i < numberOfSteps; i++) {
            const Point currentPoint = start + directionalVector * (static_cast<double>(i) / numberOfSteps);
            const Voxel currentVoxel(static_cast<uint>(currentPoint.x / voxelMesh.sideX), static_cast<uint>(currentPoint.y / voxelMesh.sideY), static_cast<uint>(currentPoint.z / voxelMesh.sideZ), true);
            const Voxel currentVoxelNonNormalized(currentVoxel.x * voxelMesh.sideX, currentVoxel.y * voxelMesh.sideY, currentVoxel.z * voxelMesh.sideZ);
            const Point oppositePoint = plane.calculatePointAcrossPlane(currentVoxelNonNormalized.centerPoint(voxelMesh), voxelMesh);
            const Voxel oppositeVoxelNonNormalized = Voxel(static_cast<uint>(oppositePoint.x / voxelMesh.sideX) * voxelMesh.sideX, static_cast<uint>(oppositePoint.y / voxelMesh.sideY) * voxelMesh.sideY, static_cast<uint>(oppositePoint.z / voxelMesh.sideZ) * voxelMesh.sideZ);
            const Voxel oppositeVoxel = Voxel(static_cast<uint>(oppositePoint.x / voxelMesh.sideX), static_cast<uint>(oppositePoint.y / voxelMesh.sideY), static_cast<uint>(oppositePoint.z / voxelMesh.sideZ));

            const std::vector<Voxel> symmetryVoxels = symmetry.getSymmetryVoxels();
            if (std::find(symmetryVoxels.begin(), symmetryVoxels.end(), currentVoxel) == symmetryVoxels.end()) {
                symmetry.addSymmetryVoxel(currentVoxel);
            }
            if (oppositePoint.x >= 0.0 && oppositePoint.x < voxelMesh.maxX &&
                oppositePoint.y >= 0.0 && oppositePoint.y < voxelMesh.maxY &&
                oppositePoint.z >= 0.0 && oppositePoint.z < voxelMesh.maxZ &&
                std::find(symmetryVoxels.begin(), symmetryVoxels.end(), oppositeVoxel) == symmetryVoxels.end()
            )
            {
                symmetry.addSymmetryVoxel(oppositeVoxel);
            }
        }
    }
    const Plane& accurateSymmetryPlane = accuratePlane(symmetry.getSymmetryVoxels(), symmetry.getVoxelMesh());
    const Plane rail1 = accurateSymmetryPlane.parallelPlane(-72);
    const Plane rail2 = accurateSymmetryPlane.parallelPlane(72);
    symmetry.setPlane(accurateSymmetryPlane);
    symmetry.setSymmetryVoxels(std::vector<Voxel>(), std::vector<Voxel>(), std::vector<Voxel>());

    // Calculation of line segment vector and other parameters.
    const Point start(0.0, rail1.calculateYCoordinateAtX(0.0), lineSegments[0].p1.z);
    const Point end(voxelMesh.maxX, rail1.calculateYCoordinateAtX(voxelMesh.maxX), lineSegments[0].p1.z);
    const Vector3d directionalVector = end - start;
    const Point start2(0.0, rail2.calculateYCoordinateAtX(0.0), lineSegments[0].p1.z);
    const Point end2(voxelMesh.maxX, rail2.calculateYCoordinateAtX(voxelMesh.maxX), lineSegments[0].p1.z);
    const Vector3d directionalVector2 = end2 - start2;
    const uint numberOfSteps = 1000;

    // Detection of rail voxels.
    int startIndex = -1;
    int endIndex = -1;

    std::vector<Voxel> rail1Voxels;
    std::vector<Voxel> rail2Voxels;
    for (uint i = 0; i < numberOfSteps; i++) {
        const Point currentPoint = start + directionalVector * (static_cast<double>(i) / numberOfSteps);
        const Voxel currentVoxel(static_cast<uint>(currentPoint.x / voxelMesh.sideX), static_cast<uint>(currentPoint.y / voxelMesh.sideY), static_cast<uint>(currentPoint.z / voxelMesh.sideZ), true);
        const Voxel currentVoxelNonNormalized(currentVoxel.x * voxelMesh.sideX, currentVoxel.y * voxelMesh.sideY, currentVoxel.z * voxelMesh.sideZ);
        
        const Point oppositePoint = start2 + directionalVector2 * (static_cast<double>(i) / numberOfSteps);;
        const Voxel oppositeVoxelNonNormalized = Voxel(static_cast<uint>(oppositePoint.x / voxelMesh.sideX) * voxelMesh.sideX, static_cast<uint>(oppositePoint.y / voxelMesh.sideY) * voxelMesh.sideY, static_cast<uint>(oppositePoint.z / voxelMesh.sideZ) * voxelMesh.sideZ);
        const Voxel oppositeVoxel = Voxel(static_cast<uint>(oppositePoint.x / voxelMesh.sideX), static_cast<uint>(oppositePoint.y / voxelMesh.sideY), static_cast<uint>(oppositePoint.z / voxelMesh.sideZ));

        const uint currentVoxelFull = symmetry.numberOfPointsInVoxel(currentVoxelNonNormalized, voxelMesh, m_Points) > 0;
        const uint oppositeVoxelFull = symmetry.numberOfPointsInVoxel(oppositeVoxelNonNormalized, voxelMesh, m_Points) > 0;

        if (startIndex == -1 && (currentVoxelFull || oppositeVoxelFull)) {
            startIndex = rail1Voxels.size();
        }
        else if ((currentVoxelFull && !oppositeVoxelFull) || (!currentVoxelFull && oppositeVoxelFull) || (currentVoxelFull && oppositeVoxelFull)) {
            endIndex = rail1Voxels.size();
        }

        if (std::find(rail1Voxels.begin(), rail1Voxels.end(), currentVoxel) == rail1Voxels.end()) {
            rail1Voxels.push_back(currentVoxel);
        }
        if (std::find(rail2Voxels.begin(), rail2Voxels.end(), oppositeVoxel) == rail2Voxels.end()) {
            rail2Voxels.push_back(oppositeVoxel);
        }
    }

    std::vector<Voxel> railVoxels;
    railVoxels.insert(railVoxels.end(), rail1Voxels.begin() + startIndex, rail1Voxels.begin() + endIndex);
    railVoxels.insert(railVoxels.end(), rail2Voxels.begin() + startIndex, rail2Voxels.begin() + endIndex);

    symmetry.setSymmetryVoxels(railVoxels, {}, {});

    return symmetry.getSymmetryVoxels();
}

std::vector<Point> RailObstacleDetector::obtainPointsForAccuratePlaneDetection(std::vector<Voxel> railVoxels, const VoxelMesh& voxelMesh) const {
    std::vector<Point> pointsForMSE;

    // Transformation to real-life coordinates and sorting of voxels.
    std::transform(railVoxels.begin(), railVoxels.end(), railVoxels.begin(), [&voxelMesh](const Voxel& v) { return Voxel(v.x * voxelMesh.sideX, v.y * voxelMesh.sideY, v.z * voxelMesh.sideZ); });
    std::sort(railVoxels.begin(), railVoxels.end(), [](const Voxel& v1, const Voxel& v2) { return (v1.x < v2.x) || (!(v1.x > v2.x) && v1.y < v2.y) || (!(v1.x > v2.x) && !(v1.y > v2.y) && v1.z < v2.z); });

    // Splitting voxels into groups with same X coordinates.
    std::vector<std::vector<Voxel>> splitRailVoxels;
    double currentX = -1.0;
    for (const Voxel& railVoxel : railVoxels) {
        if (railVoxel.x > currentX) {
            splitRailVoxels.push_back(std::vector<Voxel>({ railVoxel }));
        }
        else {
            splitRailVoxels.back().push_back(railVoxel);
        }

        currentX = railVoxel.x;
    }

    // Calculation of middle voxel points.
    for (const std::vector<Voxel>& railVoxelGroup : splitRailVoxels) {
        for (uint i = 0; i < railVoxelGroup.size() - 1; i++) {
            for (uint j = 1; j < railVoxelGroup.size(); j++) {
                pointsForMSE.push_back((railVoxelGroup[i].centerPoint(voxelMesh) + railVoxelGroup[j].centerPoint(voxelMesh)) / 2);
            }
        }
    }

    return pointsForMSE;
}

Plane RailObstacleDetector::accuratePlane(const std::vector<Voxel>& railVoxels, const VoxelMesh& voxelMesh) const {
    // Retrieving points for Mean Square Error method.
    const std::vector<Point> pointsForMSE = obtainPointsForAccuratePlaneDetection(railVoxels, voxelMesh);

    // Linear regression for accurate plane detection.
    LinearRegression regression(pointsForMSE);
    const auto& [coefficient, constant] = regression.fit();

    return Plane(Point(0, constant, 0), Vector3d(1.0, coefficient, 0.0).normalize(), Vector3d(0, 0, 1));
}

std::tuple<std::vector<Voxel>, std::vector<Voxel>, std::vector<Voxel>> RailObstacleDetector::detectObstaclesOnRailway(ReflectionSymmetry& symmetry) {
    // Retrieving basic symmetry properties.
    const std::vector<LineSegment> lineSegments = symmetry.getLineSegments();
    const VoxelMesh voxelMesh = symmetry.getVoxelMesh();
    const std::vector<Voxel> symmetryVoxels = symmetry.getSymmetryVoxels();
    const Plane symmetryPlane = symmetry.getPlane();
    const VoxelVector voxelVector = symmetry.getVoxelVector();
    
    std::vector<Voxel> railVoxels(symmetry.getSymmetryVoxels());
    std::vector<Voxel> obstacleVoxels;
    std::vector<Voxel> voxelsAcrossObstacle;

    for (const Voxel& basicSymmetryVoxel : symmetryVoxels) {
        // Detection of the potential obstacle above the current voxel.
        bool obstacleFlag = false;
        for (uint j = 2; j <= 5; j++) {
            if (basicSymmetryVoxel.z + j >= voxelVector.size() || basicSymmetryVoxel.y >= voxelVector[0].size() || basicSymmetryVoxel.x >= voxelVector[0][0].size()) {
                break;
            }
            if (voxelVector[basicSymmetryVoxel.z + j][basicSymmetryVoxel.y][basicSymmetryVoxel.x].material) {
                obstacleFlag = true;
                break;
            }
        }

        // If obstacle is detected, the other rail is also checked for the obstacle.
        if (obstacleFlag) {
            const Voxel unnormalizedVoxel(basicSymmetryVoxel.x * voxelMesh.sideX, basicSymmetryVoxel.y * voxelMesh.sideY, basicSymmetryVoxel.z * voxelMesh.sideZ);
            const Point oppositePoint = symmetryPlane.calculatePointAcrossPlane(unnormalizedVoxel.centerPoint(voxelMesh), voxelMesh);
            const Voxel oppositeVoxelNonNormalized = Voxel(static_cast<uint>(oppositePoint.x / voxelMesh.sideX) * voxelMesh.sideX, static_cast<uint>(oppositePoint.y / voxelMesh.sideY) * voxelMesh.sideY, static_cast<uint>(oppositePoint.z / voxelMesh.sideZ) * voxelMesh.sideZ);
            const Voxel oppositeVoxel = Voxel(static_cast<uint>(oppositePoint.x / voxelMesh.sideX), static_cast<uint>(oppositePoint.y / voxelMesh.sideY), static_cast<uint>(oppositePoint.z / voxelMesh.sideZ));

            // Detection of the potential obstacle above the opposite voxel.
            bool obstacleFlagAcross = false;
            for (uint j = 2; j <= 5; j++) {
                if (basicSymmetryVoxel.z + j >= voxelVector.size() || oppositeVoxel.y > voxelVector[0].size() || oppositeVoxel.x > voxelVector[0][0].size()) {
                    break;
                }
                // If the opposite voxel is not material, there is no rail, therefore, checking for obstacle is utter non-sense.
                if (!voxelVector[oppositeVoxel.z][oppositeVoxel.y][oppositeVoxel.x].material) {
                    break;
                }
                if (voxelVector[oppositeVoxel.z + j][oppositeVoxel.y][oppositeVoxel.x].material) {
                    obstacleFlagAcross = true;
                    break;
                }
            }

            // Removing the two voxels from rails.
            railVoxels.erase(std::remove(railVoxels.begin(), railVoxels.end(), basicSymmetryVoxel), railVoxels.end());
            railVoxels.erase(std::remove(railVoxels.begin(), railVoxels.end(), oppositeVoxel), railVoxels.end());
            
            // Appending the obstacle and opposite obstacle/non-obstacle voxels.
            obstacleVoxels.push_back(basicSymmetryVoxel);
            if (obstacleFlagAcross) {
                obstacleVoxels.push_back(oppositeVoxel);
            }
            else {
                voxelsAcrossObstacle.push_back(oppositeVoxel);
            }
        }
    }

    // Setting the voxels in the symmetry.
    symmetry.setSymmetryVoxels(railVoxels, obstacleVoxels, voxelsAcrossObstacle);

    return { railVoxels, obstacleVoxels, voxelsAcrossObstacle };
}

void RailObstacleDetector::addOriginalSceneMaterialVoxels(ReflectionSymmetry& symmetry) {
    const VoxelMesh voxelMesh = symmetry.getVoxelMesh();

    // Detecting original material voxels.
    LocalSymmetry ls;
    ls.setPoints(m_Points);
    ls.calculateVoxelMeshByVoxelSideLength(voxelMesh.sideX, voxelMesh.sideY, voxelMesh.sideZ);
    ls.buildVoxelVector();
    auto [originalVoxelVector, originalInterestingVoxels, originalMaterialVoxels] = ls.findInterestingVoxels(1);

    // Setting new material and interesting voxel properties.
    VoxelVector newVoxelVector = symmetry.getVoxelVector();
    for (uint z = 0; z < newVoxelVector.size(); z++) {
        for (uint y = 0; y < newVoxelVector[0].size(); y++) {
            for (uint x = 0; x < newVoxelVector[0][0].size(); x++) {
                newVoxelVector[z][y][x].interesting = originalVoxelVector[z][y][x].interesting;
                newVoxelVector[z][y][x].material = originalVoxelVector[z][y][x].material;
            }
        }
    }

    symmetry.setVoxelVector(newVoxelVector);
}



RailObstacleDetector::RailObstacleDetector(const std::vector<ReflectionSymmetry>& symmetries, const std::vector<Point>& points, const std::vector<Point>& filteredPoints) : m_Symmetries(symmetries), m_Points(points), m_FilteredPoints(filteredPoints)
{}

std::vector<ReflectionSymmetry> RailObstacleDetector::detectRails(const bool mergeAll) {
    std::vector<ReflectionSymmetry> symmetries;

    if (mergeAll) {
        mergeAllSymmetries();
    }

    // Iterating through each symmetry and detecting rails.
    for (ReflectionSymmetry& symmetry : m_Symmetries) {
        addOriginalSceneMaterialVoxels(symmetry);  // Adding the original scene material voxels back to symmetry.
        basicRailDetection(symmetry);              // Basic detection of rails in the symmetry.
        detectObstaclesOnRailway(symmetry);        // Detection of potential obstacles on the railway.
  
        symmetries.push_back(symmetry);
    }

    // Sorting railways according to number of voxels in the decreasing order.
    std::sort(symmetries.begin(), symmetries.end(), [](const ReflectionSymmetry& s1, const ReflectionSymmetry& s2) { return s1.getSymmetryVoxels().size() > s2.getSymmetryVoxels().size(); });

    return symmetries;
}