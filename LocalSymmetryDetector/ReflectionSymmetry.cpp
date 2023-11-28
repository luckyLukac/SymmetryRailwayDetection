#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>

#include "Structs/Tolerance.hpp"
#include "ReflectionSymmetry.hpp"


using namespace Symmetry;



std::pair<ReflectionSymmetry, bool> ReflectionSymmetry::detectTrivialReflectionSymmetry(const LineSegment& ls1, const LineSegment& ls2, const double tolerance, const double angleTolerance, const std::vector<std::pair<double, double>>& forbiddenDistances) {
    // At first, we have to connect the two line segments.
    LineSegment conn1(ls1.p1, ls2.p1);
    LineSegment conn2(ls1.p2, ls2.p2);

    // As we only search symmetries that lie in the same layer, a
    // false pair is returned if the two Z coordinates are not equal.
    if (!Tolerance::isInTolerance(ls1.p1.z, ls2.p1.z, tolerance)) {
        return { ReflectionSymmetry(), false };
    }

    // If the connections intersect with each other, it is necessary
    // to reconnect the both line segments with other two points.
    if (LineSegmentFunctions::doLineSegmentsIntersect(conn1, conn2, tolerance, true)) {
        conn1.p1 = ls1.p1;
        conn1.p2 = ls2.p2;
        conn2.p1 = ls1.p2;
        conn2.p2 = ls2.p1;
    }

    // Checking for forbidden distance ranges.
    const double dist1 = PointFunctions::distance(ls1.p1, ls2.p1);
    const double dist2 = PointFunctions::distance(ls1.p2, ls2.p2);
    for (const auto& forbiddenDistanceRange : forbiddenDistances) {
        if ((dist1 > forbiddenDistanceRange.first && dist1 < forbiddenDistanceRange.second) || (dist2 > forbiddenDistanceRange.first && dist2 < forbiddenDistanceRange.second)) {
            return { ReflectionSymmetry(), false };
        }
    }

    // Calculating center points of the both line segments.
    const Point center1 = conn1.midpoint();
    const Point center2 = conn2.midpoint();

    // Calculating vectors.
    const Vector3d axis = (center1 - center2).normalize();
    const Vector3d V1 = (conn1.p1 - conn1.p2).normalize();
    const Vector3d V2 = (conn2.p1 - conn2.p2).normalize();

    // If axis cannot be calculated, there is no symmetry.
    if (std::isnan(axis.x) || std::isnan(axis.y) || std::isnan(axis.z)) {
        return { ReflectionSymmetry(), false };
    }

    // Calculating the angles between the center line segment and the connections.
    const double angle1 = std::acos(axis.dot(V1));
    const double angle2 = std::acos(axis.dot(V2));

    // If both angles are equal (with tolerance) to PI/2, a symmetry is found.
    if ((std::isnan(angle1) || Tolerance::isInTolerance(angle1, PI / 2, angleTolerance * PI / 180)) && (std::isnan(angle2) || Tolerance::isInTolerance(angle2, PI / 2, angleTolerance * PI / 180))) {
        // Creating a plane and a vector of line segments.
        Vector3d point(center1.x, center1.y, center1.z);
        if (point.x == 0 && point.y == 0 && point.z == 0) {
            point.x = center2.x;
            point.y = center2.z;
            point.z = center2.z;
        }

        ReflectionSymmetry newSymmetry(Plane(point, axis, Vector3d(0, 0, 1)), { ls1, ls2 });
        newSymmetry.m_Points = m_Points;
        newSymmetry.m_VoxelMesh = m_VoxelMesh;
        newSymmetry.m_VoxelVector = m_VoxelVector;
        newSymmetry.m_InterestingVoxels = m_InterestingVoxels;
        newSymmetry.m_MaterialVoxels = m_MaterialVoxels;

        return { newSymmetry, true };
    }

    return { ReflectionSymmetry(), false };
}

std::vector<ReflectionSymmetry> ReflectionSymmetry::detectTrivialSymmetries(const std::vector<LineSegment>& lineSegments, const double tolerance, const double angleTolerance, const std::vector<std::pair<double, double>>& forbiddenDistances) {
    std::vector<ReflectionSymmetry> symmetries;

    // Splitting line segments into parts.
    const SplitLineSegments splitLineSegments = splitLineSegmentVectorIntoParts(lineSegments, tolerance);

    // Iterating through all line segment vector parts.
    for (const std::vector<LineSegment>& part : splitLineSegments) {
        if (part.empty()) {
            continue;
        }
        for (uint i = 0; i < part.size() - 1; i++) {
            for (uint j = i + 1; j < part.size(); j++) {
                const auto [symmetry, found] = detectTrivialReflectionSymmetry(part[i], part[j], tolerance, angleTolerance, forbiddenDistances);
                if (found) {
                    symmetries.push_back(symmetry);
                }
            }
        }
    }

    return symmetries;
}

std::vector<ReflectionSymmetry> ReflectionSymmetry::mergeTrivialSymmetries(std::vector<ReflectionSymmetry>& reflectionSymmetries, const double tolerance, const double angleTolerance) {
    std::vector<ReflectionSymmetry> mergedSymmetries;

    // Merging the found symmetries by Z axis and center point.
    for (ReflectionSymmetry& symmetry : reflectionSymmetries) {
        // Searching for a potential symmetry in the list, where X and
        // Y coordinates are the same as in the current symmetry, the
        // symmetry level is also the same.
        auto itr = std::find_if(
            mergedSymmetries.begin(),
            mergedSymmetries.end(),
            [&symmetry, &angleTolerance](ReflectionSymmetry& curSymmetry) {
                const Vector3d v1 = symmetry.getPlane().normalVector();
                const Vector3d v2 = curSymmetry.getPlane().normalVector();

                return (
                    Tolerance::isInTolerance(VectorFunctions::angle(v1, v2), 0, angleTolerance * PI / 180) ||
                    Tolerance::isInTolerance(VectorFunctions::angle(v1, v2), 2 * PI, angleTolerance * PI / 180)
                ) && Tolerance::isInTolerance(symmetry.getPlane().d, curSymmetry.getPlane().d, 100);
            }
        );

        // Adding the symmetry to the list if not exists.
        if (itr == mergedSymmetries.end()) {
            mergedSymmetries.push_back(symmetry);
        }
        // Updating the Symmetry object with new limits.
        else {
            const u128 i = std::distance(mergedSymmetries.begin(), itr);
            mergedSymmetries[i] = mergedSymmetries[i] &= symmetry;
        }
    }

    // Merging the found symmetries by Z axis and center point.
    std::vector<ReflectionSymmetry> mergedSymmetries2;
    for (ReflectionSymmetry& symmetry : mergedSymmetries) {
        // Searching for a potential symmetry in the list, where X and
        // Y coordinates are the same as in the current symmetry, the
        // symmetry level is also the same.
        auto itr = std::find_if(
            mergedSymmetries2.begin(),
            mergedSymmetries2.end(),
            [&symmetry, &angleTolerance](ReflectionSymmetry& curSymmetry) {
                const Vector3d v1 = symmetry.getPlane().normalVector();
                const Vector3d v2 = curSymmetry.getPlane().normalVector();
                const double angle = VectorFunctions::angle(v1, v2);

                return (
                    Tolerance::isInTolerance(VectorFunctions::angle(v1, v2), 0, angleTolerance * PI / 180) ||
                    Tolerance::isInTolerance(VectorFunctions::angle(v1, v2), 2 * PI, angleTolerance * PI / 180)
                    ) && Tolerance::isInTolerance(symmetry.getPlane().d, curSymmetry.getPlane().d, 100);
            }
        );

        // Adding the symmetry to the list if not exists.
        if (itr == mergedSymmetries2.end()) {
            mergedSymmetries2.push_back(symmetry);
        }
        // Updating the Symmetry object with new limits.
        else {
            const u128 i = std::distance(mergedSymmetries2.begin(), itr);
            mergedSymmetries2[i] = mergedSymmetries2[i] &= symmetry;
        }
    }

    return mergedSymmetries2;
}

void ReflectionSymmetry::addSymmetryVoxels() {
    m_SymmetryVoxels.clear();  // Clearing the previous possible artefacts in the list. Better safe than sorry.

    // Adding voxels according to reflection symmetry line segments.
    for (const LineSegment& ls : m_LineSegments) {
        // First line segment point.
        {
            // Getting the line segment point and calculating the voxel coordinates.
            const Point p = ls.p1;
            const uint x = static_cast<uint>(p.x / m_VoxelMesh.sideX) * m_VoxelMesh.sideX;
            const uint y = static_cast<uint>(p.y / m_VoxelMesh.sideY) * m_VoxelMesh.sideY;
            const uint z = static_cast<uint>(p.z / m_VoxelMesh.sideZ) * m_VoxelMesh.sideZ;

            // Creating a new voxel.
            Voxel v(x, y, z, true, true);
            v.inSymmetry = true;

            // If the voxel is not already present in the symmetry, it is added to the list.
            if (std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), v) == m_SymmetryVoxels.end()) {
                m_SymmetryVoxels.push_back(Voxel(x, y, z, true, true));
            }
        }
        // Second line segment point.
        {
            // Getting the line segment point and calculating the voxel coordinates.
            Point p = ls.p2;
            const uint x = static_cast<uint>(p.x / m_VoxelMesh.sideX) * m_VoxelMesh.sideX;
            const uint y = static_cast<uint>(p.y / m_VoxelMesh.sideY) * m_VoxelMesh.sideY;
            const uint z = static_cast<uint>(p.z / m_VoxelMesh.sideZ) * m_VoxelMesh.sideZ;

            // Creating a new voxel.
            Voxel v(x, y, z, true, true);
            v.inSymmetry = true;

            // If the voxel is not already present in the symmetry, it is added to the list.
            if (std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), v) == m_SymmetryVoxels.end()) {
                m_SymmetryVoxels.push_back(Voxel(x, y, z, true, true));
            }
        }
    }
}

void ReflectionSymmetry::addMaterialVoxels(const double maxDistanceFromPlane) {
    // Adding interesting voxels for each reflection symmetry.
    std::vector<std::tuple<double, Point, Voxel>> distances;  // Distances from the plane to interesting voxels and plane projection points.

    // Checking voxels by each layer.
    for (uint z = 0; z < m_VoxelMesh.zSize; z++) {
        distances.clear();  // Clearing distances from the previous iteration.

        for (uint y = 0; y < m_VoxelMesh.ySize; y++) {
            for (uint x = 0; x < m_VoxelMesh.xSize; x++) {
                // Getting a voxel from the mesh according to indices.
                const Voxel& v = m_VoxelVector[z][y][x];
                const double distanceToPlane = m_Plane.distanceToPoint(v.centerPoint(m_VoxelMesh));
                if (distanceToPlane > maxDistanceFromPlane) {
                    continue;
                }

                // If the voxel is not interesting, we move to the next one as soon as possible.
                if (!v.material) {
                    continue;
                }

                // Checking whether the voxel is already present in the list.
                bool exists = std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), v) != m_SymmetryVoxels.end();

                // If the voxel is (super) interesting and is not a part of the symmetry,
                // it is checked if the voxel has its counterpart voxel on the other side
                // of the symmetry plane.
                const Point voxelPoint = m_VoxelVector[z][y][x].centerPoint(m_VoxelMesh);
                const Point opposite = m_Plane.calculatePointAcrossPlane(voxelPoint, m_VoxelMesh);

                // Calculation projection point voxel indices.
                const Voxel oppositeVoxel = PointFunctions::voxelPositionFromPoint(opposite, m_VoxelMesh);
                const auto [xOpposite, yOpposite, zOpposite] = oppositeVoxel.getNormalizedCoordinates(m_VoxelMesh);

                // Checking whether the voxel is already present in the list.
                bool oppositeExists = std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), oppositeVoxel) != m_SymmetryVoxels.end();

                if (VoxelFunctions::isVoxelInsideVoxelMesh(oppositeVoxel, m_VoxelMesh) && m_VoxelVector[zOpposite][yOpposite][xOpposite].material) {
                    if (!exists) {
                        m_SymmetryVoxels.push_back(v);
                    }
                    if (v != oppositeVoxel && !oppositeExists) {
                        m_SymmetryVoxels.push_back(oppositeVoxel);
                    }

                    const Point oppositeLeftUpper = opposite - Vector3d(-0.5 * m_VoxelMesh.sideX, -0.5 * m_VoxelMesh.sideY, 0);
                    const Voxel oppositeLeftUpperVoxel = PointFunctions::voxelPositionFromPoint(oppositeLeftUpper, m_VoxelMesh);
                    bool oppositeLeftUpperExists = std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), oppositeLeftUpperVoxel) != m_SymmetryVoxels.end();
                    const auto [xOppositeLeftUpper, yOppositeLeftUpper, zOppositeLeftUpper] = oppositeLeftUpperVoxel.getNormalizedCoordinates(m_VoxelMesh);
                    if (VoxelFunctions::isVoxelInsideVoxelMesh(oppositeLeftUpperVoxel, m_VoxelMesh) && m_VoxelVector[zOppositeLeftUpper][yOppositeLeftUpper][xOppositeLeftUpper].material && !oppositeLeftUpperExists) {
                        m_SymmetryVoxels.push_back(oppositeLeftUpperVoxel);
                    }

                    const Point oppositeRightUpper = opposite - Vector3d(0.5 * m_VoxelMesh.sideX, -0.5 * m_VoxelMesh.sideY, 0);
                    const Voxel oppositeRightUpperVoxel = PointFunctions::voxelPositionFromPoint(oppositeRightUpper, m_VoxelMesh);
                    bool oppositeRightUpperExists = std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), oppositeRightUpperVoxel) != m_SymmetryVoxels.end();
                    const auto [xOppositeRightUpper, yOppositeRightUpper, zOppositeRightUpper] = oppositeRightUpperVoxel.getNormalizedCoordinates(m_VoxelMesh);
                    if (VoxelFunctions::isVoxelInsideVoxelMesh(oppositeRightUpperVoxel, m_VoxelMesh) && m_VoxelVector[zOppositeRightUpper][yOppositeRightUpper][xOppositeRightUpper].material && !oppositeRightUpperExists) {
                        m_SymmetryVoxels.push_back(oppositeRightUpperVoxel);
                    }

                    const Point oppositeLeftLower = opposite - Vector3d(-0.5 * m_VoxelMesh.sideX, 0.5 * m_VoxelMesh.sideY, 0);
                    const Voxel oppositeLeftLowerVoxel = PointFunctions::voxelPositionFromPoint(oppositeLeftLower, m_VoxelMesh);
                    bool oppositeLeftLowerExists = std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), oppositeLeftLowerVoxel) != m_SymmetryVoxels.end();
                    const auto [xOppositeLeftLower, yOppositeLeftLower, zOppositeLeftLower] = oppositeLeftLowerVoxel.getNormalizedCoordinates(m_VoxelMesh);
                    if (VoxelFunctions::isVoxelInsideVoxelMesh(oppositeLeftLowerVoxel, m_VoxelMesh) && m_VoxelVector[zOppositeLeftLower][yOppositeLeftLower][xOppositeLeftLower].material && !oppositeLeftLowerExists) {
                        m_SymmetryVoxels.push_back(oppositeLeftLowerVoxel);
                    }

                    const Point oppositeRightLower = opposite - Vector3d(0.5 * m_VoxelMesh.sideX, 0.5 * m_VoxelMesh.sideY, 0);
                    const Voxel oppositeRightLowerVoxel = PointFunctions::voxelPositionFromPoint(oppositeRightLower, m_VoxelMesh);
                    bool oppositeRightLowerExists = std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), oppositeRightLowerVoxel) != m_SymmetryVoxels.end();
                    const auto [xOppositeRightLower, yOppositeRightLower, zOppositeRightLower] = oppositeRightLowerVoxel.getNormalizedCoordinates(m_VoxelMesh);
                    if (VoxelFunctions::isVoxelInsideVoxelMesh(oppositeRightLowerVoxel, m_VoxelMesh) && m_VoxelVector[zOppositeRightLower][yOppositeRightLower][xOppositeRightLower].material && !oppositeRightLowerExists) {
                        m_SymmetryVoxels.push_back(oppositeRightLowerVoxel);
                    }
                }
            }
        }

        // If there are no distances in the list, no pairs have to be checked. Yippee ki-yay!!!
        if (distances.empty()) {
            continue;
        }
    }
}

void ReflectionSymmetry::removeSmallClusters(const uint minimumClusterSize) {
    // If the minimum cluster size equals 1, no clusters have to be removed. Yaaayyy!
    if (minimumClusterSize == 1) {
        return;
    }

    // Finding clusters for the current symmetry.
    std::vector<std::set<Voxel>> clusters = findClusters();

    // Removing clusters that contain too few voxels.
    for (uint i = 0; i < clusters.size(); i++) {
        const u128 clusterSize = clusters[i].size();  // Reading the current cluster size (number of voxels).

        // If the current cluster size is smaller than the minimum cluster size, the cluster is removed from the symmetry.
        if (clusterSize < minimumClusterSize) {
            for (const Voxel& v : clusters[i]) {
                const auto [x, y, z] = v.getNormalizedCoordinates(m_VoxelMesh);
                std::vector<Voxel>::iterator position = std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), Voxel(x, y, z, true, true));
                if (position != m_SymmetryVoxels.end()) {
                    m_SymmetryVoxels.erase(position);
                }
            }
        }
    }
}

void ReflectionSymmetry::addInterestingClusterItem(std::set<Voxel>& cluster, VoxelVector& voxelVector, const VoxelMesh& voxelMesh, std::stack<Voxel>& stack, const uint x, const uint y, const uint z) {
    // Checking whether the point is in voxel mesh bounds.
    if (x >= 0 && x < voxelMesh.xSize &&
        y >= 0 && y < voxelMesh.ySize &&
        z >= 0 && z < voxelMesh.zSize &&
        !voxelVector[z][y][x].checked && voxelVector[z][y][x].inSymmetry
        ) {
        voxelVector[z][y][x].checked = true;    // Setting the current voxel to checked. 
        stack.push(voxelVector[z][y][x]);       // Pushing a voxel to the stack.
        cluster.emplace(Voxel(x, y, z, true));  // Adding a voxel to the current cluster.
    }
}

void ReflectionSymmetry::clusterRecursiveStep(std::set<Voxel>& cluster, VoxelVector& voxelVector, const uint x, const uint y, const uint z) const {
    // If the voxel coordinates lie outside of the voxel mesh or the
    // voxel has been already checked, the recursion unfolds.
    if (x < 0 || x >= m_VoxelMesh.xSize ||
        y < 0 || y >= m_VoxelMesh.ySize ||
        z < 0 || z >= m_VoxelMesh.zSize ||
        voxelVector[z][y][x].checked
    )
    {
        return;
    }

    // Setting the current voxel to checked.
    voxelVector[z][y][x].checked = true;

    // If the current voxel is not in symmetry, the recursion starts to unfold.
    if (!voxelVector[z][y][x].inSymmetry) {
        return;
    }

    // Inserting the current voxel to the cluster.
    cluster.insert(Voxel(x, y, z, true));

    // 26 recursive steps to for depth-first search of the neighborhood.
    // 9 in the upper layer, 8 (all but current) in the current layer, 9 in the lower layer.
    clusterRecursiveStep(cluster, voxelVector, x - 1, y - 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y - 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y - 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 0, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 0, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 0, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y - 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y - 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y - 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 0, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 0, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 0, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y - 1, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y - 1, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y - 1, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 0, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 0, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 0, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 1, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 1, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 1, z - 1);
}

std::vector<std::set<Voxel>> ReflectionSymmetry::findClusters() {
    std::vector<std::set<Voxel>> clusters;

    // Cleaning the checked property.
    for (uint z = 0; z < m_VoxelVector.size(); z++) {
        for (uint y = 0; y < m_VoxelVector[z].size(); y++) {
            for (uint x = 0; x < m_VoxelVector[z][y].size(); x++) {
                m_VoxelVector[z][y][x].checked = false;
            }
        }
    }

    // Getting the 3D voxel vector.
    for (const Voxel& v : m_SymmetryVoxels) {
        // Setting all voxels in symmetry to true.
        const auto [x, y, z] = v.getNormalizedCoordinates(m_VoxelMesh);
        m_VoxelVector[z][y][x].inSymmetry = true;
    }

    // Searching for clusters in each voxel in symmetry (depth-first-search).
    for (const Voxel& v : m_SymmetryVoxels) {
        // Getting the coordinates of the voxel.
        const auto [x, y, z] = v.getNormalizedCoordinates(m_VoxelMesh);

        // If the current voxel has not been checked and is in symmetry, a new cluster has been found.
        if (!m_VoxelVector[z][y][x].checked) {
            std::set<Voxel> cluster;         // Creating a new cluster.
            std::stack<Voxel> stack({ v });  // Creating a stack for flood fill.

            // While stack is not empty, we push the new items on it.
            while (!stack.empty()) {
                // Getting the top item from the stack.
                const Voxel top = stack.top();
                stack.pop();

                // Casting the coordinates of the top voxel.
                const auto [topX, topY, topZ] = top.getNormalizedCoordinates(m_VoxelMesh);

                // 26 steps to for depth-first search of the neighborhood.
                // 9 in the upper layer, 8 (all but current) in the current layer, 9 in the lower layer.
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX - 1, topY - 1, topZ + 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 0, topY - 1, topZ + 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 1, topY - 1, topZ + 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX - 1, topY + 0, topZ + 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 0, topY + 0, topZ + 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 1, topY + 0, topZ + 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX - 1, topY + 1, topZ + 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 0, topY + 1, topZ + 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 1, topY + 1, topZ + 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX - 1, topY - 1, topZ + 0);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 0, topY - 1, topZ + 0);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 1, topY - 1, topZ + 0);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX - 1, topY + 0, topZ + 0);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 0, topY + 0, topZ + 0);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 1, topY + 0, topZ + 0);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX - 1, topY + 1, topZ + 0);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 0, topY + 1, topZ + 0);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 1, topY + 1, topZ + 0);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX - 1, topY - 1, topZ - 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 0, topY - 1, topZ - 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 1, topY - 1, topZ - 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX - 1, topY + 0, topZ - 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 0, topY + 0, topZ - 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 1, topY + 0, topZ - 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX - 1, topY + 1, topZ - 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 0, topY + 1, topZ - 1);
                addInterestingClusterItem(cluster, m_VoxelVector, m_VoxelMesh, stack, topX + 1, topY + 1, topZ - 1);
            }

            // Adding a new cluster to the list.
            clusters.push_back(cluster);
        }
    }

    return clusters;
}

void ReflectionSymmetry::processPointsOnThePlaneAndVoxelEdge() {
    const std::vector<uint>& pointsOnPlaneAndVoxelEdge = m_Plane.getPointsIndicesOnPlaneAndVoxelEdge(m_Points, m_VoxelMesh);

    // Processing each point.
    for (const uint pointIndex : pointsOnPlaneAndVoxelEdge) {
        const Point& p = m_Points[pointIndex];  // Retrieving the point.

        // Getting voxel coordinates.
        uint voxelX = static_cast<uint>(p.x / m_VoxelMesh.sideX);
        uint voxelY = static_cast<uint>(p.y / m_VoxelMesh.sideY);
        uint voxelZ = static_cast<uint>(p.z / m_VoxelMesh.sideZ);
        if (voxelX == m_VoxelMesh.xSize) {
            voxelX--;
        }
        if (voxelY == m_VoxelMesh.ySize) {
            voxelY--;
        }
        if (voxelZ == m_VoxelMesh.zSize) {
            voxelX--;
        }

        if (Tolerance::isInTolerance(fmod(voxelX, 1.0), 0.0, 0.0001) && Tolerance::isInTolerance(fmod(voxelY, 1.0), 0.0, 0.0001) && m_Plane.a == 0) {
            if (!m_VoxelVector[voxelZ][voxelY][voxelX - 1].material && !m_VoxelVector[voxelZ][voxelY - 1][voxelX - 1].material) {
                m_SymmetryVoxels.push_back(Voxel(voxelX - 1, voxelY - 1, voxelZ));
                m_SymmetryVoxels.push_back(Voxel(voxelX - 1, voxelY, voxelZ));
            }

            if (!m_VoxelVector[voxelZ][voxelY][voxelX].material && !m_VoxelVector[voxelZ][voxelY - 1][voxelX].material) {
                m_SymmetryVoxels.push_back(Voxel(voxelX, voxelY - 1, voxelZ));
                m_SymmetryVoxels.push_back(Voxel(voxelX, voxelY, voxelZ));
            }
        }
        else if (Tolerance::isInTolerance(fmod(voxelX, 1.0), 0.0, 0.0001) && Tolerance::isInTolerance(fmod(voxelY, 1.0), 0.0, 0.0001) && m_Plane.a == 1) {
            if (!m_VoxelVector[voxelZ][voxelY - 1][voxelX - 1].material && !m_VoxelVector[voxelZ][voxelY - 1][voxelX].material) {
                m_SymmetryVoxels.push_back(Voxel(voxelX - 1, voxelY - 1, voxelZ));
                m_SymmetryVoxels.push_back(Voxel(voxelX, voxelY - 1, voxelZ));
            }

            if (!m_VoxelVector[voxelZ][voxelY][voxelX - 1].material && !m_VoxelVector[voxelZ][voxelY][voxelX].material) {
                m_SymmetryVoxels.push_back(Voxel(voxelX - 1, voxelY, voxelZ));
                m_SymmetryVoxels.push_back(Voxel(voxelX, voxelY, voxelZ));
            }
        }
        else if (Tolerance::isInTolerance(fmod(voxelX, 1.0), 0.0, 0.0001) && !m_VoxelVector[voxelZ][voxelY][voxelX].material && !m_VoxelVector[voxelZ][voxelY][voxelX - 1].material) {
            m_SymmetryVoxels.push_back(Voxel(voxelX - 1, voxelY, voxelZ));
            m_SymmetryVoxels.push_back(Voxel(voxelX, voxelY, voxelZ));
        }
        else if (Tolerance::isInTolerance(fmod(voxelY, 1.0), 0.0, 0.0001) && !m_VoxelVector[voxelZ][voxelY][voxelX].material && !m_VoxelVector[voxelZ][voxelY - 1][voxelX].material) {
            m_SymmetryVoxels.push_back(Voxel(voxelX - 1, voxelY, voxelZ));
            m_SymmetryVoxels.push_back(Voxel(voxelX, voxelY, voxelZ));
        }
    }
}

void ReflectionSymmetry::postprocessSymmetries(std::vector<ReflectionSymmetry>& symmetries) {
    // Removing empty symmetries.
    for (uint i = 0; i < symmetries.size(); i++) {
        if (symmetries[i].getVoxelCount() <= 0) {
            symmetries.erase(symmetries.begin() + i);
            i--;
        }
    }

    // Sorting reflection symmetries by voxel count in each symmetry.
    std::sort(
        symmetries.begin(),
        symmetries.end(),
        [](const ReflectionSymmetry& s1, const ReflectionSymmetry& s2) {
            return s1.getVoxelCount() > s2.getVoxelCount();
        }
    );

    // Setting the index for each symmetry.
    for (uint i = 0; i < symmetries.size(); i++) {
        symmetries[i].setIndex(i + 1);
    }
}



ReflectionSymmetry::ReflectionSymmetry() :
    m_Index(0),
    m_Plane(Plane()),
    m_LineSegments(std::vector<LineSegment>())
{}

ReflectionSymmetry::ReflectionSymmetry(const Plane& plane, const std::vector<LineSegment>& lineSegments) :
    m_Index(0),
    m_Plane(plane),
    m_LineSegments(lineSegments)
{}

ReflectionSymmetry::ReflectionSymmetry(const Plane& plane, const std::vector<Voxel>& voxels) :
    m_Index(0),
    m_Plane(plane),
    m_SymmetryVoxels(voxels)
{}

ReflectionSymmetry::ReflectionSymmetry(const ReflectionSymmetry& symmetry1, const ReflectionSymmetry& symmetry2) : m_Index(0) {
    const Plane averagePlane(symmetry1.getPlane(), symmetry2.getPlane());
    m_Plane = averagePlane;

    std::vector<LineSegment> lineSegments;
    const std::vector<LineSegment> ls1 = symmetry1.getLineSegments();
    const std::vector<LineSegment> ls2 = symmetry2.getLineSegments();
    lineSegments.insert(lineSegments.end(), ls1.begin(), ls1.end());
    lineSegments.insert(lineSegments.end(), ls2.begin(), ls2.end());
    m_LineSegments = lineSegments;

    std::vector<Voxel> symmetryVoxels;
    const std::vector<Voxel> v1 = symmetry1.getSymmetryVoxels();
    const std::vector<Voxel> v2 = symmetry2.getSymmetryVoxels();
    symmetryVoxels.insert(symmetryVoxels.end(), v1.begin(), v1.end());
    symmetryVoxels.insert(symmetryVoxels.end(), v2.begin(), v2.end());
    m_SymmetryVoxels = symmetryVoxels;

    m_VoxelVector = symmetry1.m_VoxelVector;

    m_VoxelMesh = symmetry1.m_VoxelMesh;
}



bool ReflectionSymmetry::operator == (const ReflectionSymmetry& symmetry) const {
    // First condition: same planes.
    if (m_Plane == symmetry.m_Plane) {
        // Second condition: same voxel count.
        if (m_SymmetryVoxels.size() == symmetry.m_SymmetryVoxels.size()) {
            // Third condition: same voxels.
            for (uint i = 0; i < m_SymmetryVoxels.size(); i++) {
                if (m_SymmetryVoxels[i] == symmetry.m_SymmetryVoxels[i]) {
                    continue;
                }
                
                return false;
            }

            return true;
        }
    }

    return false;
}

ReflectionSymmetry ReflectionSymmetry::operator &= (const ReflectionSymmetry& symmetry) const {
    ReflectionSymmetry sym(*this, symmetry);
    return sym;
}



uint ReflectionSymmetry::getIndex() const {
    return m_Index;
}

void ReflectionSymmetry::setIndex(const uint index) {
    m_Index = index;
}

Plane ReflectionSymmetry::getPlane() const {
    return m_Plane;
}

void ReflectionSymmetry::setPlane(const Plane& plane) {
    m_Plane = plane;
}

std::vector<LineSegment> ReflectionSymmetry::getLineSegments() const {
    return m_LineSegments;
}

u128 ReflectionSymmetry::getLineSegmentCount() const {
    return m_LineSegments.size();
}

void ReflectionSymmetry::addLineSegment(const LineSegment& lineSegment) {
    m_LineSegments.push_back(lineSegment);
}

std::vector<Voxel> ReflectionSymmetry::getSymmetryVoxels() const {
    return m_SymmetryVoxels;
}

std::vector<Voxel> ReflectionSymmetry::getObstacleVoxels() const {
    return m_AsymmetryVoxels;
}

std::vector<Voxel> ReflectionSymmetry::getVoxelsAcrossObstacle() const {
    return m_AsymmetryVoxelsAcross;
}

void ReflectionSymmetry::setSymmetryVoxels(const std::vector<Voxel>& voxels, const std::vector<Voxel>& aVoxels, const std::vector<Voxel>& aVoxelsAcross) {
    m_SymmetryVoxels = voxels;
    m_AsymmetryVoxels = aVoxels;
    m_AsymmetryVoxelsAcross = aVoxelsAcross;
}

u128 ReflectionSymmetry::getVoxelCount() const {
    return m_SymmetryVoxels.size();
}

void ReflectionSymmetry::setVoxelVector(const VoxelVector& voxelVector) {
    m_VoxelVector = voxelVector;
}

void ReflectionSymmetry::addSymmetryVoxel(const Voxel& voxel) {
    m_SymmetryVoxels.push_back(voxel);
}

bool ReflectionSymmetry::doesVoxelContainPoint(const Voxel& voxel, const VoxelMesh& vm, const std::vector<Point>& points) {
    for (const Point& p : points) {
        if (p.x >= voxel.x && p.x < voxel.x + vm.sideX && p.y >= voxel.y && p.y < voxel.y + vm.sideY && p.z >= voxel.z && p.z < voxel.z + vm.sideZ) {
            return true;
        }
    }

    return false;
}

uint ReflectionSymmetry::numberOfPointsInVoxel(const Voxel& voxel, const VoxelMesh& vm, const std::vector<Point>& points) {
    uint number = 0;
    
    for (const Point& p : points) {
        if (p.x >= voxel.x && p.x < voxel.x + vm.sideX && p.y >= voxel.y && p.y < voxel.y + vm.sideY && p.z >= voxel.z && p.z < voxel.z + vm.sideZ) {
            number++;
        }
    }

    return number;
}



std::vector<ReflectionSymmetry> ReflectionSymmetry::detectReflectionSymmetries(const uint voxelSideX, const uint voxelSideY, const uint voxelSideZ, const double tolerance, const double angleTolerance, const uint minimumClusterSize, const double minimumLineSegmentDistance, const bool postprocessing) {
    // Voxelization procedures.
    calculateVoxelMeshByVoxelSideLength(voxelSideX, voxelSideY, voxelSideZ);
    buildVoxelVector();
    findInterestingVoxels(minimumClusterSize);
    
    // Symmetry detection.
    std::vector<Point> voxelPoints = pointsFromInterestingVoxels();
    std::vector<LineSegment> lineSegments = calculateLineSegmentsBetweenPoints(voxelPoints, tolerance, minimumLineSegmentDistance, 0.6, {{5, 175}, {185, 355}});
    std::vector<ReflectionSymmetry> trivialSymmetries = detectTrivialSymmetries(lineSegments, tolerance, angleTolerance, {{0, 130}, {155, MAX}});
    std::vector<ReflectionSymmetry> mergedSymmetries = mergeTrivialSymmetries(trivialSymmetries, tolerance, angleTolerance);

    // Processing of each symmetry.
    if (postprocessing) {
        for (ReflectionSymmetry& symmetry : mergedSymmetries) {
            symmetry.addSymmetryVoxels();
            symmetry.addMaterialVoxels();
            symmetry.removeSmallClusters(minimumClusterSize);
        }

        // Postprocessing of symmetries.
        postprocessSymmetries(mergedSymmetries);
    }

    return mergedSymmetries;
}

std::vector<ReflectionSymmetry> ReflectionSymmetry::detectLocalSymmetries(std::vector<ReflectionSymmetry>& reflectionSymmetries, const uint minSymmetrySize) {
    std::vector<ReflectionSymmetry> partialSymmetries;

    // Iterating through all the reflection symmetries.
    for (ReflectionSymmetry& symmetry : reflectionSymmetries) {
        const std::vector<std::set<Voxel>> clusters = symmetry.findClusters();
        const Plane plane = symmetry.getPlane();

        // Iterating through all the clusters.
        for (const std::set<Voxel>& cluster : clusters) {
            // If the cluster size is equal or bigger than the minimal
            // symmetry size, a new symmetry is added to the list.
            if (cluster.size() >= minSymmetrySize) {
                std::vector<Voxel> voxels(cluster.begin(), cluster.end());

                for (const Voxel& voxel : voxels) {
                    // Calculation of voxel and projection points.
                    const double voxelX = (voxel.x + 0.5) * m_VoxelMesh.sideX + m_VoxelMesh.minX;
                    const double voxelY = (voxel.y + 0.5) * m_VoxelMesh.sideY + m_VoxelMesh.minY;
                    const double voxelZ = (voxel.z + 0.5) * m_VoxelMesh.sideZ + m_VoxelMesh.minZ;
                    const Point voxelPoint(voxelX, voxelY, voxelZ);
                    const Point projection = plane.calculateProjectionPoint(voxelPoint, m_VoxelMesh);

                    // Calculating the projection point voxel coordinates.
                    const uint xp = static_cast<uint>(floor((projection.x - m_VoxelMesh.minX) / m_VoxelMesh.sideX));  // X index.
                    const uint yp = static_cast<uint>(floor((projection.y - m_VoxelMesh.minY) / m_VoxelMesh.sideY));  // Y index.
                    const uint zp = static_cast<uint>(floor((projection.z - m_VoxelMesh.minZ) / m_VoxelMesh.sideZ));  // Z index.

                    // If the voxel and the projection point voxel are the same,
                    // a reflection symmetry is also a partial symmetry.
                    if (voxel.x == xp && voxel.y == yp && voxel.z && zp) {
                        ReflectionSymmetry sym(plane, voxels);
                        sym.setIndex(static_cast<uint>(partialSymmetries.size()));
                        partialSymmetries.push_back(sym);
                        break;
                    }
                }
            }
        }
    }

    std::sort(partialSymmetries.begin(), partialSymmetries.end(), [](const ReflectionSymmetry& s1, const ReflectionSymmetry& s2) { return s1.getVoxelCount() > s2.getVoxelCount(); });

    return partialSymmetries;
}

std::vector<ReflectionSymmetry> ReflectionSymmetry::calculateNearestSymmetriesByPlane(const std::vector<ReflectionSymmetry>& symmetries, const Plane& desiredPlane, const VoxelMesh&) {
    std::vector<ReflectionSymmetry> nearestSymmetries(symmetries.size());
    std::vector<std::pair<double, uint>> areas(symmetries.size());

    for (int i = 0; i < symmetries.size(); i++) {
        const Plane plane = symmetries[i].getPlane();
        
        double angle = 180 * VectorFunctions::angle(plane.parallelVector(), desiredPlane.parallelVector()) / PI;
        if (angle > 180) {
            angle = 180 - Difference::difference(90, angle);
        }
        else if (angle > 90) {
            angle = 90 - Difference::difference(90, angle);
        }
        double diff = Difference::difference(plane.d, desiredPlane.d);
        areas[i] = std::make_pair(angle + diff, i);
    }

    std::sort(
        areas.begin(),
        areas.end(),
        [](auto area1, auto area2) {
            auto [distance1, index1] = area1;
            auto [distance2, index2] = area2;

            return distance1 < distance2;
        }
    );

    //for (auto area : areas) {
    for (int i = 0; i < areas.size(); i++) {
        const auto [distance, index] = areas[i];

        nearestSymmetries[i] = symmetries[index];
    }

    return nearestSymmetries;
}

std::vector<PointPosition> ReflectionSymmetry::calculatePositionsFromPlane(const std::vector<Point>& points, const VoxelVector& voxelVector, const VoxelMesh& voxelMesh, const Plane& plane, const bool validPlane, const bool normalize) const {
    const u128 numberOfPoints = points.size();  // Getting the point vector size.
    const uint normalizeFactorX = normalize ? voxelMesh.sideX : 1;
    const uint normalizeFactorY = normalize ? voxelMesh.sideY : 1;
    const uint normalizeFactorZ = normalize ? voxelMesh.sideZ : 1;
    const Vector3d normalizeVector(normalizeFactorX, normalizeFactorY, normalizeFactorZ);

    // If there is no plane, all points are undefined.
    if (!validPlane) {
        return std::vector<PointPosition>(numberOfPoints, PointPosition::notInSymmetry);
    }

    // Creating a new vector with positions.
    std::vector<PointPosition> positions(numberOfPoints, PointPosition::notInSymmetry);

    for (u128 i = 0; i < numberOfPoints; i++) {
        // Getting each single point and its voxel coordinates.
        const Point& point = points[i];
        const Voxel voxel = PointFunctions::voxelPositionFromPoint(point, voxelMesh);
        const auto [x, y, z] = voxel.getNormalizedCoordinates(voxelMesh);

        if (x == 111 && y == 43 && z == 3) {
            int xxx = 4;
        }

        // If the voxel is not at least material, its position is set to undefined.
        if (x >= m_VoxelMesh.xSize || y >= m_VoxelMesh.ySize || z >= m_VoxelMesh.zSize || !voxelVector[z][y][x].material) {
            positions[i] = PointPosition::notInSymmetry;
            continue;
        }

        // If the current point lies outside of an interesting voxel, its position is set to undefined.
        if (std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), voxel.divideVoxelWithVector(normalizeVector)) == m_SymmetryVoxels.end()) {
            if (x == 0 ||
                y == 0 ||
                (
                    plane.a == 1 &&
                    (
                        !Tolerance::isInTolerance(fmod(point.y / static_cast<double>(voxelMesh.sideY), voxelMesh.sideY), 0.0, 0.0001) ||
                        std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), voxelVector[z][y - 1][x - 1].divideVoxelWithVector(normalizeVector)) == m_SymmetryVoxels.end() ||
                        std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), voxelVector[z][y - 1][x].divideVoxelWithVector(normalizeVector)) == m_SymmetryVoxels.end()
                        )
                    )
                ) {
                positions[i] = PointPosition::notInSymmetry;
                continue;
            }
            else if (
                x == 0 ||
                y == 0 ||
                (
                    plane.a == 0 &&
                    (
                        !Tolerance::isInTolerance(fmod(point.x / static_cast<double>(voxelMesh.sideX), voxelMesh.sideX), 0.0, 0.0001) ||
                        std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), voxelVector[z][y - 1][x - 1].divideVoxelWithVector(normalizeVector)) == m_SymmetryVoxels.end() ||
                        std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), voxelVector[z][y][x - 1].divideVoxelWithVector(normalizeVector)) == m_SymmetryVoxels.end())
                    )
                ) {
                positions[i] = PointPosition::notInSymmetry;
                continue;
            }
            else if (plane.a != 0 && plane.a != 1) {
                if (voxelVector[z][y][x].inAsymmetry) {
                    positions[i] = PointPosition::asymmetry;
                }
                else {
                    positions[i] = PointPosition::notInSymmetry;
                }
                continue;
            }
        }

        // If the point and the projection point have the same coordinates, the
        // plane goes through the point. Therefore, its position is set to center.
        const Point pointProjection = plane.calculateProjectionPoint(point, voxelMesh);
        if (point == pointProjection) {
            positions[i] = PointPosition::center;
            continue;
        }

        // Calculating the center coordinates of the voxel where the point is located.
        const Point voxelCenter = voxel.centerPoint(voxelMesh);

        // Calculating the projection point to the symmetry plane.
        const Point projection = plane.calculateProjectionPoint(voxelCenter, voxelMesh);
        const auto [voxelProjectionX, voxelProjectionY, voxelProjectionZ] = PointFunctions::voxelPositionFromPoint(projection, voxelMesh).getNormalizedCoordinates(voxelMesh);

        positions[i] = PointPosition::rotational;
        continue;

        // If the projection point lies within the same voxel as the point, the position is set to center.
        // Note: voxels that are only touched by the symmetry plane on one edge are NOT center voxels.
        if (x == voxelProjectionX && !Tolerance::isInTolerance(voxel.x, projection.x, 0.001) &&
            y == voxelProjectionY && !Tolerance::isInTolerance(voxel.y, projection.y, 0.001) &&
            z == voxelProjectionZ && !Tolerance::isInTolerance(voxel.z, projection.z, 0.001)
            ) {
            positions[i] = PointPosition::center;
            continue;
        }

        // If the point is on the left side of the plane, left
        // position is added, right otherwise.
        if (plane.isPointOnTheLeftSide(point, voxelMesh)) {
            positions[i] = PointPosition::left;
        }
        else {
            positions[i] = PointPosition::right;
        }
    }

    return positions;
}