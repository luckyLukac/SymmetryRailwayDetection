#include <algorithm>
#include <cmath>
#include <queue>
#include <stack>
#include <tuple>

#include "Structs/Constants.hpp"
#include "Structs/Tolerance.hpp"
#include "RotationalSymmetry.hpp"

using namespace Symmetry;



// CLASS METHODS
// Splitting the line segments.
SplitLineSegmentPairs RotationalSymmetry::splitLineSegmentsInPairs(const std::vector<LineSegment>& lineSegments, const VoxelMesh& voxelMesh, const double lengthTolerance, const double angleTolerance) {
    const unsigned int maximumLevelOfRotation = 15;
    SplitLineSegmentPairs lineSegmentPairs(maximumLevelOfRotation);

    // If there are no line segments, there is nothing to do here.
    if (lineSegments.empty()) {
        return lineSegmentPairs;
    }

    // Creating a 2D vector (each vector in it will contain the
    // tuples with the angle and the two line segments).
    const double toleranceRadians = angleTolerance * (PI / 180);  // Transforming the angle tolerance from degrees to radians.
    std::vector<std::vector<double>> lineSegmentLengths(maximumLevelOfRotation);

    // Splitting line segments into parts.
    //auto splitLineSegments = splitLineSegmentVectorIntoParts(lineSegments, lengthTolerance);

    //// Iterating through all line segment vector parts.
    //for (std::vector<LineSegment>& part : splitLineSegments) {
    //    // Checking every possible combination.
    //    for (unsigned long long i = 0; i < part.size() - 1; i++) {
    //        for (unsigned long long j = i + 1; j < part.size(); j++) {
    //            // If the two line segments do not lie on the same height,
    //            // we move to the next pair as quickly as possible.
    //            if (!Tolerance::isInTolerance(part[i].p1.z, part[j].p1.z, 0.0001)) {
    //                continue;
    //            }

    //            // Calculating the angle between the two line segments.
    //            // Smaller angle in the interval of [0°, 180°] is
    //            // chosen in order to prevent programmer's hair loss.
    //            const double angle1 = part[i].angle(part[j]);
    //            const double angle2 = part[j].angle(part[i]);
    //            const double angle = std::min(angle1, angle2);

    //            // Calculating the rotation level in radians. If the level is not an integer,
    //            // the angle is not interesting for a potential rotational symmetry.
    //            // If the rotation level is bigger than the maximum level, we move on as well.
    //            const double levelFloat = 2.0 / (1 - (angle / PI));
    //            const double innerAngle = (levelFloat - 2) * PI / levelFloat;  // Calculation of the inner angle of the polygon.
    //            const unsigned int level = static_cast<unsigned int>(std::round(levelFloat));

    //            if (static_cast<unsigned int>(level) >= maximumLevelOfRotation) {
    //                continue;
    //            }

    //            if (!Tolerance::isInTolerance(angle, innerAngle, toleranceRadians) || levelFloat > maximumLevelOfRotation) {
    //                continue;
    //            }

    //            // If the current level vector is empty, a new vector of line segment pairs is added.
    //            // A new line segment pair is added to the first vector of the current vector.
    //            // Line segment pairs are divided by separate vectors on the second level based by
    //            // their length (line segment pairs with the same length belong to the same vector).
    //            if (lineSegmentPairs[level].empty()){
    //                lineSegmentPairs[level].push_back(std::vector<std::vector<LineSegmentPair>>());
    //                for (unsigned int layer = 0; layer < static_cast<unsigned int>(voxelMesh.zSize); layer++) {
    //                    lineSegmentPairs[level][0].push_back(std::vector<LineSegmentPair>());
    //                }

    //                const unsigned int layer = static_cast<unsigned int>(part[i].p1.z / voxelMesh.sideZ);
    //                lineSegmentPairs[level][0][layer].push_back(LineSegmentPair(part[i], part[j]));
    //                lineSegmentLengths[level].push_back(part[i].length());

    //                continue;
    //            }

    //            // Searching for a vector of line segment pairs with the same length as the current line segment.
    //            int foundIndex = -1;
    //            for (unsigned long k = 0; k < lineSegmentPairs[level].size(); k++) {
    //                const double d1 = lineSegmentLengths[level][k];  // Current vector line segment length.
    //                const double d2 = part[i].length();           // The current line segment length.

    //                // If the two line segment lengths are within a tolerance,
    //                // we have found the right vector for our new pair. Yaaaay!!!
    //                if (Tolerance::isInTolerance(d1, d2, 0.01)) {
    //                    foundIndex = static_cast<int>(k);
    //                    break;
    //                }
    //            }

    //            // If there is no vector with the suitable line segment length,
    //            // a new vector with a new line segment length is added.
    //            if (foundIndex == -1) {
    //                lineSegmentPairs[level].push_back(std::vector<std::vector<LineSegmentPair>>());
    //                const unsigned int index = static_cast<unsigned int>(lineSegmentPairs[level].size()) - 1;

    //                for (uint layer = 0; layer < voxelMesh.zSize; layer++) {
    //                    lineSegmentPairs[level][index].push_back(std::vector<LineSegmentPair>());
    //                }

    //                const unsigned int layer = static_cast<unsigned int>(part[i].p1.z / voxelMesh.sideZ);
    //                lineSegmentPairs[level][index][layer].push_back(LineSegmentPair(part[i], part[j]));
    //                lineSegmentLengths[level].push_back(part[i].length());
    //            }
    //            // Otherwise, a line segment pair is added to the beforefound vector.
    //            else {
    //                const unsigned int layer = static_cast<unsigned int>(part[i].p1.z / voxelMesh.sideZ);
    //                lineSegmentPairs[level][static_cast<unsigned int>(foundIndex)][layer].push_back(LineSegmentPair(part[i], part[j]));
    //            }
    //        }
    //    }
    //}

    //// Removing groups that contain too few (less than the rotation level) number of line segment pairs.
    //for (unsigned long level = 2; level < maximumLevelOfRotation; level++) {
    //    for (unsigned long group = 0; group < lineSegmentPairs[level].size(); group++) {
    //        for (unsigned long layer = 0; layer < lineSegmentPairs[level][group].size(); layer++) {
    //            // If the layer contains too few items, the whole layer is deleted.
    //            if (2 * lineSegmentPairs[level][group][layer].size() < level) {
    //                lineSegmentPairs[level][group].erase(lineSegmentPairs[level][group].begin() + layer);
    //                layer--;
    //            }
    //            else {
    //                // Getting the number of unique line segments.
    //                //std::set<LineSegment> _lineSegments;
    //                //for (unsigned long lineSegment = 0; lineSegment < lineSegmentPairs[level][group][layer].size(); lineSegment++) {
    //                //    auto& ls1 = lineSegmentPairs[level][group][layer][lineSegment].ls1;  // First line segment.
    //                //    auto& ls2 = lineSegmentPairs[level][group][layer][lineSegment].ls2;  // Second line segment.

    //                //    // Inserting the both line segments.
    //                //    _lineSegments.insert(ls1);
    //                //    _lineSegments.insert(ls2);

    //                //    // If there is equal or larger number of line segments
    //                //    // than the current level, we're all good.
    //                //    if (_lineSegments.size() >= level) {
    //                //        break;
    //                //    }
    //                //}

    //                //// If there is less line segments than the current level of rotation, the layer is removed.
    //                //if (_lineSegments.size() < level) {
    //                //    lineSegmentPairs[level][group].erase(lineSegmentPairs[level][group].begin() + layer);
    //                //    layer--;
    //                //}
    //            }
    //        }

    //        // If there is an empty group of line segment pairs, it is deleted.
    //        if (lineSegmentPairs[level][group].empty()) {
    //            lineSegmentPairs[level].erase(lineSegmentPairs[level].begin() + group);
    //            group--;
    //        }
    //    }
    //}

    return lineSegmentPairs;
}

// Getting one recursive vector combination.
void RotationalSymmetry::buildCombinationsVectorRecursively(std::vector<std::vector<int> >& combinations, std::vector<int>& currentCombination, int N, int left, int K) {
    // Adding a new combination to the vector.
    if (K == 0) {
        combinations.push_back(currentCombination);
        return;
    }

    // Adding the new elements to the current combination.
    for (int i = left; i < N; i++) {
        currentCombination.push_back(i);
        buildCombinationsVectorRecursively(combinations, currentCombination, N, i + 1, K - 1);  // Recursive build of a vector.
        currentCombination.pop_back();
    }
}

// Getting all vector combinations with size K from a vector with the size of N.
std::vector<std::vector<int>> RotationalSymmetry::allVectorCombinations(const int N, const int K) {
    std::vector<std::vector<int>> combinations;
    std::vector<int> currentCombination;

    // Recursive build of the combination vector.
    buildCombinationsVectorRecursively(combinations, currentCombination, N, 0, K);

    return combinations;
}

// Splitting line segment pairs by length between pair center points.
std::vector<LineSegmentPairGroup> RotationalSymmetry::splitPairsIntoGroups(const SplitLineSegmentPairs& pairs) {
    std::vector<LineSegmentPairGroup> groups;

    // Iterating through each pair.
    for (unsigned int level = 3; level < pairs.size(); level++) {
        for (unsigned int distanceGroup = 0; distanceGroup < pairs[level].size(); distanceGroup++) {
            for (unsigned int layer = 0; layer < pairs[level][distanceGroup].size(); layer++) {
                for (const LineSegmentPair& pair : pairs[level][distanceGroup][layer]) {
                    // Distance between pair center points calculation.
                    const double distance = PointFunctions::distance(pair.ls1.midpoint(), pair.ls2.midpoint());

                    // Finding the group with the current parameters in the vector.
                    auto it = std::find_if(
                        groups.begin(),
                        groups.end(),
                        [&level, &distanceGroup, &layer, &distance](const LineSegmentPairGroup& group) {
                            return (
                                group.rotation == level &&
                                group.level == distanceGroup &&
                                group.layer == layer &&
                                Tolerance::isInTolerance(distance, group.distanceBetweenPairCenterPoints, 0.001)
                            );
                        }
                    );

                    // If a group has been found, a new line segment pair is added.
                    if (it != groups.end()) {
                        it->addPair(pair);
                    }
                    // Otherwise, a new group is added.
                    else {
                        std::vector<LineSegmentPair> currentPairVector({pair});
                        groups.push_back(LineSegmentPairGroup(currentPairVector, level, distanceGroup, layer));
                    }
                }
            }
        }
    }

    // Deleting the groups that have too little pairs.
    groups.erase(
        std::remove_if(
            groups.begin(),
            groups.end(),
            [](const LineSegmentPairGroup& group) {
                return 2 * group.pairs.size() < group.rotation;
            }
        ),
        groups.end()
    );

    return groups;
}

// Splitting line segment groups into subgroups.
std::vector<LineSegmentSubgroup> RotationalSymmetry::splitLineSegmentGroupsIntoSubgroups(const std::vector<LineSegmentPairGroup>& groups) {
    std::vector<LineSegmentSubgroup> subgroups;

    // Iterating through all groups.
    for (const LineSegmentPairGroup& group : groups) {
        // Building a graph and getting vertices and edges.
        const LineSegmentGraph graph(group.pairs, group.rotation);      // Building a line segment graph.
        std::vector<Vertex> vertices = graph.getVertices();             // Vertices represent line segments.
        const std::vector<std::vector<bool>> edges = graph.getEdges();  // Edges represent relations between line segment pairs.

        // Flood Fill® algorithm.
        // Iterating through each start vertex.
        for (unsigned int startVertex = 0; startVertex < vertices.size(); startVertex++) {
            // If a vertex at the current index has already been checked
            // in the previous flood fills, the graph component has
            // already been found in a previous iteration.
            if (vertices[startVertex].getChecked()) {
                continue;
            }

            // Creating the data structures for Flood Fill.
            std::stack<unsigned int> stack({ startVertex });                                             // Stack for storing indices of the line segments, found by Flood Fill.
            std::vector<LineSegment> currentGraphComponent({ vertices[startVertex].getLineSegment() });  // Data structure for storing current graph component line segments.
            vertices[startVertex].setChecked();                                                          // Setting the start vertex to checked.

            // While there are line segments on the stack, the loop is continued.
            while (!stack.empty()) {
                // Getting the top element of the stack and its removal
                const unsigned int element = stack.top();
                stack.pop();

                // Getting the connected vertices of the current vertex.
                for (unsigned int foundVertex = element + 1; foundVertex < vertices.size(); foundVertex++) {
                    // If a vertex has not been checked previously and there is
                    // a connection between the current and the found element,
                    // a new connected line segment has been found.
                    if (!vertices[foundVertex].getChecked() && edges[element][foundVertex]) {
                        stack.push(foundVertex);                                                  // Adding a new line segment index to the stack.
                        currentGraphComponent.push_back(vertices[foundVertex].getLineSegment());  // Adding a new line segment to the current component.
                        vertices[foundVertex].setChecked();                                       // Setting the current vertex to checked.
                    }
                }
            }

            // Adding a new line segment subgroup with the current graph component line segments.
            subgroups.push_back(LineSegmentSubgroup(currentGraphComponent, group.rotation));
        }
    }

    return subgroups;
}

// Building line segment midpoints circumcircles.
std::vector<Circumcircle> RotationalSymmetry::buildLineSegmentMidpointsCircumcircles(const SplitLineSegmentPairs& pairs) {
    std::vector<Circumcircle> circumcircles;

    // Splitting line segment pairs into separate groups by length between pair center points.
    std::vector<LineSegmentPairGroup> groups = splitPairsIntoGroups(pairs);
    std::vector<LineSegmentSubgroup> subgroups = splitLineSegmentGroupsIntoSubgroups(groups);

    // Finding maximum number of line segments in a class.
    unsigned int maxNumberOfLineSegments = 0;
    for (const LineSegmentPairGroup& group : groups) {
        if (group.pairs.size() > maxNumberOfLineSegments) {
            maxNumberOfLineSegments = static_cast<unsigned int>(group.pairs.size());
        }
    }

    std::vector<std::vector<LineSegment>> bisections(2 * maxNumberOfLineSegments, std::vector<LineSegment>(2 * maxNumberOfLineSegments));

    for (const LineSegmentSubgroup& subgroup : subgroups) {
        std::vector<std::vector<int>> combinations = allVectorCombinations(static_cast<int>(subgroup.lineSegments.size()), 3);  // Getting all vector index triplet combinations.

        // Calculating all bisections.
        for (unsigned int i = 0; i < subgroup.lineSegments.size() - 1; i++) {
            for (unsigned int j = i + 1; j < subgroup.lineSegments.size(); j++) {
                // Retrieving the line segments at i and j indices.
                const LineSegment ls1 = subgroup.lineSegments[i];
                const LineSegment ls2 = subgroup.lineSegments[j];

                // Calculation of midpoints.
                const Point midPoint1 = ls1.midpoint();
                const Point midPoint2 = ls2.midpoint();

                const LineSegment connection(midPoint1, midPoint2);

                // Calculation of the bisection,for [i][j].
                LineSegment bisection = connection.perpendicularLineSegment();

                bisections[i][j].p1 = bisection.p1;
                bisections[i][j].p2 = bisection.p2;
            }
        }

        // Iterating through all combinations and calculating circumcircles.
        for (const std::vector<int>& combination : combinations) {
            // Retrieving the line segments at combination indices.
            const LineSegment ls1 = subgroup.lineSegments[static_cast<unsigned int>(combination[0])];
            const LineSegment ls2 = subgroup.lineSegments[static_cast<unsigned int>(combination[1])];
            const LineSegment ls3 = subgroup.lineSegments[static_cast<unsigned int>(combination[2])];

            // Calculation of the two bisections, the center point and the radius.
            const LineSegment bisection1 = bisections[static_cast<unsigned int>(combination[0])][static_cast<unsigned int>(combination[1])];
            const LineSegment bisection2 = bisections[static_cast<unsigned int>(combination[0])][static_cast<unsigned int>(combination[2])];

            // Calculation of the intersection between the both bisections.
            const auto [center, validCenter] = LineSegmentFunctions::calculateIntersection(bisection1, bisection2);
            if (!validCenter) {
                continue;
            }
            const double radius = PointFunctions::distance(center, ls1.midpoint());

            // Finding the currrent circumcircle in the vector.
            auto it = std::find_if(
                circumcircles.begin(),
                circumcircles.end(),
                [&center, &radius, &subgroup](const Circumcircle& circumcircle) {
                    return (
                        circumcircle.rotation == subgroup.rotation &&
                        center == circumcircle.center &&
                        Tolerance::isInTolerance(radius, circumcircle.radius, 0.0001)
                    );
                }
            );

            // If a circumcircle is already present in the
            // vector, we only have to add the line segments.
            if (it != circumcircles.end()) {
                // Adding the first line segment.
                //if (!it->lineSegmentsInCircumcircle[static_cast<unsigned int>(combination[0])]) {
                    //it->lineSegmentsInCircumcircle[static_cast<unsigned int>(combination[0])] = true;
                    it->addLineSegment(ls1);
                //}

                // Adding the second line segment.
                //if (!it->lineSegmentsInCircumcircle[static_cast<unsigned int>(combination[1])]) {
                    //it->lineSegmentsInCircumcircle[static_cast<unsigned int>(combination[1])] = true;
                    it->addLineSegment(ls2);
                //}

                // Adding the third line segment.
                //if (!it->lineSegmentsInCircumcircle[static_cast<unsigned int>(combination[2])]) {
                    //it->lineSegmentsInCircumcircle[static_cast<unsigned int>(combination[2])] = true;
                    it->addLineSegment(ls3);
                //}
            }
            // Otherwise, a new circumcircle is added to the vector.
            else {
                std::vector<bool> lineSegmentsInCircumcircle(subgroup.lineSegments.size(), false);
                lineSegmentsInCircumcircle[static_cast<unsigned int>(combination[0])] = true;
                lineSegmentsInCircumcircle[static_cast<unsigned int>(combination[1])] = true;
                lineSegmentsInCircumcircle[static_cast<unsigned int>(combination[2])] = true;

                circumcircles.push_back(Circumcircle(subgroup.rotation, center, radius, lineSegmentsInCircumcircle, std::vector<LineSegment>({ ls1, ls2, ls3 })));
            }
        }
    }

    // Removing circumcircles that contain less line segments than the rotation level.
    circumcircles.erase(
        std::remove_if(
            circumcircles.begin(),
            circumcircles.end(),
            [](const Circumcircle& circumcircle) {
                return circumcircle.lineSegmentsWithAngles.size() < circumcircle.rotation;
            }
        ),
        circumcircles.end()
    );

    // Sorting the circumcircle line segments by angle.
    for (unsigned int i = 0; i < circumcircles.size(); i++) {
        std::sort(
            circumcircles[i].lineSegmentsWithAngles.begin(),
            circumcircles[i].lineSegmentsWithAngles.end(),
            [](const std::pair<LineSegment, double>& pair1, const std::pair<LineSegment, double>& pair2) {
                return pair1.second < pair2.second;
            }
        );
    }

    return circumcircles;
}

// Finding simple rotational symmetries from circumcircles.
std::vector<RotationalSymmetry> RotationalSymmetry::findSimpleSymmetries(const std::vector<Circumcircle>& circumcircles, const double angleTolerance) {
    std::vector<RotationalSymmetry> symmetries;

    // Iterating through each circumcircle.
    for (const Circumcircle& circumcircle : circumcircles) {
        const std::vector<std::pair<LineSegment, double>> lsPairs(circumcircle.lineSegmentsWithAngles);
        std::vector<bool> checked(lsPairs.size(), false);     // Vector with checked property for each line segment.
        const unsigned int rotation = circumcircle.rotation;  // Symmetry rotation level.
        const double desiredAngle = 360.0 / rotation;         // Desired angle in degrees.

        std::vector<unsigned int> currentSequence;

        // Iterating through all lineSegment pairs of line
        // segment and the angle from the vector (1, 0, 0).
        for (unsigned int i = 0; i < lsPairs.size() - 2; i++) {
            // If a line segment has already been checked,
            // there is no need to check it twice.
            if (checked[i]) {
                continue;
            }

            // Adding the new element to the sequence.
            currentSequence.clear();
            currentSequence.push_back(i);
            checked[i] = true;

            // Iterating through the other line segments.
            for (unsigned int j = i + 1; j < lsPairs.size(); j++) {
                // Getting the last stored line segment in the current sequence.
                const unsigned int lastElement = currentSequence.back();

                // Calculating the angles between the last stored and the current line segments.
                double lsAngleFirst = (180 - (180 * VectorFunctions::angle(lsPairs[j].first.p2 - lsPairs[j].first.p1, lsPairs[lastElement].first.p2 - lsPairs[lastElement].first.p1) / PI));
                double lsAngleSecond = 180 - (180 * VectorFunctions::angle(lsPairs[lastElement].first.p2 - lsPairs[lastElement].first.p1, lsPairs[j].first.p2 - lsPairs[j].first.p1) / PI);

                // Our desire is to have angles in range [0, 180]°:
                if (lsAngleFirst < 0) {
                    lsAngleFirst += 180;
                }
                if (lsAngleSecond < 0) {
                    lsAngleSecond += 180;
                }

                // If the angles are alright, a new line segment is added to the current sequence.
                if (Tolerance::isInTolerance(lsPairs[j].second - lsPairs[lastElement].second, desiredAngle, angleTolerance) &&
                    (Tolerance::isInTolerance(lsAngleFirst, desiredAngle, angleTolerance) || Tolerance::isInTolerance(lsAngleSecond, desiredAngle, angleTolerance))
                )
                {
                    currentSequence.push_back(j);  // Adding a new line segment.

                    // If the current sequence size equals the rotation level,
                    // a new symmetry is formed.
                    if (currentSequence.size() == rotation) {
                        // Adding line segments to the symmetry.
                        std::vector<LineSegment> lineSegments;
                        for (const unsigned int element : currentSequence) {
                            lineSegments.push_back(lsPairs[element].first);
                            checked[element] = true;
                        }

                        // Adding a new symmetry.
                        symmetries.push_back(RotationalSymmetry(circumcircle.center, rotation, lineSegments));
                        break;
                    }
                }
            }

            // Clearing the current sequence.
            currentSequence.clear();
        }
    }

    return symmetries;
}

// Merging symmetries by layers, where symmetries with the same rotation axis and rotation level are merged.
std::vector<RotationalSymmetry> RotationalSymmetry::mergeSymmetries(std::vector<RotationalSymmetry>& symmetries) {
    std::vector<RotationalSymmetry> mergedSymmetries;

    // Merging the found symmetries by Z axis and center point.
    for (RotationalSymmetry& symmetry : symmetries) {
        // Searching for a potential symmetry in the list, where X and
        // Y coordinates are the same as in the current symmetry, the
        // symmetry level is also the same.
        auto itr = std::find_if(
            mergedSymmetries.begin(),
            mergedSymmetries.end(),
            [&symmetry](const RotationalSymmetry& curSymmetry) {
                Point axis1 = symmetry.getSymmetryAxis();
                Point axis2 = curSymmetry.getSymmetryAxis();
                return (
                    axis1.x == axis2.x &&
                    axis1.y == axis2.y &&
                    symmetry.getRotation() == curSymmetry.getRotation()
                );
            }
        );

        // Adding the symmetry to the list if not exists.
        if (itr == mergedSymmetries.end()) {
            mergedSymmetries.push_back(symmetry);
        }
        // Updating the Symmetry object with new limits.
        else {
            const unsigned int index = static_cast<unsigned int>(std::distance(mergedSymmetries.begin(), itr));
            mergedSymmetries[index] &= symmetry;
        }
    }

    // Sorting the symmetries according to the rotation level.
    std::sort(
        mergedSymmetries.begin(),
        mergedSymmetries.end(),
        [](const RotationalSymmetry& s1, const RotationalSymmetry& s2) {
            return s1.getRotation() > s2.getRotation();
        }
    );

    return mergedSymmetries;
}

// Getting voxels for each rotational symmetry.
void RotationalSymmetry::getVoxelsInSymmetries(std::vector<RotationalSymmetry>& rotationalSymmetries, const VoxelVector& voxelVector, const VoxelMesh& voxelMesh) {
    // Getting voxels in each rotational symmetry.
    for (unsigned long long i = 0; i < rotationalSymmetries.size(); i++) {
        rotationalSymmetries[i].voxels.clear();  // Clearing the previous possible artefacts in the list. Better safe than sorry.

        // Adding voxels according to reflection symmetry line segments.
        for (const LineSegment& ls : rotationalSymmetries[i].lineSegments) {
            // First line segment point.
            {
                // Getting the line segment point and calculating the voxel coordinates.
                const Point p = ls.p1;
                const int x = static_cast<int>(floor((p.x - voxelMesh.minX) / voxelMesh.sideX));
                const int y = static_cast<int>(floor((p.y - voxelMesh.minY) / voxelMesh.sideY));
                const int z = static_cast<int>(floor((p.z - voxelMesh.minZ) / voxelMesh.sideZ));

                // Creating a new voxel.
                Voxel v(x, y, z, true, true);
                v.inSymmetry = true;

                // If the voxel is not already present in the symmetry, it is added to the list.
                if (std::find(rotationalSymmetries[i].voxels.begin(), rotationalSymmetries[i].voxels.end(), v) == rotationalSymmetries[i].voxels.end()) {
                    rotationalSymmetries[i].voxels.push_back(Voxel(x, y, z, true, true));
                }
            }
            // Second line segment point.
            {
                // Getting the line segment point and calculating the voxel coordinates.
                Point p = ls.p2;
                const int x = static_cast<int>(floor((p.x - voxelMesh.minX) / voxelMesh.sideX));
                const int y = static_cast<int>(floor((p.y - voxelMesh.minY) / voxelMesh.sideY));
                const int z = static_cast<int>(floor((p.z - voxelMesh.minZ) / voxelMesh.sideZ));

                // Creating a new voxel.
                Voxel v(x, y, z, true, true);
                v.inSymmetry = true;

                // If the voxel is not already present in the symmetry, it is added to the list.
                if (std::find(rotationalSymmetries[i].voxels.begin(), rotationalSymmetries[i].voxels.end(), v) == rotationalSymmetries[i].voxels.end()) {
                    rotationalSymmetries[i].voxels.push_back(Voxel(x, y, z, true, true));
                }
            }
        }

        // Adding voxels that lie at the axis.
        const Point p = rotationalSymmetries[i].symmetryAxis;
        const unsigned int x = static_cast<unsigned int>(floor((p.x - voxelMesh.minX) / voxelMesh.sideX));
        const unsigned int y = static_cast<unsigned int>(floor((p.y - voxelMesh.minY) / voxelMesh.sideY));

        for (unsigned int z = 0; z < static_cast<unsigned int>(voxelMesh.zSize); z++) {
            if (!voxelVector[z][y][x].material) {
                continue;
            }

            // Creating a new voxel.
            Voxel v(static_cast<int>(x), static_cast<int>(y), static_cast<int>(z), true, true);
            v.inSymmetry = true;

            // If the voxel is not already present in the symmetry, it is added to the list.
            if (std::find(rotationalSymmetries[i].voxels.begin(), rotationalSymmetries[i].voxels.end(), v) == rotationalSymmetries[i].voxels.end()) {
                rotationalSymmetries[i].voxels.push_back(Voxel(static_cast<int>(x), static_cast<int>(y), static_cast<int>(z), true, true));
            }
        }
    }

    // Sorting the symmetries according to the voxel count.
    std::sort(
        rotationalSymmetries.begin(),
        rotationalSymmetries.end(),
        [](const RotationalSymmetry& s1, const RotationalSymmetry& s2) {
            return s1.voxels.size() > s2.voxels.size();
        }
    );
}



// OBJECT METHODS
// Main constructor.
RotationalSymmetry::RotationalSymmetry(const Point& symmetryAxis, const unsigned int rotation, const std::vector<LineSegment>& lineSegments) :
    rotation(rotation),
    symmetryAxis(symmetryAxis),
    lineSegments(lineSegments)
{}

// Calculating the point positions according to the axis.
std::vector<PointPosition> RotationalSymmetry::calculatePositionsFromAxis(const std::vector<Point>& points, const VoxelMesh& voxelMesh, const bool validAxis) {
    const unsigned long long numberOfPoints = points.size();  // Getting the point vector size.

    // If there is no plane, all points are undefined.
    if (!validAxis) {
        return std::vector<PointPosition>(numberOfPoints, PointPosition::notInSymmetry);
    }

    // Creating a new vector with positions.
    std::vector<PointPosition> positions(numberOfPoints, PointPosition::notInSymmetry);
    for (unsigned long long i = 0; i < numberOfPoints; i++) {
        // Getting each single point and its voxel coordinates.
        const Point& point = points[i];
        const int voxelX = static_cast<int>(floor((point.x - voxelMesh.minX) / static_cast<float>(voxelMesh.sideX)));
        const int voxelY = static_cast<int>(floor((point.y - voxelMesh.minY) / static_cast<float>(voxelMesh.sideY)));
        const int voxelZ = static_cast<int>(floor((point.z - voxelMesh.minZ) / static_cast<float>(voxelMesh.sideZ)));

        // If the current point lies outside of an symmetry voxel, its position is set to undefined.
        if (std::find(voxels.begin(), voxels.end(), Voxel(voxelX, voxelY, voxelZ)) == voxels.end()) {
             positions[i] = PointPosition::notInSymmetry;
             continue;
        }

        // Otherwise, the position is set to rotational.
        positions[i] = PointPosition::rotational;
    }

    return positions;
}



// GETTERS
// Symmetry axis getter.
Point RotationalSymmetry::getSymmetryAxis() const {
    return symmetryAxis;
}

// Rotation getter.
unsigned int RotationalSymmetry::getRotation() const {
    return rotation;
}

// Line segments getter.
const std::vector<LineSegment>& RotationalSymmetry::getLineSegments() const {
    return lineSegments;
}

// Voxel getter.
const std::vector<Voxel>& RotationalSymmetry::getVoxels() const {
    return voxels;
}



// OVERLOADED OPERATORS
// Operator &= is used for merging the current symmetry with another.
RotationalSymmetry RotationalSymmetry::operator &= (RotationalSymmetry& s) {
    lineSegments.insert(lineSegments.end(), s.lineSegments.begin(), s.lineSegments.end());
    voxels.insert(voxels.end(), s.voxels.begin(), s.voxels.end());

    return *this;
}



// CLASS METHODS
// Calculating rotational symmetries.
std::vector<RotationalSymmetry> RotationalSymmetry::calculateRotationalSymmetries(
    const std::vector<Point>& points,
    const VoxelVector& voxelVector,
    const VoxelMesh& voxelMesh,
    const double lengthTolerance,
    const double angleTolerance,
    const double minimumLineSegmentDistance
)
{
    std::vector<Point> voxelPoints = pointsFromInterestingVoxels();
    std::vector<LineSegment> lineSegments = calculateLineSegmentsBetweenPoints(voxelPoints, lengthTolerance, minimumLineSegmentDistance, 0.5);
    SplitLineSegmentPairs splitLineSegments = splitLineSegmentsInPairs(lineSegments, voxelMesh, lengthTolerance, angleTolerance);
    std::vector<Circumcircle> circumcircles = buildLineSegmentMidpointsCircumcircles(splitLineSegments);
    std::vector<RotationalSymmetry> symmetries = findSimpleSymmetries(circumcircles, angleTolerance);
    std::vector<RotationalSymmetry> mergedSymmetries = mergeSymmetries(symmetries);
    getVoxelsInSymmetries(mergedSymmetries, voxelVector, voxelMesh);

    return mergedSymmetries;
}





/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
// LINESEGMENTGRAPH
// CLASS METHODS
// Main constructor.
LineSegmentGraph::LineSegmentGraph(const std::vector<LineSegmentPair>& lineSegmentPairs, const unsigned int rotation) :
    rotation(static_cast<int>(rotation))
{
    buildGraphFromLineSegmentsPairs(lineSegmentPairs);  // Building a graph.
}

// Building a graph from line segment pairs.
void LineSegmentGraph::buildGraphFromLineSegmentsPairs(const std::vector<LineSegmentPair>& lineSegmentPairs) {
    // Initializing helper structures for building a graph.
    std::vector<LineSegment> lineSegments;  // Declaring a line segment vector (future vertices).
    std::vector<std::vector<bool>> connections(2 * lineSegmentPairs.size(), std::vector<bool>(2 * lineSegmentPairs.size(), false));  // Initializing helper connection 2D vector (future edges).

    // Iterating through each pair of line segments.
    for (const LineSegmentPair& pair : lineSegmentPairs) {
        // Getting the both line segments from a pair.
        const LineSegment firstLineSegment = pair.ls1;
        const LineSegment secondLineSegment = pair.ls2;

        // Searching whether each separate line segment is already present in the vector.
        auto firstIterator = std::find(lineSegments.begin(), lineSegments.end(), firstLineSegment);
        auto secondIterator = std::find(lineSegments.begin(), lineSegments.end(), secondLineSegment);

        // Checking whether the iterators are at the end.
        // If not, something has been found :)
        bool firstFound = firstIterator != lineSegments.end();
        bool secondFound = secondIterator != lineSegments.end();

        // Getting the indices of the both line segments in the vector.
        unsigned int firstIndex = static_cast<unsigned int>(std::distance(lineSegments.begin(), firstIterator));
        unsigned int secondIndex = static_cast<unsigned int>(std::distance(lineSegments.begin(), secondIterator));

        // If the first line segment has not been found, we add it
        // to the vector and update the first index accordingly.
        if (!firstFound) {
            lineSegments.push_back(firstLineSegment);
            firstIndex = static_cast<unsigned int>(lineSegments.size()) - 1;
        }

        // If the second line segment has not been found, we add it
        // to the vector and update the second index accordingly.
        if (!secondFound) {
            lineSegments.push_back(secondLineSegment);
            secondIndex = static_cast<unsigned int>(lineSegments.size()) - 1;
        }

        // Setting the connection (relation or edge)
        // between the two line segments to True.
        connections[firstIndex][secondIndex] = true;
        connections[secondIndex][firstIndex] = true;
    }

    // After having completed the unique line segment search,
    // graph vertices are created from those line segments.
    for (unsigned int i = 0; i < lineSegments.size(); i++) {
        vertices.push_back(Vertex(i, lineSegments[i]));
    }

    // After creating graph vertices, edges are set. In this case,
    // edges represent relations between pairs of line segments.
    edges = std::vector<std::vector<bool>>(lineSegments.size(), std::vector<bool>(lineSegments.size(), false));
    for (unsigned int i = 0; i < lineSegments.size(); i++) {
        for (unsigned int j = 0; j < lineSegments.size(); j++) {
            edges[i][j] = connections[i][j];
        }
    }
}


// Vertices getter.
std::vector<Vertex> LineSegmentGraph::getVertices() const {
    return vertices;
}

// Edges getter.
std::vector<std::vector<bool>> LineSegmentGraph::getEdges() const {
    return edges;
}

// Rotation getter.
int LineSegmentGraph::getRotation() const {
    return rotation;
}


/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
// VERTEX
// OBJECT METHODS
// Main constructor.
Vertex::Vertex(const unsigned int index, const LineSegment& lineSegment) :
    index(index),
    checked(false),
    lineSegment(lineSegment)
{}


// Index getter.
unsigned int Vertex::getIndex() const {
    return index;
}

// Checked getter.
bool Vertex::getChecked() const {
    return checked;
}

// Checked setter (true).
void Vertex::setChecked() {
    checked = true;
}

// Line segment getter.
LineSegment Vertex::getLineSegment() const {
    return lineSegment;
}

// Inverted line segment getter.
LineSegment Vertex::getInvertedLineSegment() const {
    auto P2 = lineSegment.p2;
    auto P1 = lineSegment.p1;

    return LineSegment(P2, P1);
}


/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
// CIRCUMCIRCLE
// Main constructor.
Circumcircle::Circumcircle(const unsigned int rotation, const Point& center, const double radius, const std::vector<bool>& lineSegmentsInCircumcircle, const std::vector<LineSegment>& lineSegments) :
    rotation(rotation),
    radius(radius),
    center(center),
    lineSegmentsInCircumcircle(lineSegmentsInCircumcircle)
{
    // Adding line segments.
    for (const LineSegment& ls : lineSegments) {
        addLineSegment(ls);
    }
}

// Adding a new line segment.
void Circumcircle::addLineSegment(const LineSegment& ls) {
    const Vector3d v = (ls.midpoint() - center).normalize();
    const Vector3d vb(1, 0, 0);

    const double angle = 180 * VectorFunctions::angle(vb, v) / PI;

    lineSegmentsWithAngles.push_back(std::make_pair(ls, angle));
}


/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
// LINE SEGMENT PAIR
// Main constructor.
LineSegmentPair::LineSegmentPair(const LineSegment& ls1, const LineSegment& ls2) :
    ls1(ls1),
    ls2(ls2)
{}


/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
// LINE SEGMENT PAIR GROUP
// Main constructor.
LineSegmentPairGroup::LineSegmentPairGroup(const std::vector<LineSegmentPair>& pairs, const unsigned int rotation, unsigned int level, unsigned int layer) :
    pairs(pairs),
    rotation(rotation),
    level(level),
    layer(layer)
{
    // If there are some pairs, the distance between the center points is calculated.
    if (!pairs.empty()) {
        distanceBetweenPairCenterPoints = PointFunctions::distance(pairs[0].ls1.midpoint(), pairs[0].ls2.midpoint());
    }
}

// Adding a new pair to the group.
void LineSegmentPairGroup::addPair(const LineSegmentPair& pair) {
    pairs.push_back(pair);
    if (distanceBetweenPairCenterPoints < 0) {
        distanceBetweenPairCenterPoints = PointFunctions::distance(pairs[0].ls1.midpoint(), pairs[0].ls2.midpoint());
    }
}


/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
/**********************************************************************************************************************************************************/
// LINE SEGMENT SUBGROUP
// Main constructor.
LineSegmentSubgroup::LineSegmentSubgroup(const std::vector<LineSegment>& lineSegments, const unsigned int rotation) :
    rotation(rotation),
    lineSegments(lineSegments)
{}
