#include <algorithm>
#include <limits>
#include <map>
#include <QString>
#include <queue>

#include "BruteRotationalSymmetry.hpp"
#include "Structs/Constants.hpp"
#include "Structs/Point.tpp"
#include "Structs/Voxel.hpp"

using namespace Symmetry;



// PRIVATE OBJECT METHODS
// Function that increments the radius of the circumference around the center point
// with the step of 1. In every step it selects all interesting voxels the circumference
// overlaps with and calls the function FindRotationalSymmetries.
// In the end a union of all interesting voxels from obtained symmetries with different
// radii is made. The function returns a list of rotational symmetries (2 through 16),
// where each rotational symmetry is present once at most.
std::vector<BruteRotationalSymmetry> BruteRotationalSymmetry::getRotationalSymmetriesAroundCenterPoint(
    std::vector<std::vector<Voxel>>& field,
    Point centerPoint,
    const int maxRadius,
    int pixelX,
    int pixelY)
{
    std::vector<BruteRotationalSymmetry> tempSymmetries, symmetries, outSymmetries;
    std::queue<Voxel*> qu;
    std::vector<Voxel> interestingVoxels;
    Voxel* tp;
    Voxel* tp1;
    std::map<Voxel, bool> voxelsInSymmetry;

    // Maximum radius is the distance from the voxel to the nearest scene edge.
    int radius;
    double minD = 0.0;
    double maxD = 0.0;

    // Radius incrementation by 1.
    for (radius = 1; radius <= maxRadius; radius++) {
        interestingVoxels.clear();

        // Adding first voxel to the queue (same Y coordinate, for radius
        // away from the center point according to X coordinate.
        tp = &field[static_cast<unsigned int>(centerPoint.y)][static_cast<unsigned int>(centerPoint.x) + radius];
        if (tp->interesting) {
            interestingVoxels.push_back(*tp);
        }
        tp->checked = true;
        qu.push(tp);

        // Four neighboring voxels of each voxel in queue are checked
        // (if not checked before and if the circumference is overlapping with them).
        // If true, voxels are added to the queue and set checked. If also interesting,
        // they are pushed to interestingVoxels.
        while (!qu.empty()) {
            tp = qu.front();

            // Left voxel.
            if (tp->x > 0) {
                tp1 = &field[static_cast<unsigned int>(tp->y)][static_cast<unsigned int>(tp->x - 1)];
                if (!tp1->checked) {
                    BruteRotationalSymmetry::getDistanceRanges(centerPoint, *tp1, minD, maxD);
                    if (minD <= double(radius) && maxD >= double(radius)) {
                        qu.push(tp1);
                        tp1->checked = true;
                        if (tp1->interesting) {
                            interestingVoxels.push_back(*tp1);
                        }
                    }
                }
            }

            // Bottom voxel.
            if (tp->y > 0) {
                tp1 = &field[static_cast<unsigned int>(tp->y - 1)][static_cast<unsigned int>(tp->x)];
                if (!tp1->checked) {
                    BruteRotationalSymmetry::getDistanceRanges(centerPoint, *tp1, minD, maxD);
                    if (minD <= double(radius) && maxD >= double(radius)) {
                        qu.push(tp1);
                        tp1->checked = true;
                        if (tp1->interesting) {
                            interestingVoxels.push_back(*tp1);
                        }
                    }
                }
            }

            // Right voxel.
            if (tp->x < static_cast<uint>(pixelX - 1)) {
                tp1 = &field[static_cast<unsigned int>(tp->y)][static_cast<unsigned int>(tp->x + 1)];
                if (!tp1->checked) {
                    BruteRotationalSymmetry::getDistanceRanges(centerPoint, *tp1, minD, maxD);
                    if (minD <= double(radius) && maxD >= double(radius)) {
                        qu.push(tp1);
                        tp1->checked = true;
                        if (tp1->interesting) {
                            interestingVoxels.push_back(*tp1);
                        }
                    }
                }
            }

            // Upper voxel.
            if (tp->y < static_cast<uint>(pixelY - 1)) {
                tp1 = &field[static_cast<unsigned int>(tp->y + 1)][static_cast<unsigned int>(tp->x)];
                if (!tp1->checked) {
                    BruteRotationalSymmetry::getDistanceRanges(centerPoint, *tp1, minD, maxD);
                    if (minD <= double(radius) && maxD >= double(radius)) {
                        qu.push(tp1);
                        tp1->checked = true;
                        if (tp1->interesting) {
                            interestingVoxels.push_back(*tp1);
                        }
                    }
                }
            }

            qu.pop();
        }

        // Resetting all checked voxels to false.
        tp = &field[static_cast<unsigned int>(centerPoint.y)][static_cast<unsigned int>(centerPoint.x) + radius];
        tp->checked= false;
        qu.push(tp);
        while (!qu.empty()) {
            tp = qu.front();
            qu.pop();

            // Left voxel.
            if (tp->x > 0) {
                tp1 = &field[static_cast<unsigned int>(tp->y)][static_cast<unsigned int>(tp->x - 1)];
                if (tp1->checked) {
                    tp1->checked = false;
                    qu.push(tp1);
                }
            }

            // Bottom voxel.
            if (tp->y > 0) {
                tp1 = &field[static_cast<unsigned int>(tp->y - 1)][static_cast<unsigned int>(tp->x)];
                if (tp1->checked) {
                    tp1->checked = false;
                    qu.push(tp1);
                }
            }

            // Right voxel.
            if (tp->x < static_cast<uint>(pixelX - 1)) {
                tp1 = &field[static_cast<unsigned int>(tp->y)][static_cast<unsigned int>(tp->x + 1)];
                if (tp1->checked) {
                    tp1->checked = false;
                    qu.push(tp1);
                }
            }

            // Upper voxel.
            if (tp->y < static_cast<uint>(pixelY - 1)) {
                tp1 = &field[static_cast<unsigned int>(tp->y + 1)][static_cast<unsigned int>(tp->x)];
                if (tp1->checked) {
                    tp1->checked = false;
                    qu.push(tp1);
                }
            }
        }

        // Test of rotation symmetries in the set of interesting
        // voxels on current radius.
        if (interestingVoxels.size() > 1) {
            tempSymmetries = getRotationalSymmetriesOnRadius(
                centerPoint,
                interestingVoxels
            );


            // Appending found symmetries to the list of all symmetries.
            for (unsigned int i = 0; i < tempSymmetries.size(); i++) {
                symmetries.push_back(tempSymmetries[i]);
            }
        }
    }

    // Sorting of symmetries by the rotation level.
    std::sort(
        symmetries.begin(),
        symmetries.end(),
        [](const auto& lhs, const auto& rhs) {
            return lhs.getRotation() < rhs.getRotation();
        }
    );

    int currentRotation = 0;
    int cosi = -1;

    // Walking through the list of all symmetries and generating a symmetry for each rotation
    // (voxel union in symmetries, as some voxels may be present also in neighbouring symmetries
    // with the same rotation (two circumferences can overlap with the same voxel)).
    for (unsigned int i = 0; i < symmetries.size(); i++) {
        if (symmetries[i].getRotation() > currentRotation) {
            // If the rotation of the current symmetry is bigger that from the previous one,
            // a new symmetry with the given rotation is created.
            voxelsInSymmetry.clear();
            currentRotation = symmetries[i].getRotation();
            outSymmetries.push_back(symmetries[i]);
            cosi++;
            for (int j = 0; j < symmetries[i].getVoxelCount(); j++) {
                voxelsInSymmetry.insert(
                    std::pair<Voxel, bool>(symmetries[i].getVoxel(j), true)
                );
            }
        }
        else {
            // Merging symmetry with the previous one.
            for (int j = 0; j < symmetries[i].getVoxelCount(); j++) {
                // If there is no voxel in the symmetry, it is added to it.
                if (voxelsInSymmetry.count(symmetries[i].getVoxel(j)) == 0) {
                    voxelsInSymmetry.insert(
                        std::pair<Voxel, bool>(symmetries[i].getVoxel(j), true)
                    );
                    outSymmetries[static_cast<unsigned int>(cosi)].addVoxel(symmetries[i].getVoxel(j));
                }
            }
        }
    }

    return outSymmetries;
}

// Merging symmetries by layers. Symmetries with the same search axis and rotation level are merged.
std::vector<BruteRotationalSymmetry> BruteRotationalSymmetry::mergeSymmetries(std::vector<BruteRotationalSymmetry>& symmetries) {
    std::vector<BruteRotationalSymmetry> mergedSymmetries;

    // Merging the found symmetries by Z axis and center point.
    for (BruteRotationalSymmetry& symmetry : symmetries) {
        // Searching for a potential symmetry in the list, where X and
        // Y coordinates are the same as in the current symmetry, the
        // symmetry level is also the same.
        auto itr = std::find_if(
            mergedSymmetries.begin(),
            mergedSymmetries.end(),
            [&symmetry](const BruteRotationalSymmetry curSymmetry) {
                Point axis1 = symmetry.getSymmetryAxis();
                Point axis2 = curSymmetry.getSymmetryAxis();
                return (
                    axis1.x == axis2.x
                    && axis1.y == axis2.y
                    && symmetry.getRotation() == curSymmetry.getRotation());
            }
        );

        // Adding the symmetry to the list if not exists.
        if (itr == mergedSymmetries.end()) {
            mergedSymmetries.push_back(symmetry);
        }
        // Updating the Symmetry object with new limits.
        else {
            const unsigned int index = std::distance(mergedSymmetries.begin(), itr);
            mergedSymmetries[index] &= symmetry;
        }
    }

    // Sorting the symmetries according to the rotation level.
    std::sort(
        mergedSymmetries.begin(),
        mergedSymmetries.end(),
        [](const BruteRotationalSymmetry& s1, const BruteRotationalSymmetry& s2) {
            return s1.getRotation() > s2.getRotation();
        }
    );

    return mergedSymmetries;
}

// A main function for searching symmetries gets a center point and a list of interesting voxels.
// Symmetries are searched for only in the set of interesting voxels that are on the same radius
// from the center point (or inside the same interval).
std::vector<BruteRotationalSymmetry> BruteRotationalSymmetry::getRotationalSymmetriesOnRadius(
    Point centerPoint,
    const std::vector<Voxel>& interestingVoxels)
{
    int maxSymmetrySize = 16;               // Maximum symmetry level of the symmetry.
    std::vector<AngleRangeEvent> angleEvents;    // Vector of angle events.
    std::vector<int> activeLevels;          // Active levels.
    int activeLevelsCount;                  // Number of active levels.
    std::map<Voxel, int> activeVoxels;    // Map of active voxels.
    std::map<Voxel, int>::iterator pitr;  // Iterator through map.
    std::vector<bool> voxelsInSymmetry;          // Vector of voxels in symmetry.
    int voxelsInSymmetryCount = 0;          // Number of voxels in symmetry.
    std::vector<BruteRotationalSymmetry> symmetries;  // List of different rotational symmetries.
    BruteRotationalSymmetry tempSymmetry;             // Current symmetry.

    // Filling the vector of active levels with zeros.
    for (int i = 0; i < maxSymmetrySize; i++) {
        activeLevels.push_back(0);
    }

    // Filling the map with the pairs of interesting voxels and zeros.
    // Currently there have been no found voxels in the symmetry.
    for (unsigned int i = 0; i < interestingVoxels.size(); i++) {
        activeVoxels.insert(std::pair<Voxel, int>(interestingVoxels[i], 0));
        voxelsInSymmetry.push_back(false);
    }

    std::vector<BruteRotationalSymmetry> tempSymmetries;                // Vector of currently found symmetries.
    std::vector<AngleRangeEvent> originalAngleEvents;    // Vector of current angle events.
    float minAngle = std::numeric_limits<float>::max();  // Minimum angle.
    float maxAngle = std::numeric_limits<float>::min();  // Maximum angle.

    // Appending angle events for all interesting voxels.
    for (unsigned int i = 0; i < interestingVoxels.size(); i++) {
        // Minimum and maximum angle calculation.
        BruteRotationalSymmetry::getAngleRanges(centerPoint, interestingVoxels[i], minAngle, maxAngle);

        // If the angle is >PI, the interval is split into 2 parts.
        if (maxAngle - minAngle > PI) {
            originalAngleEvents.push_back(
                AngleRangeEvent(0, 0, true, interestingVoxels[i])
            );
            originalAngleEvents.push_back(
                AngleRangeEvent(minAngle, 0, false, interestingVoxels[i])
            );
            originalAngleEvents.push_back(
                AngleRangeEvent(maxAngle, 0, true, interestingVoxels[i])
            );
            originalAngleEvents.push_back(
                AngleRangeEvent(2 * PI, 0, false, interestingVoxels[i])
            );
        }
        else {
            originalAngleEvents.push_back(
                AngleRangeEvent(minAngle, 0, true, interestingVoxels[i])
            );
            originalAngleEvents.push_back(
                AngleRangeEvent(maxAngle, 0, false, interestingVoxels[i])
            );
        }
    }

    // Sorting events by angle.
    std::sort(
        originalAngleEvents.begin(),
        originalAngleEvents.end(),
        [](auto& lhs, auto& rhs) {
            return lhs.angle < rhs.angle;
        }
    );

    std::map<Voxel, bool> activeInterval;  // An active interval.
    std::map<Voxel, bool>::iterator itr;   // Iterator through intervals.

    // Appending the active intervals for each interesting voxel.
    for (unsigned int i = 0; i < interestingVoxels.size(); i++) {
        activeInterval.insert({std::make_pair(interestingVoxels[i], false)});
    }

    int level = 0;     // Symmetry level.
    int newLevel = 0;  // New symmetry level.
    double rotationAngle;   // Angle of the rotation.

    // Marching through all the symmetries.
    for (int i = 2; i <= maxSymmetrySize; i++) {
        // Cleaning variables from the previous loop iteration.
        for (itr = activeInterval.begin(); itr != activeInterval.end(); ++itr) {
            itr->second = false;
        }
        angleEvents.clear();
        level = 0;
        rotationAngle = 2 * PI / i;  // Calculation of the rotation symmetry angle.

        // Creating intervals on all levels.
        for (unsigned int j = 0; j < originalAngleEvents.size(); j++) {
            newLevel = int(std::floor(originalAngleEvents[j].angle / rotationAngle));
            if (newLevel > 0
                && !originalAngleEvents[j].start
                && originalAngleEvents[j].angle / rotationAngle - double(newLevel) == 0.0)
            {
                newLevel--;
            }

            if (newLevel > level) {
                // Initialization of the new levels start.
                level = newLevel;
                for (itr = activeInterval.begin(); itr != activeInterval.end(); ++itr) {
                    if (itr->second) {
                        angleEvents.push_back(AngleRangeEvent(0, level, true, itr->first));
                    }
                }
            }

            // Adding new angle event.
            angleEvents.push_back(
                AngleRangeEvent(
                    originalAngleEvents[j].angle - (newLevel * rotationAngle),
                    newLevel,
                    originalAngleEvents[j].start,
                    originalAngleEvents[j].voxel
                )
            );
            activeInterval[originalAngleEvents[j].voxel] = originalAngleEvents[j].start;
        }

        // Sorting all angle events by angle and start.
        std::sort(angleEvents.begin(), angleEvents.end(), [](auto& lhs, auto& rhs) {
            return lhs.start < rhs.start;
        });
        std::sort(angleEvents.begin(), angleEvents.end(), [](auto& lhs, auto& rhs) {
            return lhs.angle < rhs.angle;
        });

        // Checking voxels that are present in intervals.
        for (unsigned int j = 0; j < static_cast<unsigned int>(maxSymmetrySize); j++) {
            activeLevels[j] = 0;
        }
        for (unsigned int j = 0; j < interestingVoxels.size(); j++) {
            voxelsInSymmetry[j] = false;
        }
        for (pitr = activeVoxels.begin(); pitr != activeVoxels.end(); pitr++) {
            pitr->second = 0;
        }
        activeLevelsCount = 0;
        voxelsInSymmetryCount = 0;

        for (unsigned int j = 0; j < angleEvents.size(); j++) {
            if (angleEvents[j].start) {
                activeLevels[static_cast<unsigned int>(angleEvents[j].level)]++;
                activeVoxels[angleEvents[j].voxel]++;
                if (activeLevels[static_cast<unsigned int>(angleEvents[j].level)] == 1) {
                    activeLevelsCount++;
                }
                if (activeLevelsCount == i) {
                    for (unsigned long long k = 0; k < voxelsInSymmetry.size(); k++) {
                        if (voxelsInSymmetry[k] == false)
                            voxelsInSymmetryCount++;
                        voxelsInSymmetry[k] = true;
                    }
                }
            }
            else {
                activeLevels[static_cast<unsigned int>(angleEvents[j].level)]--;
                activeVoxels[angleEvents[j].voxel]--;
                if (activeLevels[static_cast<unsigned int>(angleEvents[j].level)] == 0) {
                    activeLevelsCount--;
                }
            }
        }

        // If no symmetry, a next iteration of the loop follows.
        if (voxelsInSymmetryCount == 0) {
            continue;
        }

        // Creating the rotational symmetry and appending to the
        // list of all symmetries on current radius.
        const int layer = centerPoint.z;
        tempSymmetry = BruteRotationalSymmetry(centerPoint, i, layer, layer);

        // Adding voxels to the current symmetry.
        for (unsigned int j = 0; j < voxelsInSymmetry.size(); j++) {
            if (voxelsInSymmetry[j]) {
                tempSymmetry.addVoxel(interestingVoxels[j]);
            }
        }

        // Adding symmetry to the list of all symmetries.
        symmetries.push_back(tempSymmetry);
    }

    return symmetries;
}



// CONSTRUCTOR AND DESTRUCTOR
// Default constructor.
BruteRotationalSymmetry::BruteRotationalSymmetry() :
    rotation(0),
    upperLayer(0),
    lowerLayer(0),
    symmetryAxis(Point())
{}

// Constructor with all parameters.
BruteRotationalSymmetry::BruteRotationalSymmetry(
    const Point& symmetryAxis,
    const int& rotation,
    const int& upperLayer,
    const int& lowerLayer
) :
    rotation(rotation),
    upperLayer(upperLayer),
    lowerLayer(lowerLayer),
    symmetryAxis(symmetryAxis)
{}



// OVERLOADED OPERATORS
// Operator & is used for merging two symmetries.
BruteRotationalSymmetry BruteRotationalSymmetry::operator & (const BruteRotationalSymmetry& s) const {
    BruteRotationalSymmetry newS = BruteRotationalSymmetry(*this);   // New symmetry.

    // Searching for minimum and maximum layer inside the both symmetries.
    int minLayer = std::min(this->getLowerLayer(), s.getLowerLayer());
    int maxLayer = std::max(this->getUpperLayer(), s.getUpperLayer());

    // Setting new upper and lower layer.
    newS.setLowerLayer(minLayer);
    newS.setUpperLayer(maxLayer);

    return newS;
}

// Operator &= is used for merging the current symmetry with another.
BruteRotationalSymmetry BruteRotationalSymmetry::operator &= (BruteRotationalSymmetry& s) {
    // Searching for minimum and maximum layer inside the both symmetries.
    int minLayer = std::min(this->getLowerLayer(), s.getLowerLayer());
    int maxLayer = std::max(this->getUpperLayer(), s.getUpperLayer());

    // Setting new upper and lower layer.
    this->lowerLayer = minLayer;
    this->upperLayer = maxLayer;

    // Adding the voxels of the second symmetry.
    std::vector<Voxel> voxels = s.getVoxels();
    for (unsigned int i = 0; i < voxels.size(); i++) {
        this->addVoxel(voxels[i]);
    }

    return* this;
}



// GETTERS AND SETTERS
// Getting voxel count.
int BruteRotationalSymmetry::getVoxelCount() const {
    return (int)voxels.size();
}

// Getting all voxels.
std::vector<Voxel>& BruteRotationalSymmetry::getVoxels() {
    return voxels;
}

// Getting a voxel at a certain index.
Voxel BruteRotationalSymmetry::getVoxel(unsigned int index) const {
    return voxels[index];
}

// Adding a voxel to the vector.
void BruteRotationalSymmetry::addVoxel(Voxel v) {
    voxels.push_back(v);
}

// Getting a symmetry axis.
Point BruteRotationalSymmetry::getSymmetryAxis() const {
    return symmetryAxis;
}

// Setting a symmetry axis.
void BruteRotationalSymmetry::setSymmetryAxis(Point symmetryAxis) {
    this->symmetryAxis = symmetryAxis;
}

// Getting a rotation.
int BruteRotationalSymmetry::getRotation() const {
    return rotation;
}

// Setting a rotation.
void BruteRotationalSymmetry::setRotation(int rotation) {
    this->rotation = rotation;
}

// Getting an upper layer.
int BruteRotationalSymmetry::getUpperLayer() const {
    return upperLayer;
}

// Setting an upper layer.
void BruteRotationalSymmetry::setUpperLayer(int upperLayer) {
    this->upperLayer = upperLayer;
}

// Getting a lower layer.
int BruteRotationalSymmetry::getLowerLayer() const {
    return lowerLayer;
}

// Setting a lower layer.
void BruteRotationalSymmetry::setLowerLayer(int lowerLayer) {
    this->lowerLayer = lowerLayer;
}



// CLASS METHODS
// Symmetry search class initialization method.
void BruteRotationalSymmetry::setGeometrySearch(int x, int y) {
    // Creating 2D array for vicinity scan.
    for (int i = 0; i < x; i++) {
        BruteRotationalSymmetry::geometrySearch.push_back(std::vector<bool>());

        for (int j = 0; j < y; j++) {
            BruteRotationalSymmetry::geometrySearch[i].push_back(false);
        }
    }
    BruteRotationalSymmetry::geometrySearchSizeX = x;  // Size of scan by X.
    BruteRotationalSymmetry::geometrySearchSizeY = y;  // Size of scan by Y.
}

// Cleaning the history of geometry search.
void BruteRotationalSymmetry::clearGeometrySearch(int x, int y) {
    std::queue<int> coordsQueue;  // Queue for the search.
    geometrySearch[x][y] = false;
    coordsQueue.push(x);
    coordsQueue.push(y);

    // Scanning the vicinity.
    while (!coordsQueue.empty()) {
        x = coordsQueue.front();
        coordsQueue.pop();

        y = coordsQueue.front();
        coordsQueue.pop();

        // Left point.
        if (x > 0) {
            if (geometrySearch[x - 1][y]) {
                geometrySearch[x - 1][y] = false;
                coordsQueue.push(x - 1);
                coordsQueue.push(y);
            }
        }
        // Right point.
        if (x < geometrySearchSizeX - 1) {
            if (geometrySearch[x + 1][y]) {
                geometrySearch[x + 1][y] = false;
                coordsQueue.push(x + 1);
                coordsQueue.push(y);
            }
        }
        // Bottom point.
        if (y > 0) {
            if (geometrySearch[x][y - 1]) {
                geometrySearch[x][y - 1] = false;
                coordsQueue.push(x);
                coordsQueue.push(y - 1);
            }
        }
        // Top point.
        if (y < geometrySearchSizeY - 1) {
            if (geometrySearch[x][y + 1]) {
                geometrySearch[x][y + 1] = false;
                coordsQueue.push(x);
                coordsQueue.push(y + 1);
            }
        }
    }
}

// Calculation of angles of rotational symmetries (interval [1, 16]).
void BruteRotationalSymmetry::setSymmetryAngles() {
    symmetryAngles.push_back(FP_INFINITE);

    // Angle calculation of a rotational symmetry.
    for (int i = 1; i < 16; i++) {
        symmetryAngles.push_back((float)(2 * PI / i));
    }
}

// Calculation of the minimum and the maximum distance from the center point to the voxel (in 2D).
void BruteRotationalSymmetry::getDistanceRanges(Point centerPoint, Voxel v, double& minD, double& maxD) {
    minD = 1000000.0;
    maxD = 0.0;
    Point dp(v.x, v.y, v.z);

    double D = PointFunctions::distance(centerPoint, dp);
    if (D < minD) {
        minD = D;
    }
    if (D > maxD) {
        maxD = D;
    }
    dp.x = dp.x + 1.0;

    D = PointFunctions::distance(centerPoint, dp);
    if (D < minD) {
        minD = D;
    }
    if (D > maxD) {
        maxD = D;
    }
    dp.y = dp.y + 1.0;

    D = PointFunctions::distance(centerPoint, dp);
    if (D < minD) {
        minD = D;
    }
    if (D > maxD) {
        maxD = D;
    }
    dp.x = dp.x - 1.0;

    D = PointFunctions::distance(centerPoint, dp);
    if (D < minD) {
        minD = D;
    }
    if (D > maxD) {
        maxD = D;
    }
}

// Getting the angle between two vectors (represented as Points).
float BruteRotationalSymmetry::getAngle(const Point a, const Point b) {
    float dot = (float)(a.x * b.x) + (float)(a.y * b.y);
    float det = (float)(a.x * b.y) - (float)(a.y * b.x);
    float angle = std::atan2(det, dot);

    // Angle lies inside the interval [0, 2 * PI].
    return (float)(angle >= 0 ? angle : angle + PI * 2);
}

// Calculation of the minimum and the maximum angle between horizontal vector
// and the vectors from the center point to the vertices of the voxel (in 2D).
void BruteRotationalSymmetry::getAngleRanges(
    const Point& centerPoint,
    Voxel v,
    float& minAngle,
    float& maxAngle)
{
    minAngle = (float)(2 * PI);
    maxAngle = 0.0;
    Point vertical(1.0, 0.0, 0.0);
    Point dp(v.x, v.y, v.z);


    if (double(v.x) > centerPoint.x
        && centerPoint.y - double(v.y) < 1.0
        && centerPoint.y - double(v.y) > 0.0)
    {
        maxAngle = getAngle(vertical, dp - centerPoint);
        dp.y = dp.y + 1;
        minAngle = getAngle(vertical, dp - centerPoint);
        return;
    }

    float angle = getAngle(vertical, dp - centerPoint);
    if (angle < minAngle)
        minAngle = angle;
    if (angle > maxAngle)
        maxAngle = angle;
    dp.x = dp.x + 1.0;

    angle = getAngle(vertical, dp - centerPoint);
    if (angle < minAngle)
        minAngle = angle;
    if (angle > maxAngle)
        maxAngle = angle;
    dp.y = dp.y + 1.0;

    angle = getAngle(vertical, dp - centerPoint);
    if (angle < minAngle)
        minAngle = angle;
    if (angle > maxAngle)
        maxAngle = angle;
    dp.x = dp.x - 1.0;

    angle = getAngle(vertical, dp - centerPoint);
    if (angle < minAngle)
        minAngle = angle;
    if (angle > maxAngle)
        maxAngle = angle;
}



// OBJECT METHODS
// Finding rotational symmetries.
std::vector<BruteRotationalSymmetry> BruteRotationalSymmetry::getRotationalSymmetries(
    std::vector<std::vector<std::vector<Voxel>>> voxels,
    VoxelMesh& voxelMesh)
{
    std::vector<BruteRotationalSymmetry> symmetries;      // Vector of all symmetries.
    int maxRadius;                              // Maximum radius for searching rotational symmetries.
    int step = 10 * 1;                               // Step in voxel mesh.
    std::vector<BruteRotationalSymmetry> tempSymmetries;  // Vector of symmetries on the radius.

    // Marching the voxel mesh of all layers and finding symmetries.
    for (size_t layer = voxels.size() - 1; layer >= 0; layer--) {
        for (uint i = 15; i < 10 * (voxelMesh.xSize - 1); i += step) {
            for (uint j = 15; j < 10 * (voxelMesh.ySize - 1); j += step) {
                float x = i / 10.0;
                float y = j / 10.0;

                // Setting the coordinates of a new center
                // point, where rotational symmetries are searched.
                Point centerPoint(x, y, layer);

                // Maximum radius where symmetries are searched for, equals
                // the distance to the nearest edge of the bounding box.
                maxRadius = (int)std::min(
                    std::min(voxelMesh.xSize - x, voxelMesh.ySize - y),
                    std::min(x, y)
                );

                // Searching rotational symmetries around the current point.
                tempSymmetries = getRotationalSymmetriesAroundCenterPoint(
                    voxels[layer],
                    centerPoint,
                    maxRadius,
                    voxelMesh.xSize,
                    voxelMesh.ySize
                );

                // Inserting the found symmetries into the vector of all symmetries.
                symmetries.insert(
                    std::end(symmetries),
                    std::begin(tempSymmetries),
                    std::end(tempSymmetries)
                );
            }
        }
    }

    // Merging symmetries according to the axis and the rotation.
    std::vector<BruteRotationalSymmetry> mergedSymmetries = mergeSymmetries(symmetries);

    return mergedSymmetries;
}

// Calculating the point positions according to the axis.
std::vector<PointPosition> BruteRotationalSymmetry::calculatePositionsFromAxis(const Point* axis) {
    const unsigned long long numberOfPoints = m_Points.size();  // Getting the point vector size.

    // If there is no plane, all points are undefined.
    if (axis == nullptr) {
        return std::vector<PointPosition>(numberOfPoints, PointPosition::notInSymmetry);
    }

    // Creating a new vector with positions.
    std::vector<PointPosition> positions(numberOfPoints, PointPosition::notInSymmetry);
    for (unsigned long long i = 0; i < numberOfPoints; i++) {
        // Getting each single point and its voxel coordinates.
        const Point& point = m_Points[i];
        const int voxelX = floor((point.x - m_VoxelMesh.minX) / m_VoxelMesh.sideX);
        const int voxelY = floor((point.y - m_VoxelMesh.minY) / m_VoxelMesh.sideY);
        const int voxelZ = floor((point.z - m_VoxelMesh.minZ) / m_VoxelMesh.sideZ);

        // If the current point lies outside of an interesting voxel, its position is set to undefined.
        if (std::find(voxels.begin(), voxels.end(), Voxel(voxelX, voxelY, voxelZ)) == voxels.end()) {
             positions[i] = PointPosition::notInSymmetry;
             continue;
        }

        // Otherwise, the position is set to rotational.
        positions[i] = PointPosition::rotational;
    }

    return positions;
}
