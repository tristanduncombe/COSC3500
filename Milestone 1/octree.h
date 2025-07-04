#include <vector>
#include <memory>
#include <array>
#include <cmath>
#include <algorithm>
#include "simulation.h"

/*
* Octree class, AI was used in the development of this file.
*/
class Octree {

public:
    // Constructor, set the empty tree
    Octree(const Vector3& center, float halfSize)
        : center(center), halfSize(halfSize), isLeaf(true) {
        for (auto& child : children) {
            child = nullptr;
        }
    }
    // Insert into the octree
    void insert(const Vector3& point) {
        // If we are at leaf node we need to increase the size of the tree
        if (isLeaf) {
            if (points.size() < maxPoints || halfSize < minSize) {
                points.push_back(point);
            } else {
                subdivide();
                for (const auto& p : points) {
                    insertIntoChildren(p);
                }
                points.clear();
                insertIntoChildren(point);
            }
        } else {
            insertIntoChildren(point);
        }
    }
    // Find elements within a range around a point
    void query(const Vector3& point, float range, std::vector<Vector3>& result) const {
        if (!intersects(point, range)) return;

        if (isLeaf) {
            for (const auto& p : points) {
                if (distanceSquared(p, point) <= range * range) {
                    result.push_back(p);
                }
            }
        } else {
            for (const auto& child : children) {
                if (child) {
                    child->query(point, range, result);
                }
            }
        }
    }
// Methods used for the tree that should not be public
private:
    // variables used to defined the tree
    static constexpr int maxPoints = 4;
    static constexpr float minSize = 1.0f;

    Vector3 center;
    float halfSize;
    bool isLeaf;
    std::vector<Vector3> points;
    std::array<std::unique_ptr<Octree>, 8> children;

    // Divide the current node into 8 child nodes.
    void subdivide() {
        isLeaf = false;
        float quarterSize = halfSize / 2.0f;
        for (int i = 0; i < 8; ++i) {
            Vector3 newCenter = center;
            newCenter.x += (i & 1 ? quarterSize : -quarterSize);
            newCenter.y += (i & 2 ? quarterSize : -quarterSize);
            newCenter.z += (i & 4 ? quarterSize : -quarterSize);
            children[i] = std::unique_ptr<Octree>(new Octree(newCenter, quarterSize));
        }
    }

    // Inserts a vector point into the right child node
    void insertIntoChildren(const Vector3& point) {
        for (const auto& child : children) {
            if (child && child->contains(point)) {
                child->insert(point);
                return;
            }
        }
    }

    // Check if the node contains a specific point
    bool contains(const Vector3& point) const {
        return point.x >= center.x - halfSize && point.x <= center.x + halfSize &&
               point.y >= center.y - halfSize && point.y <= center.y + halfSize &&
               point.z >= center.z - halfSize && point.z <= center.z + halfSize;
    }

    // Check if the current point intersects with a given point and range.
    bool intersects(const Vector3& point, float range) const {
        float dx = std::max(center.x - halfSize - point.x, 0.0f) + std::max(point.x - center.x - halfSize, 0.0f);
        float dy = std::max(center.y - halfSize - point.y, 0.0f) + std::max(point.y - center.y - halfSize, 0.0f);
        float dz = std::max(center.z - halfSize - point.z, 0.0f) + std::max(point.z - center.z - halfSize, 0.0f);
        return dx * dx + dy * dy + dz * dz <= range * range;
    }

    // Get the distance squared
    float distanceSquared(const Vector3& a, const Vector3& b) const {
        float dx = a.x - b.x;
        float dy = a.y - b.y;
        float dz = a.z - b.z;
        return dx * dx + dy * dy + dz * dz;
    }
};