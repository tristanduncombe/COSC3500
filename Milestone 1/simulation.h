/*
    Common functions for both the simulation and tree. These are 
*/
struct Vector3 {
    float x, y, z;

    bool operator==(const Vector3& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

const float k = 8.99e9;
const float e = 1.6e-19;
const float m = 9.11e-31;
const float t = 1e-4; 