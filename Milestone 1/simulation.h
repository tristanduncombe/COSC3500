/*
    Common functions for both the simulation and tree. These are simple structures and constants
*/
struct Vector3 {
    float x, y, z;

    bool operator==(const Vector3& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};
// Coloumb's Constant
const float k = 8.99e9;
// Charge of an electron
const float e = 1.6e-19;
// Mass of an electron
const float m = 9.11e-31;
// Time step size
const float t = 1e-4; 