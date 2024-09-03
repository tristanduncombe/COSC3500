#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <vector> // Include the vector header
#include <cstring>
#include "H5Cpp.h"

using namespace H5;
using namespace std::chrono;

const H5std_string FILE_NAME("electron_positions.h5");
const H5std_string DATASET_NAME("positions");

const float k = 8.99e9; // Coulomb's constant
const float e = 1.6e-19; // Elementary charge
const float m = 9.11e-31; // Electron mass
const float t = 1e-4; // Time step
const int numFrames = 2500;
const int numElectrons = 1000;
const int gridDimX = 10;
const int gridDimY = 10;
const int gridDimZ = 10;
const int numGrids = gridDimX * gridDimY * gridDimZ;
const float gridSize = 1.0;

struct Vector3 {
    float x, y, z;
};

std::vector<Vector3> electronPos(numElectrons);
std::vector<Vector3> electronVel(numElectrons);

int getCubeIndex(float x, float y, float z, float gridSize, int gridDimX, int gridDimY, int gridDimZ) {
    // Shift the range of x, y, z from [-10, 10] to [0, 20]
    x += 10;
    y += 10;
    z += 10;
    int ix = static_cast<int>(x / gridSize);
    int iy = static_cast<int>(y / gridSize);
    int iz = static_cast<int>(z / gridSize);
    return ix + iy * gridDimX + iz * gridDimX * gridDimY;
}

void initializeElectrons() {
    for (int i = 0; i < numElectrons; ++i) {
        electronPos[i] = {static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) * 20 - 10,
                          static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) * 20 - 10,
                          static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) * 20 - 10};
        electronVel[i] = {0.0f, 0.0f, 0.0f};
    }
}

void updateElectrons() {
    for (int i = 0; i < numElectrons; ++i) {
        Vector3 force = {0.0f, 0.0f, 0.0f};

        for (int j = 0; j < numElectrons; ++j) {
            if (i == j) continue;

            float dx = electronPos[j].x - electronPos[i].x;
            float dy = electronPos[j].y - electronPos[i].y;
            float dz = electronPos[j].z - electronPos[i].z;
            float distance = sqrt(dx * dx + dy * dy + dz * dz);

            float forceMagnitude = -(k * (e * e) / (distance * distance));
            force.x += forceMagnitude * dx / distance;
            force.y += forceMagnitude * dy / distance;
            force.z += forceMagnitude * dz / distance;
        }
        // Update velocity: v = v + (F/m) * t
        electronVel[i].x += (force.x / m) * t;
        electronVel[i].y += (force.y / m) * t;
        electronVel[i].z += (force.z / m) * t;

        // Update position: x = x + v * t
        electronPos[i].x += electronVel[i].x * t;
        electronPos[i].y += electronVel[i].y * t;
        electronPos[i].z += electronVel[i].z * t;

        // Ensure electrons stay within bounds
        electronPos[i].x = std::max(-10.0f, std::min(electronPos[i].x, 10.0f));
        electronPos[i].y = std::max(-10.0f, std::min(electronPos[i].y, 10.0f));
        electronPos[i].z = std::max(-10.0f, std::min(electronPos[i].z, 10.0f));
    }
}

int main(int argc, char* argv[]) {
    bool timingMode = false;

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-timing") == 0) {
            timingMode = true;
        }
    }

    std::srand(std::time(0));
    initializeElectrons();

    // Create HDF5 file
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    hsize_t dims[3] = {static_cast<hsize_t>(numFrames), static_cast<hsize_t>(numElectrons), 3};
    DataSpace dataspace(3, dims);
    DataSet dataset = file.createDataSet(DATASET_NAME, PredType::NATIVE_FLOAT, dataspace);

    float* frameData = new float[numElectrons * 3];

    high_resolution_clock::time_point t1;
    if (timingMode) {
        std::cout << "Begin calculating forces" << std::endl;
        t1 = high_resolution_clock::now();
    }
    for (int f = 0; f < numFrames; ++f) {
        high_resolution_clock::time_point t3;
        if (timingMode) {
            t3 = high_resolution_clock::now();
        }
        updateElectrons();

        for (int i = 0; i < numElectrons; ++i) {
            frameData[i * 3] = electronPos[i].x;
            frameData[i * 3 + 1] = electronPos[i].y;
            frameData[i * 3 + 2] = electronPos[i].z;
        }

        // Write frame data to HDF5 file
        hsize_t offset[3] = {static_cast<hsize_t>(f), 0, 0};
        hsize_t count[3] = {1, static_cast<hsize_t>(numElectrons), 3};
        DataSpace memspace(3, count);
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
        dataset.write(frameData, PredType::NATIVE_FLOAT, memspace, dataspace);

        high_resolution_clock::time_point loop_end = high_resolution_clock::now();
        high_resolution_clock::time_point t4;
        duration<double> time_span_34;
        if (timingMode) {
            high_resolution_clock::time_point t4 = high_resolution_clock::now();
            duration<double> time_span_34 = duration_cast<duration<double>>(t4 - t3);
            std::cout << std::to_string(time_span_34.count()) << std::endl;
        }
    }

    delete[] frameData;

    high_resolution_clock::time_point t2;
    duration<double> time_span;
    if (timingMode) {
        t2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
        std::cout << std::to_string(time_span.count()) << std::endl;
    }

    return 0;
}