#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cstring>
#include <immintrin.h>
#include "H5Cpp.h"
#include "octree.h"
#include <omp.h>

using namespace H5;
using namespace std::chrono;

const int numFrames = 2500;
const int numElectrons = 1000;

std::vector<Vector3> electronPos(numElectrons);
std::vector<Vector3> electronVel(numElectrons);

/*
* Simple function to set position for each electron
*/
void initializeElectrons() {
    #pragma omp parallel for
    for (int i = 0; i < numElectrons; ++i) {
        electronPos[i] = {static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) * 20 - 10,
                          static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) * 20 - 10,
                          static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) * 20 - 10};
        electronVel[i] = {0.0f, 0.0f, 0.0f};
    }
}

/*
* Main loop logic, populate octree with data and calculate force
*/
void updateElectrons() {
    Octree octree({0.0f, 0.0f, 0.0f}, 10.0f);

    // Add all electrons to octree
    for (int i = 0; i < numElectrons; ++i) {
        octree.insert(electronPos[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < numElectrons; ++i) {
        Vector3 force = {0.0f, 0.0f, 0.0f};
        std::vector<Vector3> nearbyElectrons;
        octree.query(electronPos[i], 5.0f, nearbyElectrons);

        for (const auto& pos : nearbyElectrons) {
            // Don't consider the same electron
            if (pos.x == electronPos[i].x && pos.y == electronPos[i].y && pos.z == electronPos[i].z) continue;

            // Calculate the distance
            float dx = pos.x - electronPos[i].x;
            float dy = pos.y - electronPos[i].y;
            float dz = pos.z - electronPos[i].z;
            float distance = sqrt(dx * dx + dy * dy + dz * dz);
            // Calculate the force
            float forceMagnitude = -(k * (e * e) / (distance * distance));
            force.x += forceMagnitude * dx / distance;
            force.y += forceMagnitude * dy / distance;
            force.z += forceMagnitude * dz / distance;
        }
        // Apply force
        electronVel[i].x += (force.x / m) * t;
        electronVel[i].y += (force.y / m) * t;
        electronVel[i].z += (force.z / m) * t;

        electronPos[i].x += electronVel[i].x * t;
        electronPos[i].y += electronVel[i].y * t;
        electronPos[i].z += electronVel[i].z * t;

        electronPos[i].x = std::max(-10.0f, std::min(electronPos[i].x, 10.0f));
        electronPos[i].y = std::max(-10.0f, std::min(electronPos[i].y, 10.0f));
        electronPos[i].z = std::max(-10.0f, std::min(electronPos[i].z, 10.0f));

    }
}

int main(int argc, char* argv[]) {
    // Initialise variables
    bool timingMode = false;
    const H5std_string FILE_NAME("electron_positions.h5");
    const H5std_string DATASET_NAME("positions");

    // Parse command-line arguments
    // Accuracy mode was removed.
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-timing") == 0) {
            timingMode = true;
        }
    }

    // Randomise and initialise the electrons
    std::srand(std::time(0));
    initializeElectrons();

    // Handle HDF5 file.
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

    // Main loop
    for (int f = 0; f < numFrames; ++f) {
        high_resolution_clock::time_point t3;
        if (timingMode) {
            t3 = high_resolution_clock::now();
        }
        updateElectrons();

        #pragma omp parallel for
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
        if (timingMode) {
            high_resolution_clock::time_point t4 = high_resolution_clock::now();
            duration<double> time_span_34 = duration_cast<duration<double>>(t4 - t3);
            std::cout << std::to_string(time_span_34.count()) << std::endl;
        }
    }

    delete[] frameData;

    if (timingMode) {
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
        std::cout << std::to_string(time_span.count()) << std::endl;
    }

    return 0;
}