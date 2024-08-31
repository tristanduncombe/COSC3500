#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include "H5Cpp.h"

using namespace H5;
using namespace std::chrono;

const H5std_string FILE_NAME("electron_positions.h5");
const H5std_string DATASET_NAME("positions");

const float k = 8.9875517873681764e9; // Coulomb's constant
const float e = 1.602176634e-19; // Elementary charge
const float m = 9.10938356e-31; // Electron mass
const float t = 1e-2; // Time step

int getCubeIndex(float x, float y, float z, float gridSize) {
    int cubeX = static_cast<int>(x / gridSize);
    int cubeY = static_cast<int>(y / gridSize);
    int cubeZ = static_cast<int>(z / gridSize);

    // Calculate the index of the cube in the gridForce array
    int cubeIndex = cubeX + cubeY * 3 + cubeZ * 9;

    return cubeIndex;
}


int main(int argc, char* argv[]) {
    int numElectrons = 1000;
    int numFrames = 2500;
    bool timing = false;
    bool accuracy = false;
    bool random = true;

    // Define grid dimensions
    const int gridDimX = 3;
    const int gridDimY = 3;
    const int gridDimZ = 3;
    const int numGrids = gridDimX * gridDimY * gridDimZ;
    const float gridSize = 20.0 / gridDimX;

    float* gridForce = (float*)malloc(sizeof(float) * numGrids);

    // Initialize gridForce to zero
    for (int i = 0; i < numGrids; ++i) {
        gridForce[i] = 0.0;
    }

    // Check for -timing flag
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-timing") {
            timing = true;
        }
        if (std::string(argv[i]) == "-accuracy") {
            accuracy = true;
        }
        if (std::string(argv[i]) == "-non-random") {
            random = false;
        }
    }
    
    // Initialize random seed
    if (random) {
        std::srand(std::time(0));
    }

    // Creating start positions  
    float** electronPos = (float**)malloc(sizeof(float*) * numElectrons);
    float* electronInaccuracy;
    if (accuracy) {
        electronInaccuracy = (float*)malloc(sizeof(float) * numElectrons);
    }
    float* electronVel = (float*)malloc(sizeof(float) * numElectrons);
    for (int i = 0; i < numElectrons; ++i) {
        electronPos[i] = (float*)malloc(sizeof(float) * 3);
        electronPos[i][0] = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
        electronPos[i][1] = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
        electronPos[i][2] = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
        electronVel[i] = 0;

        gridForce[getCubeIndex(electronPos[i][0], electronPos[i][1], electronPos[i][2], gridSize)] += e;
    }

    // Create HDF5 file thank you chat gpt
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    hsize_t dims[3] = {static_cast<hsize_t>(numFrames), static_cast<hsize_t>(numElectrons), 3};
    DataSpace dataspace(3, dims);
    DataSet dataset = file.createDataSet(DATASET_NAME, PredType::NATIVE_FLOAT, dataspace);

    float* frameData = new float[numElectrons * 3];

    high_resolution_clock::time_point t1;
    if (timing) {
        std::cout << "Begin calculating forces" << std::endl;
        t1 = high_resolution_clock::now();
    }
    for (int f = 0; f < numFrames; ++f) {
        if (electronPos == nullptr || electronVel == nullptr) {
            std::cerr << "Memory allocation failed!" << std::endl;
            return -1;
        }
        high_resolution_clock::time_point t3;
        if (timing) {
            t3 = high_resolution_clock::now();
        }
        for (int i = 0; i < numElectrons; ++i) {
            if (electronPos[i] == nullptr) {
                std::cerr << "Memory allocation failed!" << std::endl;
                return -1;
            }
            float xComponent = 0;
            float yComponent = 0;
            float zComponent = 0; 
            int curGrid = getCubeIndex(electronPos[i][0], electronPos[i][1], electronPos[i][2], gridSize);
            for (int j = 0; j < numGrids; ++j) {
                const float xDiff = (j % gridDimX) * gridSize - electronPos[i][0];
                const float yDiff = ((j / gridDimX) % gridDimY) * gridSize - electronPos[i][1];
                const float zDiff = (j / (gridDimX * gridDimY)) * gridSize - electronPos[i][2];
                // std::cout << xDiff << "," << yDiff << "," << zDiff << std::endl;
                const float distance = sqrt(xDiff * xDiff + yDiff * yDiff + zDiff * zDiff);

                float forceMag = k * ((e * gridForce[j]) / pow(distance == 0 ? std::numeric_limits<float>::infinity() : distance, 2));
                // std::cout << forceMag << std::endl;
                const float angle = atan2(yDiff, xDiff);
                const float angleZ = atan2(zDiff, sqrt(xDiff * xDiff + yDiff * yDiff));
                const float acceleration = forceMag / m;
                const float v = electronVel[i] + acceleration * t;
                const float s = electronVel[i] * t + 0.5 * acceleration * pow(t, 2);
                xComponent += -cos(angle) * s;
                yComponent += -sin(angle) * s;
                zComponent += -sin(angleZ) * s;
            }
            if (xComponent == 0 && yComponent == 0 && zComponent == 0) {
                std::cout << "All 0" << std::endl;
            }
            electronPos[i][0] = std::max((static_cast<float>(-10)), std::min(electronPos[i][0] + xComponent, (static_cast<float>(10))));
            electronPos[i][1] = std::max((static_cast<float>(-10)), std::min(electronPos[i][1] + yComponent, (static_cast<float>(10))));
            electronPos[i][2] = std::max((static_cast<float>(-10)), std::min(electronPos[i][2] + zComponent, (static_cast<float>(10))));
            int newGrid = getCubeIndex(electronPos[i][0], electronPos[i][1], electronPos[i][2], gridSize);
            gridForce[curGrid] -= e;
            gridForce[newGrid] += e;
            frameData[i * 3] = electronPos[i][0];
            frameData[i * 3 + 1] = electronPos[i][1];
            frameData[i * 3 + 2] = electronPos[i][2];
        }
        hsize_t start[3] = {static_cast<hsize_t>(f), 0, 0};
        hsize_t count[3] = {1, static_cast<hsize_t>(numElectrons), 3};
        dataspace.selectHyperslab(H5S_SELECT_SET, count, start);

        DataSpace memspace(3, count);

        dataset.write(frameData, PredType::NATIVE_FLOAT, memspace, dataspace);
        high_resolution_clock::time_point t4;
        duration<double> time_span_34;
        if (timing) {
            t4 = high_resolution_clock::now();
            time_span_34 = duration_cast<duration<double>>(t4 - t3);
            std::cout << std::to_string(time_span_34.count()) << std::endl;
        }
    }

    if (timing) {
        auto t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(t2 - t1).count();
        std::cout << "Total duration: " << duration << " microseconds" << std::endl;
    }

    // Free allocated memory
    for (int i = 0; i < numElectrons; ++i) {
        free(electronPos[i]);
    }
    free(electronPos);
    free(electronVel);
    free(electronInaccuracy);
    delete[] frameData;
    free(gridForce);

    return 0;
}