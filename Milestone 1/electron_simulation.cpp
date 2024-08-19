#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <random>
#include "H5Cpp.h"
using namespace H5;
#include <chrono>

const H5std_string FILE_NAME("electron_positions.h5");
const H5std_string DATASET_NAME("positions");

int randomPos() {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(-10, 10);

    return dist6(rng);
}

int main() {
    int numElectrons = 100;
    int numFrames = 1000;
    
    // creating start positions  
    float** electronPos = (float**)malloc(sizeof(float*) * numElectrons);
    float* electronVel = (float*)malloc(sizeof(float) * numElectrons);
    std::cout << "Initialising electron positions" << std::endl;
    for (int i = 0; i < numElectrons; ++i) {
        electronPos[i] = (float*)malloc(sizeof(float) * 3);
        electronPos[i][0] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        electronPos[i][1] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        electronPos[i][2] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        electronVel[i] = 0;
    }

    // Create HDF5 file thank you chat gpt
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    hsize_t dims[3] = {static_cast<hsize_t>(numFrames), static_cast<hsize_t>(numElectrons), 3};
    DataSpace dataspace(3, dims);
    DataSet dataset = file.createDataSet(DATASET_NAME, PredType::NATIVE_FLOAT, dataspace);

    float* frameData = new float[numElectrons * 3];

    std::cout << "Begin calculating forces" << std::endl;
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for (int f = 0; f < numFrames; ++f) {
        if (electronPos == nullptr || electronVel == nullptr) {
            std::cerr << "Memory allocation failed!" << std::endl;
            return -1;
        }
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        for (int i = 0; i < numElectrons; ++i) {
            if (electronPos[i] == nullptr) {
                std::cerr << "Memory allocation failed!" << std::endl;
                return -1;
            }
            float xComponent = 0;
            float yComponent = 0;
            float zComponent = 0;
            for (int j = 0; j < numElectrons; ++ j) {
                // std::cout << i << " " << j << std::endl;
                if (i == j) {
                    continue;
                }
                const float xDiff = electronPos[j][0] - electronPos[i][0];
                const float yDiff = electronPos[j][1] - electronPos[i][1];
                const float zDiff = electronPos[j][2] - electronPos[i][2];

                const float distance = sqrt(xDiff * xDiff + yDiff * yDiff + zDiff * zDiff);
                // copilot used to fix my bad math
                // and liv for figuring out that my force is 0ed
                //pow(10,-31)
                const float m = 9.1093837 * pow(10,-31);
                const float k = 8.987 * pow(10, 9);
                const float e = 1.602 * pow(10, -19);
                const float t = 0.01;
                float forceMag = k * (pow(e, 2) / pow(distance == 0 ? std::numeric_limits<float>::infinity() : (distance), 2));

                const float angle = atan2(yDiff, xDiff);
                const float angleZ = atan2(zDiff, sqrt(xDiff * xDiff + yDiff * yDiff));
                const float acceleration = forceMag / m;
                // std::cout << "accel: " << acceleration << std::endl;
                const float v = electronVel[i] + acceleration * t;
                const float s = electronVel[i] * t + 0.5 * acceleration * pow(t, 2);
                // std::cout << "s: " << s << std::endl;
                // std::cout << "v: " << v << std::endl;
                xComponent += -cos(angle) * s;
                yComponent += -sin(angle) * s;
                zComponent += -sin(angleZ) * s;
                // std::cout << "z " << zComponent << std::endl;
                //std::cout << "Force magnitude: " << forceMag << std::endl;
                //std::cout << "X Component: " << xComponent << " Y Component: " << yComponent << std::endl;
            }
            // std::cout << "Force magnitude: " << forceMag << std::endl;
            // std::cout << "X Component: " << xComponent << " Y Component: " << yComponent << std::endl;
            electronPos[i][0] = std::max((static_cast<float>(-10)), std::min(electronPos[i][0] + xComponent, (static_cast<float>(10))));
            electronPos[i][1] = std::max((static_cast<float>(-10)), std::min(electronPos[i][1] + yComponent, (static_cast<float>(10))));
            electronPos[i][2] = std::max((static_cast<float>(-10)), std::min(electronPos[i][2] + zComponent, (static_cast<float>(10))));
            // electronPos[i][0] += xComponent;
            // electronPos[i][1] += yComponent;
            // electronPos[i][2] += zComponent;
            // std::cout << "Set X Component: " << electronPos[i][0] << " Set Y Component: " << electronPos[i][1] << std::endl;
            frameData[i * 3] = electronPos[i][0];
            frameData[i * 3 + 1] = electronPos[i][1];
            frameData[i * 3 + 2] = electronPos[i][2];
        }
        hsize_t start[3] = {static_cast<hsize_t>(f), 0, 0};
        hsize_t count[3] = {1, static_cast<hsize_t>(numElectrons), 3};
        dataspace.selectHyperslab(H5S_SELECT_SET, count, start);

        DataSpace memspace(3, count);

        dataset.write(frameData, PredType::NATIVE_FLOAT, memspace, dataspace);
        high_resolution_clock::time_point t4 = high_resolution_clock::now();
        duration<double> time_span_34 = duration_cast<duration<double>>(t4 - t3);
        std::cout << std::to_string(time_span_34.count()) << std::endl;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << std::to_string(time_span.count()) << std::endl;
    // Clean up
    delete[] frameData;
    for (int i = 0; i < numElectrons; ++i) {
        free(electronPos[i]);
    }
    free(electronPos);
    free(electronVel);

    return 0;
}
