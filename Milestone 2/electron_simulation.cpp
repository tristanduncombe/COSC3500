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
#include <immintrin.h>

using namespace H5;
using namespace std::chrono;

const int numFrames = 2500;
const int numElectrons = 5000;

alignas(32) Vector3 electronPos[numElectrons];
alignas(32) Vector3 electronVel[numElectrons];

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

// Not available in AVX, but neatens our code
inline float _mm256_reduce_add_ps(__m256 x) {
    __m128 lo = _mm256_castps256_ps128(x);
    __m128 hi = _mm256_extractf128_ps(x, 1);
    lo = _mm_add_ps(lo, hi);
    __m128 shuf = _mm_movehdup_ps(lo);
    __m128 sums = _mm_add_ps(lo, shuf);
    shuf = _mm_movehl_ps(shuf, sums);
    sums = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
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

        int numNearby = nearbyElectrons.size();
        int j = 0;

        // Constants for AVX calculations
        __m256 posIX = _mm256_set1_ps(electronPos[i].x);
        __m256 posIY = _mm256_set1_ps(electronPos[i].y);
        __m256 posIZ = _mm256_set1_ps(electronPos[i].z);
        __m256 k_e2 = _mm256_set1_ps(-(k * (e * e)));
        __m256 one = _mm256_set1_ps(1.0f);
        __m256 zero = _mm256_setzero_ps();

        for (; j <= numNearby - 8; j += 8) {
            // Load positions of nearby electrons
            __m256 posJX = _mm256_loadu_ps(reinterpret_cast<const float*>(&nearbyElectrons[j].x));
            __m256 posJY = _mm256_loadu_ps(reinterpret_cast<const float*>(&nearbyElectrons[j].y));
            __m256 posJZ = _mm256_loadu_ps(reinterpret_cast<const float*>(&nearbyElectrons[j].z));

            // Compute distance components
            __m256 dx = _mm256_sub_ps(posJX, posIX);
            __m256 dy = _mm256_sub_ps(posJY, posIY);
            __m256 dz = _mm256_sub_ps(posJZ, posIZ);

            // Compute squared distance
            __m256 dx2 = _mm256_mul_ps(dx, dx);
            __m256 dy2 = _mm256_mul_ps(dy, dy);
            __m256 dz2 = _mm256_mul_ps(dz, dz);
            __m256 distSqr = _mm256_add_ps(_mm256_add_ps(dx2, dy2), dz2);

            __m256 distance = _mm256_sqrt_ps(distSqr);

            // Avoid division by zero
            __m256 mask = _mm256_cmp_ps(distance, zero, _CMP_GT_OQ);
            __m256 invDist = _mm256_div_ps(one, distance);

            // Compute force magnitude
            __m256 forceMag = _mm256_div_ps(k_e2, distSqr);

            // Compute force components
            __m256 forceX = _mm256_mul_ps(forceMag, dx);
            __m256 forceY = _mm256_mul_ps(forceMag, dy);
            __m256 forceZ = _mm256_mul_ps(forceMag, dz);

            // Apply mask to ignore division by zero results
            forceX = _mm256_blendv_ps(zero, forceX, mask);
            forceY = _mm256_blendv_ps(zero, forceY, mask);
            forceZ = _mm256_blendv_ps(zero, forceZ, mask);

            // Sum the force components
            float fx[8], fy[8], fz[8];
            _mm256_storeu_ps(fx, forceX);
            _mm256_storeu_ps(fy, forceY);
            _mm256_storeu_ps(fz, forceZ);

            for (int k = 0; k < 8; ++k) {
                force.x += fx[k];
                force.y += fy[k];
                force.z += fz[k];
            }
        }

        // Handle remaining electrons
        for (; j < numNearby; ++j) {
            float dx = nearbyElectrons[j].x - electronPos[i].x;
            float dy = nearbyElectrons[j].y - electronPos[i].y;
            float dz = nearbyElectrons[j].z - electronPos[i].z;
            float distSqr = dx * dx + dy * dy + dz * dz;

            if (distSqr > 0.0f) {
                float forceMag = -(k * (e * e)) / distSqr;
                force.x += forceMag * dx;
                force.y += forceMag * dy;
                force.z += forceMag * dz;
            }
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
    // std::srand(std::time(0));
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