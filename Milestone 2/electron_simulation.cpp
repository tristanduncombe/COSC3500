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
#include <mpi.h>
#include "H5Cpp.h"
#include "octree.h"
#include <omp.h>

using namespace H5;
using namespace std::chrono;

const int numFrames = 2500;
const int numElectrons = 5000;

// Function to initialize electrons
void initializeElectrons(int numLocalElectrons, Vector3* electronPosLocal, Vector3* electronVelLocal) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand(time(NULL) + rank);
    #pragma omp parallel for
    for (int i = 0; i < numLocalElectrons; ++i) {
        electronPosLocal[i].x = static_cast<float>(rand()) / RAND_MAX * 20 - 10;
        electronPosLocal[i].y = static_cast<float>(rand()) / RAND_MAX * 20 - 10;
        electronPosLocal[i].z = static_cast<float>(rand()) / RAND_MAX * 20 - 10;
        electronVelLocal[i] = {0.0f, 0.0f, 0.0f};
    }
}

// Helper function for horizontal addition
inline float _mm256_reduce_add_ps(__m256 x) {
    __m128 hi = _mm256_extractf128_ps(x, 1);
    __m128 lo = _mm256_castps256_ps128(x);
    lo = _mm_add_ps(lo, hi);
    __m128 shuf = _mm_movehdup_ps(lo);
    __m128 sums = _mm_add_ps(lo, shuf);
    shuf = _mm_movehl_ps(shuf, sums);
    sums = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
}

// Function to update electrons
void updateElectrons(int numLocalElectrons, int numTotalElectrons, Vector3* electronPosLocal, Vector3* electronVelLocal, Vector3* allElectronPos) {
    Octree octree({0.0f, 0.0f, 0.0f}, 20.0f);

    // Build octree with all electrons
    for (int i = 0; i < numTotalElectrons; ++i) {
        octree.insert(allElectronPos[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < numLocalElectrons; ++i) {
        Vector3 force = {0.0f, 0.0f, 0.0f};
        std::vector<Vector3> nearbyElectrons;
        octree.query(electronPosLocal[i], 5.0f, nearbyElectrons);

        int numNearby = nearbyElectrons.size();
        int j = 0;

        // Constants for AVX calculations
        __m256 posIX = _mm256_set1_ps(electronPosLocal[i].x);
        __m256 posIY = _mm256_set1_ps(electronPosLocal[i].y);
        __m256 posIZ = _mm256_set1_ps(electronPosLocal[i].z);
        __m256 k_e2 = _mm256_set1_ps(-(k * (e * e)));
        __m256 one = _mm256_set1_ps(1.0f);
        __m256 zero = _mm256_setzero_ps();

        for (; j <= numNearby - 8; j += 8) {
            // Load positions of nearby electrons
            float posJX_f[8], posJY_f[8], posJZ_f[8];
            for (int n = 0; n < 8; ++n) {
                posJX_f[n] = nearbyElectrons[j + n].x;
                posJY_f[n] = nearbyElectrons[j + n].y;
                posJZ_f[n] = nearbyElectrons[j + n].z;
            }
            __m256 posJX = _mm256_loadu_ps(posJX_f);
            __m256 posJY = _mm256_loadu_ps(posJY_f);
            __m256 posJZ = _mm256_loadu_ps(posJZ_f);

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
            force.x += _mm256_reduce_add_ps(forceX);
            force.y += _mm256_reduce_add_ps(forceY);
            force.z += _mm256_reduce_add_ps(forceZ);
        }

        // Handle remaining electrons
        for (; j < numNearby; ++j) {
            float dx = nearbyElectrons[j].x - electronPosLocal[i].x;
            float dy = nearbyElectrons[j].y - electronPosLocal[i].y;
            float dz = nearbyElectrons[j].z - electronPosLocal[i].z;
            float distSqr = dx * dx + dy * dy + dz * dz;

            if (distSqr > 0.0f) {
                float forceMag = -(k * (e * e)) / distSqr;
                force.x += forceMag * dx;
                force.y += forceMag * dy;
                force.z += forceMag * dz;
            }
        }
        // Apply force
        electronVelLocal[i].x += (force.x / m) * t;
        electronVelLocal[i].y += (force.y / m) * t;
        electronVelLocal[i].z += (force.z / m) * t;

        electronPosLocal[i].x += electronVelLocal[i].x * t;
        electronPosLocal[i].y += electronVelLocal[i].y * t;
        electronPosLocal[i].z += electronVelLocal[i].z * t;

        electronPosLocal[i].x = std::max(-10.0f, std::min(electronPosLocal[i].x, 10.0f));
        electronPosLocal[i].y = std::max(-10.0f, std::min(electronPosLocal[i].y, 10.0f));
        electronPosLocal[i].z = std::max(-10.0f, std::min(electronPosLocal[i].z, 10.0f));
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    bool timingMode = false;

    const H5std_string FILE_NAME("electron_positions.h5");
    const H5std_string DATASET_NAME("positions");

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-timing") == 0) {
            timingMode = true;
        }
    }

    // Determine the number of electrons per process
    int baseNumLocalElectrons = numElectrons / size;
    int remainder = numElectrons % size;
    int numLocalElectrons = baseNumLocalElectrons + (rank < remainder ? 1 : 0);

    // Allocate local arrays
    Vector3* electronPosLocal = new Vector3[numLocalElectrons];
    Vector3* electronVelLocal = new Vector3[numLocalElectrons];

    // Initialize local electrons
    initializeElectrons(numLocalElectrons, electronPosLocal, electronVelLocal);

    // Allocate array to gather all electron positions
    Vector3* allElectronPos = new Vector3[numElectrons];

    // File handling (only by root process)
    H5File* file = nullptr;
    DataSet* dataset = nullptr;
    DataSpace* dataspace = nullptr;
    float* frameData = nullptr;

    if (rank == 0) {
        file = new H5File(FILE_NAME, H5F_ACC_TRUNC);
        hsize_t dims[3] = {static_cast<hsize_t>(numFrames), static_cast<hsize_t>(numElectrons), 3};
        dataspace = new DataSpace(3, dims);
        dataset = new DataSet(file->createDataSet(DATASET_NAME, PredType::NATIVE_FLOAT, *dataspace));
        frameData = new float[numElectrons * 3];
    }

    high_resolution_clock::time_point t1;
    if (rank == 0 && timingMode) {
        std::cout << "Begin calculating forces" << std::endl;
        t1 = high_resolution_clock::now();
    }

    // Main loop
    for (int f = 0; f < numFrames; ++f) {
        high_resolution_clock::time_point t3;
    if (timingMode) {
        t3 = high_resolution_clock::now();
    }
        int* recvCounts = new int[size];
        int* displs = new int[size];
        for (int i = 0; i < size; ++i) {
            int count = baseNumLocalElectrons + (i < remainder ? 1 : 0);
            recvCounts[i] = count * sizeof(Vector3);
            displs[i] = (baseNumLocalElectrons * i + std::min(i, remainder)) * sizeof(Vector3);
        }

        MPI_Allgatherv(electronPosLocal, numLocalElectrons * sizeof(Vector3), MPI_BYTE,
                       allElectronPos, recvCounts, displs, MPI_BYTE,
                       MPI_COMM_WORLD);

        // Update electrons
        updateElectrons(numLocalElectrons, numElectrons, electronPosLocal, electronVelLocal, allElectronPos);

        // Optionally write data to file (only by root process)
        if (rank == 0) {
            // Gather all electron positions from all processes
            MPI_Gatherv(MPI_IN_PLACE, numLocalElectrons * sizeof(Vector3), MPI_BYTE,
                        allElectronPos, recvCounts, displs, MPI_BYTE,
                        0, MPI_COMM_WORLD);

            // Prepare data for HDF5
            #pragma omp parallel for
            for (int i = 0; i < numElectrons; ++i) {
                frameData[i * 3] = allElectronPos[i].x;
                frameData[i * 3 + 1] = allElectronPos[i].y;
                frameData[i * 3 + 2] = allElectronPos[i].z;
            }

            // Write frame data to HDF5 file
            hsize_t offset[3] = {static_cast<hsize_t>(f), 0, 0};
            hsize_t count[3] = {1, static_cast<hsize_t>(numElectrons), 3};
            DataSpace memspace(3, count);
            dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
            dataset->write(frameData, PredType::NATIVE_FLOAT, memspace, *dataspace);

            high_resolution_clock::time_point loop_end = high_resolution_clock::now();
            if (rank == 0 && timingMode) {
                high_resolution_clock::time_point t4 = high_resolution_clock::now();
                duration<double> time_span_34 = duration_cast<duration<double>>(t4 - t3);
                std::cout << std::to_string(time_span_34.count()) << std::endl;
            }
        } else {
            // Non-root processes participate in MPI_Gatherv
            MPI_Gatherv(electronPosLocal, numLocalElectrons * sizeof(Vector3), MPI_BYTE,
                        nullptr, nullptr, nullptr, MPI_BYTE,
                        0, MPI_COMM_WORLD);
        }

        delete[] recvCounts;
        delete[] displs;
    }

    // Clean up
    delete[] electronPosLocal;
    delete[] electronVelLocal;
    delete[] allElectronPos;

    if (rank == 0) {
        delete[] frameData;
        delete dataset;
        delete dataspace;
        delete file;
    }

    if (rank == 0 && timingMode) {
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
        std::cout << std::to_string(time_span.count()) << std::endl;
    }

    MPI_Finalize();
    return 0;
}