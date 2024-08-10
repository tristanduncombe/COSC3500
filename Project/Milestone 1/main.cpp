#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <random>

int randomPos() {
    // https://stackoverflow.com/questions/13445688/how-to-generate-a-random-number-in-c

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(-10,10);

    return dist6(rng);
}

int main() {

    // assemble grid
    // add initial atom state
    int numElectrons = 100;
    srand (static_cast <unsigned> (time(0)));

    float** electronPos = (float**) malloc(sizeof (float**) * 100);
    float* electronVel = (float*) malloc(sizeof (float*) * 100);
    std::cout << "[" << std::endl; 
    // creating start positions  
    for (int i = 0; i < numElectrons; ++i) {
        electronPos[i] = (float *) malloc(sizeof (float*) * 3);
        electronPos[i][0] = static_cast<float>(randomPos());
        electronPos[i][1] = static_cast<float>(randomPos());
        electronPos[i][2] = static_cast<float>(randomPos());
        electronVel[i] = 0;
    }

    for (int i = 0; i < numElectrons; ++i) {
            std::cout << "\t";
            std::cout << "[";
            std::cout << std::setprecision(64) << electronPos[i][0];
            std::cout << ", ";
            std::cout << std::setprecision(64) << electronPos[i][1];
            std::cout << ", ";
            std::cout << std::setprecision(64) << electronPos[i][2];
            if (numElectrons - i == 1) {
                std::cout << "]" <<std::endl;
            }
            else {
                std::cout << "],"<<std::endl;
            }
    }
    std::cout << "]," << std::endl;
    int numFrames = 1000;
    for (int f = 0; f < numFrames; ++f) {

        // calculate each probability of each atom
        // if thingy change
        // calculate energy???????
        for (int i = 0; i < numElectrons; ++i) {
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
                const float s = electronVel[i] * 0.01 + 0.5 * acceleration * pow(t, 2);
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
        } 
        std::cout << "[" << std::endl;
        for (int i = 0; i < numElectrons; ++i) {
            std::cout << "\t";
            std::cout << "[";
            std::cout << std::setprecision(64) << electronPos[i][0];
            std::cout << ", ";
            std::cout << std::setprecision(64) << electronPos[i][1];
            std::cout << ", ";
            std::cout << std::setprecision(64) << electronPos[i][2];
            if (numElectrons - i == 1) {
                std::cout << "]" <<std::endl;
            }
            else {
                std::cout << "],"<<std::endl;
            }
        }
        if (numFrames - f == 1) {
            std::cout << "]" << std::endl;
        }
        else {
            std::cout << "]," << std::endl;
        }
    }
    
    // return it
    return 0;
}