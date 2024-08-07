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
    std::uniform_int_distribution<std::mt19937::result_type> dist6(1,6);

    return dist6(rng);
}
 

int main() {

    // assemble grid
    // add initial atom state
    int numElectrons = 100;

    float** electronPos = (float**) malloc(sizeof (float**) * 100);

    // creating start positions  
    for (int i = 0; i < numElectrons; ++i) {
        electronPos[i] = (float *) malloc(sizeof (float*) * 2);
        electronPos[i][0] = randomPos();
        electronPos[i][1] = randomPos();
    }

    for (int i = 0; i < numElectrons; ++i) {
        std::cout << "(";
        std::cout << std::setprecision(15) << electronPos[i][0];
        std::cout << ",";
        std::cout << std::setprecision(15) << electronPos[i][1];
        std::cout << ")"<<std::endl;
    }
    std::cout << "moving" << std::endl;

    for (int f = 0; f < 100; ++f) {

        // calculate each probability of each atom
        // if thingy change
        // calculate energy???????
        for (int i = 0; i < numElectrons; ++i) {
            float xComponent = 0;
            float yComponent = 0;
            for (int j = 0; j < numElectrons; ++ j) {
                // std::cout << i << " " << j << std::endl;
                if (i == j) {
                    continue;
                }
                const float xDiff = electronPos[j][0] - electronPos[i][0];
                const float yDiff = electronPos[j][1] - electronPos[i][1];

                const float distance = hypot(xDiff, yDiff);
                // copilot used to fix my bad math
                // and liv for figuring out that my force is 0ed

                const float k = 8.987 * pow(10, 9);
                const float e = 1.602 * pow(10, -19);
                float forceMag = k * (pow(e, 2) / pow(distance == 0 ? std::numeric_limits<float>::infinity() : distance, 2));

                const float angle = atan2(yDiff, xDiff);
                xComponent += cos(angle) * forceMag;
                yComponent += sin(angle) * forceMag;
                //std::cout << "Force magnitude: " << forceMag << std::endl;
                //std::cout << "X Component: " << xComponent << " Y Component: " << yComponent << std::endl;
            }
            // std::cout << "Force magnitude: " << forceMag << std::endl;
            std::cout << "X Component: " << xComponent << " Y Component: " << yComponent << std::endl;

            electronPos[i][0] += xComponent;
            electronPos[i][1] += yComponent;
            std::cout << "Set X Component: " << electronPos[i][0] << " Set Y Component: " << electronPos[i][1] << std::endl;
        } 

    }

    for (int i = 0; i < numElectrons; ++i) {
        std::cout << "(";
        std::cout << std::setprecision(15) << electronPos[i][0];
        std::cout << ",";
        std::cout << std::setprecision(15) << electronPos[i][1];
        std::cout << ")"<<std::endl;
    }
    
    // return it
    return 0;
}