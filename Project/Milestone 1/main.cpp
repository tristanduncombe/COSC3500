#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <random>

// double **grid(int rows, int columns) {
//     double **assembledGrid;
//     assembledGrid = malloc(sizeof (double) * rows);

//     for (int i = 0; i < rows; ++i) {
//         assembledGrid[i] = malloc(sizeof (double) * columns);
//     }

//     for (int i = 0; i < rows; ++i) {
//         for (int j = 0; i < columns; ++i) {
//             assembledGrid[i][j] = 0.;
//         }
//     }

//     return assembledGrid;
// }

int randomPos() {
    // https://stackoverflow.com/questions/13445688/how-to-generate-a-random-number-in-c

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(1,6);

    return dist6(rng);
}
 
int main(int width, int height, int loops) {

    // assemble grid
    // add initial atom state
    int numberAtoms = 100;

    float** electronPos = (float**) malloc(sizeof (float**) * 100);

    // creating start positions  
    for (int i = 0; i < numberAtoms; ++i) {
        electronPos[i] = (float *) malloc(sizeof (float*) * 2);
        electronPos[i][0] = randomPos();
        electronPos[i][1] = randomPos();
    }

    for (int i = 0; i < loops; ++i) {

        // calculate each probability of each atom
        // if thingy change
        // calculate energy???????
        for (int j = 0; j < numberAtoms; ++j) {
            int xComponent = 0;
            int yComponent = 0;
            for (int k = 0; k < numberAtoms; ++ k) {
                if (k == j) {
                    continue;
                }

                const float distance = hypot(electronPos[j][0] - electronPos[i][0], electronPos[j][0] - electronPos[j][0]);

                const float forceMag = 8.987 * pow(10,19) * (pow((1.602*pow(10, -19)), 2) / distance);

                xComponent = xComponent + cos(x2 - x1 / distance) * forceMag;
                yComponent = yComponent + sin(y2 - y1 / distance) * forceMag;
            }
            const resultantForce = hypot(xComponent, yComponent);

            // move the fucker here
        } 

    }

    // return it
    return 1;
}