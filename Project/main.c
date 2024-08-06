#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double **grid(int rows, int columns) {
    double **assembledGrid;
    assembledGrid = malloc(sizeof (double) * rows);

    for (int i = 0; i < rows; ++i) {
        assembledGrid[i] = malloc(sizeof (double) * columns);
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; i < columns; ++i) {
            assembledGrid[i][j] = 0.;
        }
    }

    return assembledGrid;
}
 
int main(int width, int height, int loops) {

    // assemble grid
    // add initial atom state
    int numberAtoms = 4;
    double** atomGrid = grid(10, 10);

    printf("%f", atomGrid);

    float electronPos[100] = []

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
                    // don't need to consider it's own force on its self
                }

                const distance = hypot(x2 - x1, y2 - y1);

                const forceMag = 8.987 * pow(10,19) (pow((1.602*pow(10, -19), 2)) / distance);

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