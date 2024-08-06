//#define ENABLE_HDF

#ifdef ENABLE_HDF
    #include "hdf5.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N               100  //number of nuclei
#define T               1000 //number of timesteps
#define DECAY_CHANCE    0.001 //probability of decay per time step

void printTimeSeries(double *particles, int timeCount, int particleCount) 
{
    //We assume it's a square grid
    int gs = (int)sqrt(1.0*particleCount);
    for (int t=0; t<T; t++) 
    {
        printf("Time step %d\n", t);
        for (int i=0; i<gs; i++) 
        {
            for (int j=0; j<gs; j++) 
            {
                int idx = (t*particleCount+i*gs + j);
                printf("%.4f ", particles[idx]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

#ifdef ENABLE_HDF
int mat2hdf5 (double* wdata, int timeCount, int particleCount, const char *FILE, const char *DATASET) 
{
    hid_t       file, space, dset;          /* Handles */
    herr_t      status;
    hsize_t     dims[2] = {timeCount,particleCount};

    /*
     * Create a new file using the default properties.
     */
    file = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    space = H5Screate_simple (2, dims, NULL);

    /*
     * Create the dataset and write the floating point data to it.  In
     * this example we will save the data as 64 bit little endian IEEE
     * floating point numbers, regardless of the native type.  The HDF5
     * library automatically converts between different floating point
     * types.
     */
    dset = H5Dcreate (file, DATASET, H5T_IEEE_F64LE, space, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                wdata);

    /*
     * Close and release resources.
     */
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Fclose (file);

    return 0;
}
#endif

int main (int argc, char **argv) 
{
    //Every nuclei for every time step. Does it exist (1), or has it decayed (0).
    //This would make sense to be an int, but we'll use float64 here as you'll probably want to deal with that datatype more often in your projects etc.
    double *nucleiAll = (double*)malloc(N*T*sizeof(double));
    
    float random_float;

    //Initialise particles to exist (1)
    for (int i=0;i<N;i++)
    {
        nucleiAll[i] = 1;
    }

    for (int t=1; t<T; t++) 
    {
        double *nucleiNow = &nucleiAll[t*N];

        double *nucleiBefore= &nucleiAll[(t-1)*N];

        for (int i=0; i<N; i++) 
        {
                //random float from 0 to 1
                random_float = (float)rand()/(float)(RAND_MAX); 
                if (random_float<DECAY_CHANCE) 
                { // Now I am become Death. Decay, particle no longer exists.
                    nucleiNow[i] = 0;
                } 
                else 
                { // no decay
                    nucleiNow[i] = nucleiBefore[i];
                }
        }
    }

    #ifdef ENABLE_HDF
        mat2hdf5 (nucleiAll,T,N,"particles.h5","DS1");
    #else
        printTimeSeries(nucleiAll,T,N);
    #endif
    free(nucleiAll);
    return 0;
}
