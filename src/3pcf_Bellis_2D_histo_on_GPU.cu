#include<stdio.h> 
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include <unistd.h>

#include <cuda_runtime.h>

using namespace std;

//#define SUBMATRIX_SIZE 16384
//#define SUBMATRIX_SIZE 2048
//#define SUBMATRIX_SIZE 1024
#define SUBMATRIX_SIZE 512
//#define SUBMATRIX_SIZE 256
//#define SUBMATRIX_SIZE 128
//#define SUBMATRIX_SIZE 64
//#define SUBMATRIX_SIZE 16
//#define SUBMATRIX_SIZE 8

#define NUM_BLOCKS 1

////////////////////////////////////////////////////////////////////////
// Number of histogram bins has to be edited by hand, prior to
// copmilation.
////////////////////////////////////////////////////////////////////////

//#define NUM_BINS 254 
//#define NUM_BINS 126 
//#define NUM_BINS 62 
//#define NUM_BINS 30 
//#define NUM_BINS 64 
#define NUM_BINS 64 // This will be squared eventually.
//#define NUM_BINS 14 
//#define NUM_BINS 6 

#define CONV_FACTOR 57.2957795 // 180/pi

void getDeviceDiagnostics(int tot_Gals, int n_coords);

////////////////////////////////////////////////////////////////////////////////
// Took this code from a Dr. Dobbs example.
////////////////////////////////////////////////////////////////////////////////
void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Diagnostic kernel
////////////////////////////////////////////////////////////////////////////////
__global__ void diagnostic_kernel (int *dev_hist)
{

    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // shared memory
    __shared__ int shared_hist[NUM_BINS];

    // access thread id
    //const unsigned int tid = threadIdx.x;
    // access number of threads in this block
    const unsigned int num_threads = blockDim.x;

    // Note that we only clear things out for the first thread on each block.
    if(threadIdx.x==0)
    {
        for (int i=0;i<NUM_BINS;i++)
            shared_hist[i] = 0;
    }
    __syncthreads();
    ////////////////////////////////////////////////////////////////////////

    // FILL THE ARRAYS

    shared_hist[0] = blockIdx.x*1000 + tid;
    //shared_hist[tid] = blockIdx.x;

    __syncthreads();

    /*
    if(threadIdx.x==0)
    {
        for(int i=0;i<NUM_BINS;i++)
        {
            dev_hist[i+(blockIdx.x*(NUM_BINS))]=shared_hist[i];
        }
    }
    */

}

////////////////////////////////////////////////////////////////////////
// Decide into which bin the data go.
////////////////////////////////////////////////////////////////////////
__device__ int distance_to_bin(float dist, float hist_min, float hist_max, int nbins, float bin_width, int flag)
{
    //int bin_index = -1;
    int bin_index = 1;

    bin_index = floor((dist-hist_min)/bin_width) + 1;
    /*
       if(dist < hist_min)
       bin_index = 0;
       else if(dist >= hist_max)
       bin_index = nbins + 1;
       else
       {
       if (flag==0)
       {
       bin_index = floor((dist-hist_min)/bin_width) + 1;
       }
       else if (flag==1)// log binning
       {
       bin_index = floor((log(dist)-log(hist_min))/bin_width) + 1;
       }
       else if (flag==2)// log 10 binning
       {
       bin_index = floor((log10(dist)-log10(hist_min))/bin_width) + 1;
       }
       }
     */
    return bin_index;
}



////////////////////////////////////////////////////////////////////////
// Kernel to calculate angular distances between galaxies and histogram
// the distances.
////////////////////////////////////////////////////////////////////////
//__global__ void distance(float *x0, float *y0, float *z0, 
__global__ void distance(
        float x_pivot0, float y_pivot0, float z_pivot0,  \
        float x_pivot1, float y_pivot1, float z_pivot1,  \
        float *x2, float *y2, float *z2, \
        int xind, int yind, int zind, \
        int max_xind, int max_yind, int max_zind, \
        int *dev_hist, float hist_min, float hist_max, \
        float bin_width, int flag=0,  float conv_factor_angle=57.2957795)
//__global__ void distance(float x0)
{

    ////////////////////////////////////////////////////////////////////////////
    // Idx will keep track of which thread is being calculated within a given
    // warp.
    ////////////////////////////////////////////////////////////////////////////
    int idx = blockIdx.x * blockDim.x + threadIdx.x; // This should range to SUBMATRIX_SIZE
    int tidx = threadIdx.x; // This just runs over the number of threads in a block.

    //idx += xind;
    //tidx += xind; // Increment by where we are in the big calculation.

    ////////////////////////////////////////////////////////////////////////
    // Shared memory stuff.
    ////////////////////////////////////////////////////////////////////////
    //__shared__ int shared_hist[NUM_BINS*NUM_BINS+1];
    // Note that we only clear things out for the first thread on each block.
    //int hbins = NUM_BINS*NUM_BINS+1;
    //if(threadIdx.x==0)
    //{
    //for (int i=0;i<hbins;i++)
    //shared_hist[i] = 0;
    //}
    //__syncthreads();
    ////////////////////////////////////////////////////////////////////////
    /*

       int i=0;
       int j=0;
       int k=0;

       int i0=0; // shortest
       int i1=0; // middle
       int i2=0; // longest

       float xdiff = 0.0;
       float ydiff = 0.0;
       float zdiff = 0.0;
       float dist0 = 0.0;
       float dist1 = 0.0;
       float dist2 = 0.0;

       int totbin = 0;

       xdiff = x_pivot0-x_pivot1;
       ydiff = y_pivot0-y_pivot1;
       zdiff = z_pivot0-z_pivot1;
       dist0 = sqrtf(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

       xdiff = x_pivot0-x2[idx];
       ydiff = y_pivot0-y2[idx];
       zdiff = z_pivot0-z2[idx];
       dist1 = sqrtf(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

       xdiff = x_pivot0-x2[idx];
       ydiff = y_pivot0-y2[idx];
       zdiff = z_pivot0-z2[idx];
       dist1 = sqrtf(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

    //totbin = 0;

    i0 = distance_to_bin(dist0,hist_min,hist_max,NUM_BINS,bin_width,flag);
    i1 = distance_to_bin(dist1,hist_min,hist_max,NUM_BINS,bin_width,flag);
    i2 = distance_to_bin(dist2,hist_min,hist_max,NUM_BINS,bin_width,flag);

    //totbin = nhistbins2*i2 + nhistbins*i1 + i0;
    totbin = NUM_BINS*i1 + i2;

    shared_hist[0] = i0;

    // THIS SEEMS TO WORK HERE!!!!!!
    //if (j>idx+1 && k>j+1 && idx<max_xind)
    //{
    //int temp = shared_hist[totbin]|1;
    //shared_hist[k] = totbin;
    //shared_hist[threadIdx.x + 8*threadIdx.y + 64*threadIdx.z] = totbin;
    if (totbin<hbins-1)
    atomicAdd(&shared_hist[totbin+1],1);
    //}

    __syncthreads();


    // This works if we only have one block.
    if(threadIdx.x==0)
    {
    //for(int i=0;i<tot_hist_size;i++)
    for(int i=0;i<hbins;i++)
    {
    dev_hist[i] = shared_hist[i];
    }
    }
     */
}

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    // Needed for parsing command-line arguments.
    extern char *optarg;
    extern int optind, optopt, opterr;
    int c;
    char *outfilename = NULL;
    char defaultoutfilename[256];
    sprintf(defaultoutfilename,"default_out.dat");

    int nbins = NUM_BINS;
    int log_binning_flag = 0; // False

    float scale_factor = 1.0; // For if we need to convert input to arcsec or arcmin
    float conv_factor_angle = 57.2957795; // 180/pi // For if we need to convert arcdistance to arcsec or arcmin

    bool silent_on_GPU_testing = false;

    float hist_min = 0;
    //float hist_max = 1.8;
    float hist_max = 7000.0;
    float bin_width = (hist_max-hist_min)/NUM_BINS;
    float hist_bin_width = bin_width; // For now
    int flag = 0;

    float hist_lower_range = 0.0;
    float hist_upper_range = hist_max;

    cudaError_t error = cudaGetLastError();

    int NBINS2 = NUM_BINS*NUM_BINS;
    int NBINS2p1 = NBINS2 + 1;
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    while ((c = getopt(argc, argv, "ao:L:l:w:smS")) != -1) {
        switch(c) {
            case 'L':
                printf("L is set\n");
                hist_lower_range = atof(optarg);
                break;
            case 'w':
                hist_bin_width = atof(optarg);
                printf("Histogram bin width: %f\n",hist_bin_width);
                break;
            case 'l':
                log_binning_flag = atoi(optarg);
                printf("Will use binning %d.\n",log_binning_flag);
                break;
            case 's':
                scale_factor = 206264.0; // To convert arcseconds to radians.
                conv_factor_angle *= 3600.0; // convert radians to arcseconds.
                printf("Reading in values assuming they are arcseconds.\n");
                printf("scale_factor: %f\n",scale_factor);
                printf("conv_factor_angle: %f\n",conv_factor_angle);
                break;
            case 'm':
                scale_factor = 3437.74677; // To convert arcminutes to radians.
                conv_factor_angle *= 60.0; // convert radians to arcminutes.
                printf("scale_factor: %f\n",scale_factor);
                printf("conv_factor_angle: %f\n",conv_factor_angle);
                printf("Reading in values assuming they are arcminutes.\n");
                break;
            case 'o':
                outfilename = optarg;
                printf("Output filename is %s\n", outfilename);
                break;
            case 'S':
                printf("Silent mode - don't run the GPU test (suppresses some output)");
                silent_on_GPU_testing = true;
                break;
            case '?':
                printf("unknown arg %c\n", optopt);
                break;
        }
    }

    if (argc < 3)
    {

        printf("\nMust pass in at least three input files on command line!\n");
        printf("\nUsage: ", argv[0] );
        exit(1);
    }

    // Set a default output file name, if none was passed in on the 
    // command line.
    if (outfilename == NULL) 
    {
        outfilename = defaultoutfilename;
        printf("Output filename is %s\n", outfilename);
    }

    float *h_x[3], *h_y[3], *h_z[3];
    float *d_x[3], *d_y[3], *d_z[3];

    // Open the input files and the output file.
    FILE *infile[3], *outfile;
    for (int i=0;i<3;i++)
    {
        infile[i] = fopen(argv[optind+i],"r");
        printf("Opening input file %d: %s\n",i,argv[optind+i]);
    }
    outfile = fopen(outfilename, "w");

    // Determine which combination of files we are doing.
    // 0 - all the same (DDD or RRR)
    // 1 - first one is different (DRR or RDD)
    // 2 - middle one is different (DRD or RDR)
    // 3 - last one is different (RRD or DDR)
    bool which_three_input_files = 0;
    if (strcmp(argv[optind+0],argv[optind+1])==0 && strcmp(argv[optind+0],argv[optind+2])==0)
    {
        which_three_input_files = 0;
        printf("Using the same file! (DDD or RRR)\n");
    }
    else if (strcmp(argv[optind+0],argv[optind+1])!=0 && strcmp(argv[optind+1],argv[optind+2])==0)
    {
        which_three_input_files = 1;
        printf("Not the same file! (DRR or RDD)\n");
    }
    else if (strcmp(argv[optind+0],argv[optind+2])!=0 && strcmp(argv[optind+0],argv[optind+2])==0)
    {
        which_three_input_files = 2;
        printf("Not the same file! (DRD or RDR)\n");
    }
    else if (strcmp(argv[optind+0],argv[optind+1])==0 && strcmp(argv[optind+0],argv[optind+2])!=0)
    {
        which_three_input_files = 2;
        printf("Not the same file! (RRD or DDR)\n");
    }


    //////////////////////////////////////////////////////////////////////
    // Read in the galaxy files.
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // Read in the files.
    ////////////////////////////////////////////////////////////////////////////

    int NUM_GALAXIES[3] = {0.0,0.0,0.0};
    int size_of_galaxy_array[3] = {0.0,0.0,0.0};
    float temp0, temp1, temp2;

    for (int i=0;i<3;i++)
    {
        fscanf(infile[i], "%d", &NUM_GALAXIES[i]);
        size_of_galaxy_array[i] = NUM_GALAXIES[i] * sizeof(float);
        printf("SIZE %d # GALAXIES: %d\n",i,NUM_GALAXIES[i]);

        h_x[i] = (float*)malloc(size_of_galaxy_array[i]);
        h_y[i] = (float*)malloc(size_of_galaxy_array[i]);
        h_z[i] = (float*)malloc(size_of_galaxy_array[i]);

        for(int j=0; j<NUM_GALAXIES[i]; j++)
        {
            fscanf(infile[i], "%f %f %f", &temp0, &temp1, &temp2);
            h_x[i][j] = temp0/scale_factor;
            h_y[i][j] = temp1/scale_factor;
            h_z[i][j] = temp2/scale_factor;
            //if (j<10 || j>NUM_GALAXIES[i]-10)
            //if (i==0)
            //printf("%e %e %e\n", h_x[i][j],h_y[i][j],h_z[i][j]);
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Now get the info from the device.
    ////////////////////////////////////////////////////////////////////////////
    if (!silent_on_GPU_testing)
    {
        getDeviceDiagnostics(NUM_GALAXIES[0], NUM_GALAXIES[0]);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Define the grid and block size
    ////////////////////////////////////////////////////////////////////////////
    dim3 grid, block;
    grid.x = NUM_BLOCKS; // Is this the number of blocks?
    block.x = SUBMATRIX_SIZE; // The number of threads per block? 

    printf("# of blocks per grid:  grid.x:  %d\n",grid.x);
    printf("# of thread per block: block.x: %d\n",block.x);
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // Allocation of histogram
    ///////////////////////////////////////////////////////////////////////////
    // FROM C CODE /////////
    int *hist;
    int *dev_hist;

    ////////////////////////////////////////////////////////////////////////////
    // This is the total number of bins/bytes we need for all of the 
    // different histograms we'll be creating
    ////////////////////////////////////////////////////////////////////////////
    int tot_nbins = (NUM_BINS)*(NUM_BINS)*(NUM_BINS);
    int tot_nbins_on_gpu = (NUM_BINS)*(NUM_BINS)+1; // This will be a 2D hist

    int size_hist = grid.x*tot_nbins;
    int size_hist_bytes = size_hist*sizeof(int);

    int size_hist_on_gpu = grid.x*tot_nbins_on_gpu;
    int size_hist_bytes_on_gpu = size_hist_on_gpu*sizeof(int);

    printf("Size of histogram: %d bins\t%d bytes\n",size_hist,size_hist_bytes);
    printf("Size of histogram on GPU: %d bins\t%d bytes\n",size_hist_on_gpu,size_hist_bytes_on_gpu);

    hist = (int*)malloc(size_hist_bytes);
    memset(hist, 0, size_hist_bytes);

    cudaMalloc((void **) &dev_hist, (size_hist_bytes_on_gpu));
    error = cudaGetLastError();
    printf("ERROR: %s\n", cudaGetErrorString(error) );

    printf("dev_hist: %x\n",dev_hist);
    cudaMemset(dev_hist, 0, size_hist_bytes_on_gpu);
    error = cudaGetLastError();
    printf("ERROR: %s\n", cudaGetErrorString(error) );

    printf("dev_hist bins: %d\n",size_hist_on_gpu);
    printf("dev_hist size: %d\n",size_hist_bytes_on_gpu);

    //exit(0);

    int *summed_hist;

    int summed_hist_size = tot_nbins * sizeof(int);
    summed_hist =  (int*)malloc(summed_hist_size);
    printf("Size of summed histogram array: %d bins\t%d bytes\n",(tot_nbins),summed_hist_size);
    memset(summed_hist,0,summed_hist_size); 

    // Check to see if we allocated enough memory.
    for (int i=0;i<3;i++)
    {
        cudaMalloc((void **) &d_x[i], size_of_galaxy_array[i] );
        cudaMalloc((void **) &d_y[i], size_of_galaxy_array[i] );
        cudaMalloc((void **) &d_z[i], size_of_galaxy_array[i] );

        //printf("Size: %d\n",size_of_galaxy_array[i]);
        if (0==d_x[i] || 0==d_y[i] || 0==d_z[i])
        {
            printf("Couldn't allocate memory on the GPU for the galaxy arrays!\n");
            printf("Size: %d\n",size_of_galaxy_array[i]);
            return 1;
        }

        cudaMemset(d_x[i],0,size_of_galaxy_array[i]);
        cudaMemset(d_y[i],0,size_of_galaxy_array[i]);
        cudaMemset(d_z[i],0,size_of_galaxy_array[i]);

        //printf("gal: %f\n",h_x[i][0]);
        cudaMemcpy(d_x[i], h_x[i], size_of_galaxy_array[i], cudaMemcpyHostToDevice );
        cudaMemcpy(d_y[i], h_y[i], size_of_galaxy_array[i], cudaMemcpyHostToDevice );
        cudaMemcpy(d_z[i], h_z[i], size_of_galaxy_array[i], cudaMemcpyHostToDevice );

    } 
    checkCUDAError("mem alloc");

    if ( 0==dev_hist )
    {
        printf("Couldn't allocate memory for the histogram on the GPU!\n");
        return 1;
    }
    checkCUDAError("mem alloc dev hist");

    int X, Y;
    //int num_submatrices_x = NUM_GALAXIES0 / SUBMATRIX_SIZE;
    //int num_submatrices_y = NUM_GALAXIES1 / SUBMATRIX_SIZE;
    int num_submatrices[3] = {0, 0, 0};
    for (int i=0;i<3;i++)
    {
        //num_submatrices[i] = NUM_GALAXIES[i]/SUBMATRIX_SIZE;
        num_submatrices[i] = NUM_GALAXIES[i]/(block.x*grid.x);
        // Take care of edges of matrix.
        //if (NUM_GALAXIES[i]%SUBMATRIX_SIZE != 0)
        if (NUM_GALAXIES[i]%(block.x*grid.x) != 0)
        {
            num_submatrices[i] += 1;
        }
    }


    printf("Breaking down the calculations.\n");
    printf("Number of submatrices: %dx%dx%d\n",num_submatrices[0],num_submatrices[1],num_submatrices[2]);
    printf("Number of calculations per submatrices: %dx%d\n",SUBMATRIX_SIZE,SUBMATRIX_SIZE);

    int xind, yind, zind;
    int bin_index = 0;
    unsigned long tot0 = 0;
    //for(int i = 0; i < num_submatrices[0]; i++)
    for(int i = 0; i < NUM_GALAXIES[0]; i++)
    {
        if (i%100==0)
        {
            printf("%d\n",i);
            fflush(stdout);
        }
        //xind = i*SUBMATRIX_SIZE;
        //xind = i*(block.x*grid.x);
        xind = i;
        int jmin = 0;
        if (which_three_input_files==0) // DDD or RRR
            jmin = 0;
        else if (which_three_input_files==1) // DRR or RDD
            jmin = 0;
        else if (which_three_input_files==2) // DRD or RDR
            jmin = 0;
        else if (which_three_input_files==3) // DDR or RRD
            jmin = 0;
        //for(int j = jmin; j < num_submatrices[1]; j++)
        int max_x = NUM_GALAXIES[0];
        printf("xind: %5d %5d\n",xind,max_x);
        for(int j = jmin; j < NUM_GALAXIES[1]; j++)
        {
            //yind = j*SUBMATRIX_SIZE;
            yind = j;
            //yind = j*(block.x*grid.x);
            int kmin = 0;
            if (which_three_input_files==0)
                kmin = 0;
            else if (which_three_input_files==1)
                kmin = 0;
            else if (which_three_input_files==2)
                kmin = 0;
            else if (which_three_input_files==3)
                kmin = 0;
            int max_y = NUM_GALAXIES[1];
            //printf("yind: %5d %5d\n",yind,max_y);
            for(int k =kmin; k < num_submatrices[2]; k++)
            {
                //zind = k*SUBMATRIX_SIZE;
                zind = k*(block.x*grid.x);
                //bool do_calc = 1;
                //if (do_calc)
                {

                    // Set the histogram to all zeros each time.
                    cudaMemset(dev_hist,0,size_hist_bytes_on_gpu);
                    checkCUDAError("memset");

                    int max_z = NUM_GALAXIES[2];

                    //printf("zind: %5d %5d\n",zind,max_z);
                    //printf("nbins: %d\n",nbins);
                    //distance<<<grid,block>>>(h_x[0],h_y[0],h_z[0], 
                    //printf("%f %f %f\n",h_x[0][i],h_y[0][i],h_z[0][i]);
                    //printf("%f %f %f\n",h_x[0][j],h_y[0][j],h_z[0][j]);
                       distance<<<grid,block>>>(
                       d_x[0][i],d_y[0][i],d_z[0][i],\
                       d_x[1][j],d_y[1][j],d_z[1][j],\
                       d_x[2],d_y[2],d_z[2],\
                       xind, yind, zind, \
                       max_x, max_y, max_z,\
                       dev_hist, hist_lower_range, hist_upper_range, \
                       hist_bin_width, log_binning_flag, conv_factor_angle);
                    //diagnostic_kernel<<<grid,block>>>(dev_hist);

                    checkCUDAError("kernel");
                    //printf("there\n");
                    //printf("dev_hist: %x\n",dev_hist);

                    cudaMemcpy(hist, dev_hist, size_hist_bytes_on_gpu, cudaMemcpyDeviceToHost);
                    checkCUDAError("memcpy");

                    ////////////////////////////////////////////////////////////////////
                    // Sum up the histograms from each thread (hist).
                    ////////////////////////////////////////////////////////////////////
                    {
                        for(int m=1; m<NBINS2p1; m++)
                        {
                            //bin_index = hist[0]*NBINS2 + m;
                            //summed_hist[bin_index] += hist[m];
                            //tot0 += hist[m];
                        }    
                    }    
                }
            }  
        }
    }

    cudaMemcpy(hist, dev_hist, size_hist_bytes, cudaMemcpyDeviceToHost);

    unsigned long total = 0;
    int index = 0;
    fprintf(outfile,"%d %d %d\n",NUM_BINS,NUM_BINS,NUM_BINS);
    for(int i = 0; i < NUM_BINS; i++)
    {
        //printf("%d --------------\n",i);
        fprintf(outfile,"%d\n",i);
        for(int j = i; j < NUM_BINS; j++)
        {
            for(int k = j; k < NUM_BINS; k++)
            {
                index = i*NUM_BINS*NUM_BINS + j*NUM_BINS + k;
                //printf("%ul ",summed_hist[index]);
                fprintf(outfile,"%ul ",summed_hist[index]);
                total += summed_hist[index];
            }
            fprintf(outfile,"\n");
            //printf("\n");
        }
    }
    printf("total: %ul\n",total);
    printf("tot0: %ul\n",tot0);

    ////////////////////////////////////////////////////////////////////////////
    // Close out the files and free up the memory.
    ////////////////////////////////////////////////////////////////////////////
    fclose(infile[0]);
    fclose(infile[1]);
    fclose(infile[2]);
    fclose(outfile);

    for (int i=0;i<3;i++)
    {
        free(h_x[i]);
        free(h_y[i]);
        free(h_z[i]);

        cudaFree(d_x[i]);
        cudaFree(d_y[i]);
        cudaFree(d_z[i]);
    }
    free(hist);
    cudaFree(dev_hist);

    return 0;
}  
//////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void getDeviceDiagnostics(int tot_gals, int n_coords)
{
    printf("\n------ CUDA device diagnostics ------\n\n");

    //int tot_gals = NUM_GALAXIES[0];
    //int tot_gals = 1000;
    int nx = SUBMATRIX_SIZE;
    //int ncalc = nx * nx;
    int ncalc = nx * nx;
    int gpu_mem_needed = int(tot_gals * sizeof(float)) * 9; // need to allocate ra, dec.
    printf("Requirements: %d calculations and %d bytes memory on the GPU \n\n", ncalc, gpu_mem_needed);

    int deviceCount = 0;
    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
    if (error_id != cudaSuccess) {
        printf( "cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id) );
    }
    // This function call returns 0 if there are no CUDA capable devices.
    if (deviceCount == 0)
        printf("There is no device supporting CUDA\n");
    else
        printf("Found %d CUDA Capable device(s)\n", deviceCount);


    int dev=0;
    for (dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
        printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

        printf("  Total amount of global memory:                 %.0f MBytes (%llu bytes)\n",
                (float)deviceProp.totalGlobalMem/1048576.0f, (int) deviceProp.totalGlobalMem);


        printf("  Warp size:                                     %d\n", deviceProp.warpSize);
        printf("  Maximum number of threads per block:           %d\n", deviceProp.maxThreadsPerBlock);
        printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
                deviceProp.maxThreadsDim[0],
                deviceProp.maxThreadsDim[1],
                deviceProp.maxThreadsDim[2]);
        printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
                deviceProp.maxGridSize[0],
                deviceProp.maxGridSize[1],
                deviceProp.maxGridSize[2]);

        // does this device have enough capcacity for the calculation?
        printf("\n*************\n");

        // check memory
        if((int) deviceProp.totalGlobalMem < gpu_mem_needed) printf(" FAILURE: Not eneough memeory on device for this calculation! \n");
        else
        {
            printf("Hurrah! This device has enough memory to perform this calculation\n");

            // check # threads

            int threadsPerBlock = deviceProp.maxThreadsPerBlock; // maximal efficiency exists if we use max # threads per block.
            int blocksPerGrid = int(ceil(ncalc / threadsPerBlock)); // need nx*nx threads total
            if(deviceProp.maxThreadsDim[0] >blocksPerGrid) printf("FAILURE: Not enough threads on the device to do this calculation!\n");
            else
            {
                printf("Hurrah! This device supports enough threads to do this calculation\n");
                // how many kernels can we run at once on this machine?
                int n_mem = floor(deviceProp.totalGlobalMem / float(gpu_mem_needed));
                int n_threads = floor(threadsPerBlock * deviceProp.maxThreadsDim[0]*deviceProp.maxThreadsDim[1] / float(ncalc) ); // max # threads possible?

                printf("%d %d  \n",  n_threads, deviceProp.maxThreadsDim[0]);

                int max_kernels = 0;
                n_mem<n_threads ? max_kernels = n_mem : max_kernels = n_threads;

                printf(" you can run %d kernels at a time on this device without overloading the resources \n", max_kernels);
            }
        }

    }

    printf("\n------ End CUDA device diagnostics ------\n\n");

}
