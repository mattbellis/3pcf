#include<stdio.h> 
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include <unistd.h>

#include <cuda_runtime.h>

using namespace std;

//#define SUBMATRIX_SIZE 16384
//#define SUBMATRIX_SIZE 4096
//#define SUBMATRIX_SIZE 2048
#define SUBMATRIX_SIZE 1024
//#define SUBMATRIX_SIZE 512
//#define SUBMATRIX_SIZE 256
//#define SUBMATRIX_SIZE 128
//#define SUBMATRIX_SIZE 16

////////////////////////////////////////////////////////////////////////
// Number of histogram bins has to be edited by hand, prior to
// copmilation.
////////////////////////////////////////////////////////////////////////

//#define DEFAULT_NBINS 254 
//#define DEFAULT_NBINS 126 
//#define DEFAULT_NBINS 62 
//#define DEFAULT_NBINS 30 
#define DEFAULT_NBINS 16
//#define DEFAULT_NBINS 8
//#define DEFAULT_NBINS 6 

#define CONV_FACTOR 57.2957795 // 180/pi

void getDeviceDiagnostics(int tot_Gals, int n_coords);

////////////////////////////////////////////////////////////////////////
// Decide into which bin the data go.
////////////////////////////////////////////////////////////////////////
__device__ int distance_to_bin(float dist, float hist_min, float hist_max, int nbins, float bin_width, int flag)
{
    int bin_index = 0;

    /*
    int bin_index = -1;
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
    //bin_index = floor((dist-hist_min)/bin_width) + 1;
    bin_index = floor((dist-hist_min)/bin_width);
    return bin_index;
}



////////////////////////////////////////////////////////////////////////
// Kernel to calculate angular distances between galaxies and histogram
// the distances.
////////////////////////////////////////////////////////////////////////
//__global__ void distance(float *x0, float *y0, float *z0, 
__global__ void distance(
        float *x0, float *y0, float *z0, \
        float *x1, float *y1, float *z1, \
        float *x2, float *y2, float *z2, \
        int xind, int yind, int zind, \
        int max_xind, int max_yind, int max_zind, \
        int *dev_hist, float hist_min, float hist_max, int nbins, \
        float bin_width, int flag=0,  float conv_factor_angle=57.2957795)
//__global__ void distance(float x0)
{

    ////////////////////////////////////////////////////////////////////////////
    // Idx will keep track of which thread is being calculated within a given 
    // warp.
    ////////////////////////////////////////////////////////////////////////////
    int idx = blockIdx.x * blockDim.x + threadIdx.x; // This should range to SUBMATRIX_SIZE

    int tot_hist_size = (DEFAULT_NBINS)*(DEFAULT_NBINS)*(DEFAULT_NBINS);
    //int tot_hist_size = ((DEFAULT_NBINS+2+2)*(DEFAULT_NBINS+2+1)*(DEFAULT_NBINS+2)/6);

    int idxorg = idx;

    idx += xind;

    ////////////////////////////////////////////////////////////////////////
    // Shared memory stuff.
    ////////////////////////////////////////////////////////////////////////
    __shared__ int shared_hist[(DEFAULT_NBINS)*(DEFAULT_NBINS)*(DEFAULT_NBINS)];
    //__shared__ int shared_hist[(DEFAULT_NBINS+2+2)*(DEFAULT_NBINS+2+1)*(DEFAULT_NBINS+2)/6];
    // Note that we only clear things out for the first thread on each block.
    if(threadIdx.x==0)
    {
        for (int i=0;i<tot_hist_size;i++)
            shared_hist[i] = 0;
    }
    __syncthreads();
    ////////////////////////////////////////////////////////////////////////

    int DEFAULT_NBINS2 = DEFAULT_NBINS*DEFAULT_NBINS;

    float x0i = 0.0;
    float y0i = 0.0;
    float z0i = 0.0;

    int i=0;
    int j=0;
    int k=0;

    //bool do_calc = 1;

    int ymax = yind + SUBMATRIX_SIZE;
    int zmax = zind + SUBMATRIX_SIZE;

    ymax = (ymax>max_yind ? max_yind : ymax);
    zmax = (zmax>max_zind ? max_zind : zmax);

    int b0 = false;
    int b1 = false;
    int b2 = false;

    int notb0 = true;
    int notb1 = true;
    int notb2 = true;

    int i0=0; // shortest
    int i1=0; // middle
    int i2=0; // longest

    int ti0=0; // shortest
    int ti1=0; // middle
    int ti2=0; // longest

    float xdiff = 0.0;
    float ydiff = 0.0;
    float zdiff = 0.0;
    float dist0 = 0.0;
    float dist1 = 0.0;
    float dist2 = 0.0;

    float xpt0,ypt0,zpt0;
    float xpt1,ypt1,zpt1;
    float xpt2,ypt2,zpt2;

    //int nhistbins = nbins+2;
    //int nhistbins2 = nhistbins*nhistbins;
    int totbin = 0;

    if (idx<max_xind)
    {
        xpt0 = x0[idx];
        ypt0 = y0[idx];
        zpt0 = z0[idx];

        # pragma unroll
        for(j=yind; j<ymax; j++)
        //for(j=idx+1; j<ymax; j++)
        {
            xpt1 = x1[j];
            ypt1 = y1[j];
            zpt1 = z1[j];

            //xdiff = x0[idx]-x1[j];
            //ydiff = y0[idx]-y1[j];
            //zdiff = z0[idx]-z1[j];
            xdiff = xpt0-xpt1;
            ydiff = ypt0-ypt1;
            zdiff = zpt0-zpt1;
            dist0 = sqrtf(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

            # pragma unroll
            for(k=zind; k<zmax; k++)
            //for(k=j+1; k<zmax; k++)
            {
                    xpt2 = x2[k];
                    ypt2 = y2[k];
                    zpt2 = z2[k];

                    //xdiff = x0[idx]-x2[k];
                    //ydiff = y0[idx]-y2[k];
                    //zdiff = z0[idx]-z2[k];
                    xdiff = xpt0-xpt2;
                    ydiff = ypt0-ypt2;
                    zdiff = zpt0-zpt2;
                    dist1 = sqrtf(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

                    //xdiff = x1[j]-x2[k];
                    //ydiff = y1[j]-y2[k];
                    //zdiff = z1[j]-z2[k];
                    xdiff = xpt1-xpt2;
                    ydiff = ypt1-ypt2;
                    zdiff = zpt1-zpt2;
                    dist2 = sqrtf(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

                    i0=0; // shortest
                    i1=0; // middlest
                    i2=0; // longest

                    totbin = 0;

                    i0 = distance_to_bin(dist0,hist_min,hist_max,nbins,bin_width,flag);
                    i1 = distance_to_bin(dist1,hist_min,hist_max,nbins,bin_width,flag);
                    i2 = distance_to_bin(dist2,hist_min,hist_max,nbins,bin_width,flag);

                    // Longest, then middle, then shortest
                    //totbin = nhistbins2*((b1*b2)*i2  + (b0*(notb1))*i1 + ((notb0)*(notb2))*i0) +  \
                    nhistbins*((notb1*b2 + b1*notb2)*i2 + (b0*b1 + (notb0*notb1))*i1 + (notb0*b2 + b0*notb2)*i0) +  \
                        ((notb1*notb2)*i2  + (notb0*(b1))*i1 + ((b0)*(b2))*i0);
                    //ti2 = ((b1*b2)*i2  + (b0*(notb1))*i1 + ((notb0)*(notb2))*i0) ;
                    //ti1 = ((notb1*b2 + b1*notb2)*i2 + (b0*b1 + (notb0*notb1))*i1 + (notb0*b2 + b0*notb2)*i0);
                    //ti0 = ((notb1*notb2)*i2  + (notb0*(b1))*i1 + ((b0)*(b2))*i0);

                    //totbin = ti2 + (ti1*(ti1+1))/2 + (ti0*(ti0+1)*(ti0+2))/6;
                    //totbin = (ti0)*(3*(nbins+2)*((nbins+2)+1)-(3*(nbins+2)+2)*(ti0+1) + (ti0+1)*(ti0+1))/6 + (ti1)*(2*(nbins+2)-(ti1+1))/2 + ti2;
                    totbin = DEFAULT_NBINS2*i2 + DEFAULT_NBINS*i1 + i0;
                    //totbin = 1.0;

                    //totbin = nhistbins2*i2 + nhistbins*i1 + i0;

                    /*
                       if (b0==true && b1==true)
                       {
                       i0 = distance_to_bin(dist0,hist_min,hist_max,nbins,bin_width,flag);
                       i1 = distance_to_bin(dist1,hist_min,hist_max,nbins,bin_width,flag);
                       i2 = distance_to_bin(dist2,hist_min,hist_max,nbins,bin_width,flag);
                       }
                       else if (b0==false && b1==false)
                       {
                       i0 = distance_to_bin(dist2,hist_min,hist_max,nbins,bin_width,flag);
                       i1 = distance_to_bin(dist1,hist_min,hist_max,nbins,bin_width,flag);
                       i2 = distance_to_bin(dist0,hist_min,hist_max,nbins,bin_width,flag);
                       }
                       else if (b0==true && b1==false && b2==true)
                       {
                       i0 = distance_to_bin(dist0,hist_min,hist_max,nbins,bin_width,flag);
                       i1 = distance_to_bin(dist2,hist_min,hist_max,nbins,bin_width,flag);
                       i2 = distance_to_bin(dist1,hist_min,hist_max,nbins,bin_width,flag);
                       }
                       else if (b0==false && b1==true && b2==true)
                       {
                       i0 = distance_to_bin(dist1,hist_min,hist_max,nbins,bin_width,flag);
                       i1 = distance_to_bin(dist0,hist_min,hist_max,nbins,bin_width,flag);
                       i2 = distance_to_bin(dist2,hist_min,hist_max,nbins,bin_width,flag);
                       }
                       else if (b0==true && b1==false && b2==false)
                       {
                       i0 = distance_to_bin(dist2,hist_min,hist_max,nbins,bin_width,flag);
                       i1 = distance_to_bin(dist0,hist_min,hist_max,nbins,bin_width,flag);
                       i2 = distance_to_bin(dist1,hist_min,hist_max,nbins,bin_width,flag);
                       }
                       else if (b0==false && b1==true && b2==false)
                       {
                       i0 = distance_to_bin(dist2,hist_min,hist_max,nbins,bin_width,flag);
                       i1 = distance_to_bin(dist1,hist_min,hist_max,nbins,bin_width,flag);
                       i2 = distance_to_bin(dist0,hist_min,hist_max,nbins,bin_width,flag);
                       }
                     */

                    //i0 = distance_to_bin(dist0,hist_min,hist_max,nbins,bin_width,flag);
                    //i1 = distance_to_bin(dist1,hist_min,hist_max,nbins,bin_width,flag);
                    //i2 = distance_to_bin(dist2,hist_min,hist_max,nbins,bin_width,flag);
                    /*
                       i0 = 1;
                       if (dist2<hist_max)
                       i0 = 0;
                       else 
                       i0 = 1;
                     */
                    //i0 = int((dist0-hist_min)/bin_width) + 1;
                    //nhistbins = nbins+2;
                    //nhistbins2 = nhistbins*nhistbins;
                    //totbin = nhistbins2*i2 + nhistbins*i1 + i0;
                    //i0 = 1;
                    //totbin = nhistbins2*i2 + nhistbins*i1 + i0;
                    //totbin = 2047;

                    //if(threadIdx.x==0 && blockIdx.x==0)
                    //if(idx==0)
                    //atomicExch(&dev_hist[0],i0);
                    //atomicAdd(&dev_hist[totbin],1);
                    //dev_hist[0]=totbin;
                    //atomicAdd(&dev_hist[0],totbin);
                    //atomicAdd(&dev_hist[totbin + (blockIdx.x*tot_hist_size)],1);
                    //dev_hist[totbin + (idxorg*tot_hist_size)]++;
                    //dev_hist[idxorg] = totbin;
                    //atomicAdd(&dev_hist[totbin],1);
                    //atomicAdd(&dev_hist[0],1);
                    //shared_hist[totbin]++;
                    //shared_hist[totbin] = 1.0;
                    //atomicAdd(&shared_hist[0],1);

                    //int temp = dev_hist[totbin + (blockIdx.x*tot_hist_size)]|1;
                    //dev_hist[totbin + (blockIdx.x*tot_hist_size)] = temp;

                    //int temp = shared_hist[totbin]|1;
                    //shared_hist[totbin] = temp;
                    //shared_hist[threadIdx.x] = totbin;
                    //shared_hist[threadIdx.x];
                    //shared_hist[totbin]++;

                    // THIS SEEMS TO WORK HERE!!!!!!
                    //if (j>idx+1 && k>j+1 && idx<max_xind)
                    if (totbin>=0 && totbin<tot_hist_size)
                    {
                        //int temp = shared_hist[totbin]|1;
                        //shared_hist[threadIdx.x] = totbin;
                        shared_hist[totbin]++;
                        //atomicAdd(&shared_hist[totbin],1);
                    }

            }
        }
    }

    __syncthreads();

    if(threadIdx.x==0)
    {
        for(int i=0;i<tot_hist_size;i++)
        {
            dev_hist[i+(blockIdx.x*tot_hist_size)]=shared_hist[i];
        }
    }
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

    int nbins = DEFAULT_NBINS;
    int log_binning_flag = 0; // False

    float scale_factor = 1.0; // For if we need to convert input to arcsec or arcmin
    float conv_factor_angle = 57.2957795; // 180/pi // For if we need to convert arcdistance to arcsec or arcmin

    bool silent_on_GPU_testing = false;

    float hist_min = 0;
    //float hist_max = 1.8;
    //float hist_max = 7000.0;
    float hist_max = sqrt(3.0*(24*24)); // For the nearest 1k in Wechsler
    //float hist_max = sqrt(3.0*(48.6*48.6)); // For the nearest 10k in Wechsler
    float bin_width = (hist_max-hist_min)/nbins;
    float hist_bin_width = bin_width; // For now
    int flag = 0;

    float hist_lower_range = 0.0;
    float hist_upper_range = hist_max;

    cudaError_t error = cudaGetLastError();
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
    // 128*4 = 512, the amount of memory needed for one histogram.
    // 8192*4 = 32768 is max memory to ask for for the histograms.
    // 8192/128 = 64, is is the right number of blocks?
    //grid.x = 8192/(tot_nbins); // Is this the number of blocks?
    //grid.x = 32; // Is this the number of blocks?
    grid.x = 8; // Is this the number of blocks?
    //grid.x = 4; // Is this the number of blocks?
    block.x = SUBMATRIX_SIZE/grid.x; // Is this the number of threads per block? NUM_GALAXIES/block.x;
    //block.x = SUBMATRIX_SIZE; // Is this the number of threads per block? NUM_GALAXIES/block.x;
    // SUBMATRIX is the number of threads per warp? Per kernel call?
    printf("# of blocks per grid:  grid.x:  %d\n",grid.x);
    printf("# of thread per block: block.x: %d\n",block.x);
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // Allocation of histogram
    ///////////////////////////////////////////////////////////////////////////
    // FROM C CODE /////////
    int *hist;
    int *dev_hist;
    printf("dev_hist: %x\n",dev_hist);

    ////////////////////////////////////////////////////////////////////////////
    // This is the total number of bins/bytes we need for all of the 
    // different histograms we'll be creating
    ////////////////////////////////////////////////////////////////////////////
    int tot_nbins = (DEFAULT_NBINS)*(DEFAULT_NBINS)*(DEFAULT_NBINS);
    //int tot_nbins = (DEFAULT_NBINS+2+2)*(DEFAULT_NBINS+2+1)*(DEFAULT_NBINS+2)/6;
    //int size_hist = SUBMATRIX_SIZE*tot_nbins;
    //int size_hist = 4*tot_nbins;
    int size_hist = grid.x*tot_nbins;
    //int size_hist = tot_nbins;
    int size_hist_bytes = size_hist*sizeof(int);

    printf("Size of histogram: %d bins\t%d bytes\n",size_hist,size_hist_bytes);

    hist = (int*)malloc(size_hist_bytes);
    memset(hist, 0, size_hist_bytes);

    cudaMalloc((void **) &dev_hist, (size_hist_bytes));
    error = cudaGetLastError();
    printf("ERROR: %s\n", cudaGetErrorString(error) );

    printf("dev_hist: %x\n",dev_hist);
    cudaMemset(dev_hist, 0, size_hist_bytes);
    error = cudaGetLastError();
    printf("ERROR: %s\n", cudaGetErrorString(error) );

    printf("dev_hist bins: %d\n",size_hist);
    printf("dev_hist size: %d\n",size_hist_bytes);

    unsigned long int *summed_hist;

    unsigned long int summed_hist_size = tot_nbins * sizeof(unsigned long int);
    summed_hist =  (unsigned long int*)malloc(summed_hist_size);
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

    if ( 0==dev_hist )
    {
        printf("Couldn't allocate memory for the histogram on the GPU!\n");
        return 1;
    }

    int x, y;
    //int num_submatrices_x = NUM_GALAXIES0 / SUBMATRIX_SIZE;
    //int num_submatrices_y = NUM_GALAXIES1 / SUBMATRIX_SIZE;
    int num_submatrices[3] = {0, 0, 0};
    for (int i=0;i<3;i++)
    {
        num_submatrices[i] = NUM_GALAXIES[i]/SUBMATRIX_SIZE;
        // Take care of edges of matrix.
        if (NUM_GALAXIES[i]%SUBMATRIX_SIZE != 0)
        {
            num_submatrices[i] += 1;
        }
    }


    printf("Breaking down the calculations.\n");
    printf("Number of submatrices: %dx%dx%d\n",num_submatrices[0],num_submatrices[1],num_submatrices[2]);
    printf("Number of calculations per submatrices: %dx%d\n",SUBMATRIX_SIZE,SUBMATRIX_SIZE);

    int xind, yind, zind;
    int bin_index = 0;
    unsigned long long tot0 = 0;
    //for(int i = 0; i < NUM_GALAXIES[0]; i++)
    for(int i = 0; i < num_submatrices[0]; i++)
    {
        if (i%100==0)
        {
            printf("%d\n",i);
            fflush(stdout);
        }
        xind = i*SUBMATRIX_SIZE;
        int jmin = 0;
        if (which_three_input_files==0) // DDD or RRR
            jmin = 0;
        else if (which_three_input_files==1) // DRR or RDD
            jmin = 0;
        else if (which_three_input_files==2) // DRD or RDR
            jmin = 0;
        else if (which_three_input_files==3) // DDR or RRD
            jmin = 0;
        for(int j = jmin; j < num_submatrices[1]; j++)
        {
            yind = j*SUBMATRIX_SIZE;
            int kmin = 0;
            if (which_three_input_files==0)
                kmin = 0;
            else if (which_three_input_files==1)
                kmin = 0;
            else if (which_three_input_files==2)
                kmin = 0;
            else if (which_three_input_files==3)
                kmin = 0;
            for(int k =kmin; k < num_submatrices[2]; k++)
            {
                zind = k*SUBMATRIX_SIZE;
                //bool do_calc = 1;
                //if (do_calc)
                {

                    // Set the histogram to all zeros each time.
                    cudaMemset(dev_hist,0,size_hist_bytes);

                    int max_x = NUM_GALAXIES[0];
                    int max_y = NUM_GALAXIES[1];
                    int max_z = NUM_GALAXIES[2];

                    //printf("here\n");
                    //printf("i: %d\n",i);
                    printf("xind: %5d %5d\n",xind,max_x);
                    printf("yind: %5d %5d\n",yind,max_y);
                    printf("zind: %5d %5d\n",zind,max_z);
                    printf("nbins: %d\n",nbins);
                    //distance<<<grid,block>>>(h_x[0],h_y[0],h_z[0], 
                    distance<<<grid,block>>>(
                            d_x[0],d_y[0],d_z[0],\
                            d_x[1],d_y[1],d_z[1],\
                            d_x[2],d_y[2],d_z[2],\
                            xind, yind, zind, \
                            max_x, max_y, max_z,\
                            dev_hist, hist_lower_range, hist_upper_range, nbins, \
                            hist_bin_width, log_binning_flag, conv_factor_angle);
                    //printf("there\n");
                    //printf("dev_hist: %x\n",dev_hist);

                    error = cudaGetLastError();
                    printf("kernel ERROR: %s\n", cudaGetErrorString(error) );

                    cudaMemcpy(hist, dev_hist, size_hist_bytes, cudaMemcpyDeviceToHost);

                    cudaError_t error = cudaGetLastError();
                    printf("memcpy ERROR: %s\n", cudaGetErrorString(error) );

                    ////////////////////////////////////////////////////////////////////
                    // Sum up the histograms from each thread (hist).
                    ////////////////////////////////////////////////////////////////////
                    for(int m=0; m<size_hist; m++)
                    {
                        bin_index = m%(tot_nbins);
                        //if (hist[m]!=0)
                        //{
                        //printf("%d %lu\n",m,hist[m]);
                        //}
                        summed_hist[bin_index] += hist[m];
                        tot0 += hist[m];
                    }    
                }
            }  
        }
    }

    cudaMemcpy(hist, dev_hist, size_hist_bytes, cudaMemcpyDeviceToHost);

    unsigned long long total = 0;
    int index = 0;
    fprintf(outfile,"%d %d %d\n",nbins,nbins,nbins);
    for(int i = 0; i < nbins; i++)
    {
        //printf("%d --------------\n",i);
        fprintf(outfile,"%d\n",i);
        //for(int j = i; j < nbins; j++)
        for(int j = 0; j < nbins; j++)
        {
            //for(int k = j; k < nbins; k++)
            for(int k = 0; k < nbins; k++)
            {
                index = (nbins)*(nbins)*k + (nbins)*j + i;
                //printf("%lu ",summed_hist[index]);
                fprintf(outfile,"%lu ",summed_hist[index]);
                total += summed_hist[index];
            }
            fprintf(outfile,"\n");
            //printf("\n");
        }
    }
    printf("total: %llu\n",total);
    printf("tot0: %llu\n",tot0);

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
