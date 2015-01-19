#include<stdio.h>
#include<string.h>
#include<stdlib.h>
//#include<cmath>
#include<math.h>
#include <unistd.h>

#include <cuda_runtime.h>

using namespace std;

#define PI 3.14159
#define INVPI 1.0/PI

///////////////////////////////////////////////////////////////////////////////
// Histogram information

/*
#define S_NBINS 60
#define S_LO 0.
#define S_HI 120.

#define Q_NBINS 16
#define Q_LO 0.9
#define Q_HI 4.1

#define THETA_NBINS 25
#define THETA_LO 0.
#define THETA_HI 1.
 */

//////////////////////////
#define S_NBINS 50
#define S_LO 2.0
#define S_HI 12.0
#define S_WIDTH (S_HI-S_LO)/S_NBINS

#define Q_NBINS 16
#define Q_LO 0.9
#define Q_HI 4.1
#define Q_WIDTH (Q_HI-Q_LO)/Q_NBINS

#define THETA_NBINS 25
#define THETA_LO 0.
#define THETA_HI 1.
#define THETA_WIDTH (THETA_HI-THETA_LO)/THETA_NBINS

///////////////////////////////////////////////////////////////////////////////

#define SUBMATRIX_SIZE 4096
//#define SUBMATRIX_SIZE 1024

#define CONV_FACTOR 57.2957795 // 180/pi

float max_dist = 0.0;
float min_dist = 999.0;

///////////////////////////////////////////////////////////////////////////////
void vox2gal(int voxel_division,int voxel_index,int ngals,int *gal_indices)
{
    int ngals_in_voxel = ngals/voxel_division;
    printf("ngals: %d\n",ngals);
    printf("ngals_in_voxel: %d\n",ngals_in_voxel);
    printf("voxel_index: %d\n",voxel_index);

    if (ngals_in_voxel==0)
    {
        printf("Voxel chunks are bigger than the number of galaxies!");
        exit(-1);
    }

    int lo = ngals_in_voxel*voxel_index;
    int hi = ngals_in_voxel*(voxel_index+1);

    printf("lo/hi: %d %d\n",lo,hi);

    if (voxel_index==voxel_division-1)
    {
        hi += ngals%voxel_division; // To account for any excess, if it is not an even division.
    }

    printf("lo/hi: %d %d\n",lo,hi);

    gal_indices[0] = lo;
    gal_indices[1] = hi;

}

///////////////////////////////////////////////////////////////////////////////

__device__ int distance_to_bin(float dist, float hist_min, float hist_max, int nbins, float bin_width)
{
    int bin_index = 0;

    if(dist < hist_min)
        return 0;
        //return -999;
    else if(dist >= hist_max)
        return 0;
        //return -999;
    else
        bin_index = floor((dist-hist_min)/bin_width);

    return bin_index;
}

////////////////////////////////////////////////////////////////////////
__global__ void distance(float *x0, float *y0, float *z0, float *x1, float *y1, float *z1,float *x2, float *y2, float *z2, \
        int xind, int yind, int zind, \
        int max_xind, int max_yind, int max_zind, \
        int *dev_hist, float hist_min, float hist_max, int nbins, \
        float bin_width, float *dev_test, int flag=0)
{

    ////////////////////////////////////////////////////////////////////////////
    // Idx will keep track of which thread is being calculated within a given
    // warp.
    ////////////////////////////////////////////////////////////////////////////
    int idx = blockIdx.x * blockDim.x + threadIdx.x; // This should range to SUBMATRIX_SIZE

    int tot_hist_size = THETA_NBINS;

    int idxorg = idx;

    idx += xind;

    ////////////////////////////////////////////////////////////////////////
    // Shared memory stuff.
    ////////////////////////////////////////////////////////////////////////
    __shared__ int shared_hist[THETA_NBINS];
    // Note that we only clear things out for the first thread on each block.
    if(threadIdx.x==0)
    {
        for (int i=0;i<tot_hist_size;i++)
            shared_hist[i] = 0;
    }
    __syncthreads();
    ////////////////////////////////////////////////////////////////////////

    float x0i = 0.0;
    float y0i = 0.0;
    float z0i = 0.0;

    float x;

    int i=0;
    int j=0;
    int k=0;

    float xpt0,ypt0,zpt0;
    float xpt1,ypt1,zpt1;
    float xpt2,ypt2,zpt2;

    float xdiff,ydiff,zdiff;
    float dist0,dist1,dist2;

    //bool do_calc = 1;

    int ymax = yind + SUBMATRIX_SIZE;
    int zmax = zind + SUBMATRIX_SIZE;

    ymax = (ymax>max_yind ? max_yind : ymax);
    zmax = (zmax>max_zind ? max_zind : zmax);

    int totbin = 0;
    float s,q,theta0,theta1,theta2;
    int i0=-1; // shortest
    int i1=-1; // middle
    int i2=-1; // longest

    float shortest,middle,longest;

    float shortest2;
    float middle2;
    float longest2;


    if (idx<max_xind)
    {
        xpt0 = x0[idx];
        ypt0 = y0[idx];
        zpt0 = z0[idx];


        //for(j=yind; j<ymax; j++)
        int jmin = 0;
        if (flag==0) // DDD or RRR
            jmin = idx;
        else if (flag==1) // DRR or RDD
            jmin = 0;
        else if (flag==2) // DRD or RDR
            jmin = 0;
        else if (flag==3) // DDR or RRD
            jmin = idx;
        for(j=jmin; j<ymax; j++)
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
            //for(k=zind; k<zmax; k++)
            int kmin = 0;
            if (flag==0)
                kmin = j;
            else if (flag==1)
                kmin = j;
            else if (flag==2)
                kmin = i;
            else if (flag==3)
                kmin = 0;
            for(k=kmin; k<zmax; k++)
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


                i0=-1; // shortest
                i1=-1; // middle
                i2=-1; // longest

                totbin=-999;

                // From Stackoverflow 
                // http://stackoverflow.com/questions/13040240/the-fastest-way-to-sort-3-values-java
                // This is if we want to sort all 3.
                if( dist0 > dist1 ){
                    if( dist0 > dist2 ){
                        longest = dist0;
                        if( dist1 > dist2 ){
                            middle = dist1;
                            shortest = dist2;
                        }else{
                            middle = dist2;
                            shortest = dist1;
                        }
                    }else{
                        middle = dist0;
                        if( dist1 > dist2 ){
                            longest = dist1;
                            shortest = dist2;
                        }else{
                            longest = dist2;
                            shortest = dist1;
                        }
                    }
                }else{
                    if( dist1 > dist2 ){
                        longest = dist1;
                        if( dist0 > dist2 ){
                            middle = dist0;
                            shortest = dist2;
                        }else{
                            middle = dist2;
                            shortest = dist0;
                        }
                    }else{
                        middle = dist1;
                        longest = dist2;
                        shortest = dist0;
                    }
                }

                //shortest = 1.0;
                //middle = 1.0;
                //longest = 1.5;

                //shortest = dist0;
                //middle = dist1;
                //longest = dist2;

                shortest2 = shortest*shortest;
                middle2 = middle*middle;
                longest2 = longest*longest;

                //dev_test[0] = xpt0;
                //dev_test[1] = xpt1;
                //dev_test[2] = xpt2;

                //if (idx==10 && j==20 && k==90)
                //{
                    //dev_test[0] = dist0;
                    //dev_test[1] = dist1;
                    //dev_test[2] = dist2;
                //}

                //dev_test[0] = xpt1;
                //dev_test[1] = xpt2;
                //dev_test[2] = xpt0;

                for (int n=0;n<3;n++)
                {
                    i0=0; 
                    i1=0; 
                    i2=0; 
                    if (n==0) {
                        s = shortest;
                        q = middle/shortest;
                        theta0 = (acosf((shortest2 + middle2 - longest2)/(2*shortest*middle)))*INVPI;
                        //theta0 = 0.5;
                        i2 = distance_to_bin(theta0,THETA_LO,THETA_HI,THETA_NBINS,THETA_WIDTH);
                    } else if (n==1){
                        s = middle;
                        q = longest/middle;
                        theta1 = (acosf((middle2 + longest2 - shortest2)/(2*middle*longest)))*INVPI;
                        i2 = distance_to_bin(theta1,THETA_LO,THETA_HI,THETA_NBINS,THETA_WIDTH);
                    } else if (n==2){
                        s = shortest;
                        q = longest/shortest;
                        //theta2 = (acosf((shortest2 + longest2 - middle2)/(2*shortest*longest)))*INVPI;
                        theta2 = 1.0 - theta0 - theta1;
                        //i2 = distance_to_bin(0.5,THETA_LO,THETA_HI,THETA_NBINS,THETA_WIDTH);
                        i2 = distance_to_bin(theta2,THETA_LO,THETA_HI,THETA_NBINS,THETA_WIDTH);
                        //i2 = 13;
                    }

                    //printf("%f %f %f\n",s,q,theta);

                    i0 = distance_to_bin(s,S_LO,S_HI,S_NBINS,S_WIDTH); //Mpc/h, delta s=0.2
                    i1 = distance_to_bin(q,Q_LO,Q_HI,Q_NBINS,Q_WIDTH); // delta q = 0.2
                    //i2 = distance_to_bin(theta,0,1.0,25,THETA_WIDTH);

                    //if (idx==1000 && j==1001 && k==1002)
                    /*
                    if ((dist0<12 && dist1<12 && dist2<12 &&
                         dist0>2 && dist1>2 && dist2>2))
                    {
                        dev_test[0] = float(i0);
                        dev_test[1] = float(i1);
                        dev_test[2] = float(i2);
                        dev_test[3] = dist0;
                        dev_test[4] = dist1;
                        dev_test[5] = dist2;
                        //dev_test[0] = s;
                        //dev_test[1] = q;
                        //dev_test[2] = theta2;
                    }
                    */
                    // Increment the ``column" of the histogram
                    /*
                    if (i0==49)
                        i0 = 1;
                    if (i1==15)
                        i1 = 1;
                    if (i2<0)
                        i2 = 0;

                    totbin = i0*i1*i2;
                    */

                    //if (i2>=0)
                    //if (1)
                    //if(totbin>0 && totbin<THETA_NBINS)
                    if (i0==49 && i1==15 && i2>=0)
                    {
                        atomicAdd(&shared_hist[i2],1);
                        //shared_hist[i2] +=1;
                        //shared_hist[totbin] +=1;
                        //atomicAdd(&shared_hist[totbin],1);
                    }

                }
            }
        }
    }

    //__syncthreads();

    if(threadIdx.x==0)
    {
        for(int i=0;i<tot_hist_size;i++)
        {
            dev_hist[i+(blockIdx.x*tot_hist_size)]=shared_hist[i];
            //dev_hist[i] = 20;
        }
    }


}

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    // Needed for parsing command-line arguments.
    extern char *optarg;
    extern int optind, optopt, opterr;
    int c;
    char *filename;
    char *outfilename = NULL;
    char defaultoutfilename[256];
    sprintf(defaultoutfilename,"default_out.dat");

    char histout[1024];
    sprintf(histout,"0");

    int log_binning_flag = 0; // False

    int flag = 0;
    int voxel_division = -999;
    int voxel_index[3] = {0,0,0};
    int vtemp = -999;

    cudaError_t error = cudaGetLastError();

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    while ((c = getopt(argc, argv, "ao:l:x:X:")) != -1) {
        switch(c) {
            case 'l':
                log_binning_flag = atoi(optarg);
                break;
            case 'o':
                outfilename = optarg;
                printf("Output filename is %s\n", outfilename);
                break;
            case 'X':
                voxel_division = atoi(optarg);
                printf("Voxel division is %d\n", voxel_division);
                break;
            case 'x':
                vtemp = atoi(optarg);
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

    if(voxel_division>0)
    {
        voxel_index[0] = vtemp/1000000;
        voxel_index[1] = (vtemp/1000)%1000;
        voxel_index[2] = (vtemp%1000);
        printf("vtemp: %d\n",vtemp);
        printf("Voxel indices are %03d %03d %03d\n",voxel_index[0],voxel_index[1],voxel_index[2]);
    }

    printf("Log binning flag: %d\n",log_binning_flag);

    float *htemp_x[3], *htemp_y[3], *htemp_z[3];
    float *h_x[3], *h_y[3], *h_z[3];
    float *d_x[3], *d_y[3], *d_z[3];

    // Open the input files and the output file.
    FILE *infile[3], *outfile;
    for (int i=0;i<3;i++)
    {
        infile[i] = fopen(argv[optind+i],"r");
        printf("Opening input file %d: %s\n",i,argv[optind+i]);
    }

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

    ////////////////////////////////////////////////////////////////////////////
    // Read in the files.
    ////////////////////////////////////////////////////////////////////////////

    int NUM_GALAXIES[3] = {0,0,0};
    int size_of_galaxy_array[3];
    int idummy;
    float temp0, temp1, temp2, dummy;

    int max_ngals = 1000000;

    for (int i=0;i<3;i++)
    {
        //fscanf(infile[i], "%d", &NUM_GALAXIES[i]);
        //size_of_galaxy_array[i] = NUM_GALAXIES[i] * sizeof(float);    
        //printf("SIZE %d # GALAXIES: %d\n",i,NUM_GALAXIES[i]);

        htemp_x[i] = (float*)malloc(max_ngals*sizeof(float));
        htemp_y[i] = (float*)malloc(max_ngals*sizeof(float));
        htemp_z[i] = (float*)malloc(max_ngals*sizeof(float));

        int j = 0;
        while(fscanf(infile[i], "%d %f %f %f %f %f %f", &idummy, &dummy, &dummy, &dummy, &temp0, &temp1, &temp2) != EOF)
        {
            htemp_x[i][j] = temp0;///scale_factor;
            htemp_y[i][j] = temp1;///scale_factor;
            htemp_z[i][j] = temp2;///scale_factor;
            if (NUM_GALAXIES[i]>=max_ngals)
            {
                printf("Exceeded max num galaxies");
                exit(-1);
            }
            //if (j<10)
            //{
                //printf("%f %f %f\n", htemp_x[i][j],htemp_y[i][j],htemp_z[i][j]);
            //}
            NUM_GALAXIES[i] += 1;
            j += 1;
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Allocation the arrays of galaxies that we actually want to run over.
    ///////////////////////////////////////////////////////////////////////////

    //int min_gals[3] = {0,0,0};
    //int ngals[3] = {500,500,500};
    //int max_gals[3] = {min_gals[0]+ngals[0],min_gals[1]+ngals[1],min_gals[2]+ngals[2]};
    int min_gals[3] = {0,0,0};
    int max_gals[3] = {NUM_GALAXIES[0],NUM_GALAXIES[1],NUM_GALAXIES[2]};
    int ngals[3] = {NUM_GALAXIES[0],NUM_GALAXIES[1],NUM_GALAXIES[2]};
    int lohi[2] = {0,0};

    if (voxel_division>0 && vtemp>=0)
    {
        for (int i=0;i<3;i++)
        {
            vox2gal(voxel_division,voxel_index[i],NUM_GALAXIES[i],lohi);
            min_gals[i] = lohi[0];
            max_gals[i] = lohi[1];
            ngals[i] = max_gals[i]-min_gals[i];
        }
    }

    for (int i=0;i<3;i++)
    {
        size_of_galaxy_array[i] = ngals[i] * sizeof(float);

        printf("Galaxy indices: %d %d %d\n",min_gals[i],max_gals[i],ngals[i]);
        h_x[i] = (float*)malloc(ngals[i]*sizeof(float));
        h_y[i] = (float*)malloc(ngals[i]*sizeof(float));
        h_z[i] = (float*)malloc(ngals[i]*sizeof(float));

        int index = 0;
        for(int j=min_gals[i];j<max_gals[i];j++)
        {
            h_x[i][index] = htemp_x[i][j];
            h_y[i][index] = htemp_y[i][j];
            h_z[i][index] = htemp_z[i][j];

            //if (index<10 || index>570)
            //{
            //printf("%d %f %f %f\n",index,h_x[i][index],h_y[i][index],h_z[i][index]);
            //}
            index++;
        }
    }

    printf("Finished filling the real galaxy arrays....\n");


    ////////////////////////////////////////////////////////////////////////////
    // Define the grid and block size
    ////////////////////////////////////////////////////////////////////////////
    dim3 grid, block;
    // 128*4 = 512, the amount of memory needed for one histogram.
    // 8192*4 = 32768 is max memory to ask for for the histograms.
    // 8192/128 = 64, is is the right number of blocks?
    //grid.x = 8192/(tot_nbins); // Is this the number of blocks?
    grid.x = 32; // Is this the number of blocks?
    //grid.x = 8; // Is this the number of blocks?
    //grid.x = 4; // Is this the number of blocks?
    block.x = SUBMATRIX_SIZE/grid.x; // Is this the number of threads per block? NUM_GALAXIES/block.x;
    //block.x = SUBMATRIX_SIZE; // Is this the number of threads per block? NUM_GALAXIES/block.x;
    // SUBMATRIX is the number of threads per warp? Per kernel call?
    printf("# of blocks per grid:  grid.x:  %d\n",grid.x);
    printf("# of thread per block: block.x: %d\n",block.x);
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // Allocation of histograms.
    ///////////////////////////////////////////////////////////////////////////
    int *hist;
    int *dev_hist;
    float *dev_test;
    float *host_test;
    unsigned long int *summed_hist;
    //int nbins;
    int log_binning=flag;

    int size_hist = grid.x*(THETA_NBINS);
    int size_hist_bytes = size_hist*sizeof(int);

    host_test = (float*)malloc(3*sizeof(float));

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

    cudaMalloc((void **) &dev_test, (6*sizeof(float)));
    cudaMemset(dev_test, 0, (6*sizeof(float)));


    cudaMemcpy(dev_hist, hist, size_hist_bytes, cudaMemcpyHostToDevice);
    error = cudaGetLastError();
    printf("BEFORE KERNEL memcpy copying hist to dev_hist ERROR: %s\n", cudaGetErrorString(error) );

    cudaMemcpy(hist, dev_hist, size_hist_bytes, cudaMemcpyDeviceToHost);
    error = cudaGetLastError();
    printf("BEFORE KERNEL memcpy copying dev_hist to hist ERROR: %s\n", cudaGetErrorString(error) );


    int summed_hist_size = THETA_NBINS * sizeof(unsigned long int);
    summed_hist =  (unsigned long int*)malloc(summed_hist_size);
    printf("Size of summed histogram array: %d bins\t%d bytes\n",(THETA_NBINS),summed_hist_size);
    memset(summed_hist,0,summed_hist_size);




    ////////////////////////////////////////////////////////////////////////////
    // Allocate the memory on the GPU
    ////////////////////////////////////////////////////////////////////////////
    // Check to see if we allocated enough memory.
    for (int i=0;i<3;i++)
    {
        cudaMalloc((void **) &d_x[i], size_of_galaxy_array[i] );
        cudaMalloc((void **) &d_y[i], size_of_galaxy_array[i] );
        cudaMalloc((void **) &d_z[i], size_of_galaxy_array[i] );

        printf("Size: %d\n",size_of_galaxy_array[i]);
        error = cudaGetLastError();
        printf("ERROR cudaMalloc coordinates: %s\n", cudaGetErrorString(error) );

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

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    int x, y;
    float dist = 0;

    bool locked = false;
    int bin_index = 0;
    int bins[3] = {0,0,0};

    int min_index[3] = {0,0,0};
    //int max_index[3] = {NUM_GALAXIES[0],NUM_GALAXIES[1],NUM_GALAXIES[2]};
    int max_index[3] = {ngals[0],ngals[1],ngals[2]};

    unsigned long long numcalc = 0;

    //int num_submatrices_x = NUM_GALAXIES0 / SUBMATRIX_SIZE;
    //int num_submatrices_y = NUM_GALAXIES1 / SUBMATRIX_SIZE;
    int num_submatrices[3] = {0, 0, 0};
    for (int i=0;i<3;i++)
    {
        num_submatrices[i] = ngals[i]/SUBMATRIX_SIZE;
        // Take care of edges of matrix.
        if (ngals[i]%SUBMATRIX_SIZE != 0)
        {
            num_submatrices[i] += 1;
        }
    }


    ////////////////////////////////////////////////////////////////////////////
    // Start the calculations
    ////////////////////////////////////////////////////////////////////////////
    printf("About to enter the loops...\n");
    printf("Breaking down the calculations.\n");
    printf("Number of submatrices: %dx%dx%d\n",num_submatrices[0],num_submatrices[1],num_submatrices[2]);
    printf("Number of calculations per submatrices: %dx%d\n",SUBMATRIX_SIZE,SUBMATRIX_SIZE);

    int xind, yind, zind;
    unsigned long long tot0 = 0;

    float hist_bin_width = (THETA_HI-THETA_LO)/THETA_NBINS;

    //for(int i=min_index[0];i<max_index[0]; i++)
    for(int i = 0; i < num_submatrices[0]; i++)
    {
        if (i%100==0)
        {
            printf("%d\n",i);
            fflush(stdout); 
        }
        //int jmin = min_index[1];
        xind = i*SUBMATRIX_SIZE;
        int jmin = 0;
        //if (which_three_input_files==0) // DDD or RRR
            //jmin = i;
        //else if (which_three_input_files==1) // DRR or RDD
            //jmin = 0;
        //else if (which_three_input_files==2) // DRD or RDR
            //jmin = 0;
        //else if (which_three_input_files==3) // DDR or RRD
            //jmin = i;
        //for(int j=jmin;j<max_index[1];j++)
        //for(int j = 0; j < NUM_GALAXIES[1]; j++)
        for(int j = jmin; j < num_submatrices[1]; j++)
        {
            yind = j*SUBMATRIX_SIZE;

            //int kmin = min_index[2];
            int kmin = 0;
            //if (which_three_input_files==0)
                //kmin = j;
            //else if (which_three_input_files==1)
                //kmin = j;
            //else if (which_three_input_files==2)
                //kmin = i;
            //else if (which_three_input_files==3)
                //kmin = 0;
            //for(int k=kmin;k<max_index[2];k++)
            //for(int k =0; k < ngals[2]; k++)
            for(int k =kmin; k < num_submatrices[2]; k++)
            {
                zind = k*SUBMATRIX_SIZE;


                int max_x = ngals[0];
                int max_y = ngals[1];
                int max_z = ngals[2];

                cudaMemset(dev_hist,0,size_hist_bytes);

                printf("xind: %5d %5d\n",xind,max_x);
                printf("yind: %5d %5d\n",yind,max_y);
                printf("zind: %5d %5d\n",zind,max_z);

                //printf("%d ",k);
                distance<<<grid,block>>>(
                        d_x[0],d_y[0],d_z[0], \
                        d_x[1],d_y[1],d_z[1], \
                        d_x[2],d_y[2],d_z[2], \
                        xind, yind, zind, \
                        max_x, max_y, max_z,\
                        dev_hist, THETA_LO, THETA_HI, THETA_NBINS, \
                        hist_bin_width, dev_test, which_three_input_files);

                //printf("bin_index: %d\n",bin_index);
                numcalc += 1;


                error = cudaGetLastError();
                printf("kernel ERROR: %s\n", cudaGetErrorString(error) );

                cudaMemcpy(hist, dev_hist, size_hist_bytes, cudaMemcpyDeviceToHost);

                error = cudaGetLastError();
                printf("memcpy copying hist to dev_hist ERROR: %s\n", cudaGetErrorString(error) );

                cudaMemcpy(host_test, dev_test, 6*sizeof(float), cudaMemcpyDeviceToHost);

                printf("host_test: %f %f %f\n",host_test[0],host_test[1],host_test[2]);
                printf("host_test: %f %f %f\n",host_test[3],host_test[4],host_test[5]);

                for(int m=0; m<size_hist; m++)
                {
                    bin_index = m%(THETA_NBINS);
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

    //exit(0);

    int index = 0;
    unsigned long long total = 0;
    outfile = fopen(outfilename, "w");
    fprintf(outfile,"%d\n",NUM_GALAXIES[0]);
    fprintf(outfile,"%d\n",NUM_GALAXIES[1]);
    fprintf(outfile,"%d\n",NUM_GALAXIES[2]);
    fprintf(outfile,"%d %d %d\n",S_NBINS,Q_NBINS,THETA_NBINS);
    fprintf(outfile,"%-6.3f %-6.3f %-6.3f\n",S_LO,S_HI,(S_HI-S_LO)/S_NBINS);
    fprintf(outfile,"%-6.3f %-6.3f %-6.3f\n",Q_LO,Q_HI,(Q_HI-Q_LO)/Q_NBINS);
    fprintf(outfile,"%-6.3f %-6.3f %-6.3f\n",THETA_LO,THETA_HI,(THETA_HI-THETA_LO)/THETA_NBINS);
    //for(int i = 0; i < S_NBINS; i++)
    {
        //printf("%d --------------\n",i);
        //fprintf(outfile,"%d\n",i);
        //for(int j = 0; j < Q_NBINS; j++)
        {
            for(int k = 0; k < THETA_NBINS; k++)
            {

                //index = (Q_NBINS)*(THETA_NBINS)*i + (THETA_NBINS)*j + k; 
                index = k;
                printf("%lu ",summed_hist[index]);
                fprintf(outfile,"%lu ",summed_hist[index]);
                total += summed_hist[index];
            }
            fprintf(outfile,"\n");
            printf("\n");
        }
    }

    printf("Total: %lu\n",total);
    printf("numcalc: %lu\n",numcalc);

    fclose(outfile);
    return 0;
}  
//////////////////////////////////////////////////////////////////////
