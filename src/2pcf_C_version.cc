#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<cmath>
#include <unistd.h>

using namespace std;

#define PI 3.14159

///////////////////////////////////////////////////////////////////////////////
// Histogram information
//////////////////////////
#define S_NBINS 700
#define S_LO 0.0
#define S_HI 70.0
///////////////////////////////////////////////////////////////////////////////

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

int distance_to_bin(float dist, float hist_min, float hist_max, int nbins, int flag)
{
    int bin_index = 0;
    float bin_width=(hist_max-hist_min)/nbins;

    if(dist < hist_min)
        return -999;
    else if(dist >= hist_max)
        return -999;
    else
    {
        if (flag==0)
        {
            //printf("%f %f %f\n",dist,hist_min,bin_width);
            //bin_index = int((dist-hist_min)/bin_width) + 1;
            bin_index = int((dist-hist_min)/bin_width);
            //printf("here! %d\n",bin_index);
        }
        else if (flag==1)// log binning
        {
            bin_index = int((log(dist)-log(hist_min))/bin_width);
        }
        else if (flag==2)// log 10 binning
        {
            bin_index = int((log10(dist)-log10(hist_min))/bin_width);
        }
    }

    return bin_index;
}

////////////////////////////////////////////////////////////////////////
int distance(float x0, float y0, float z0, float x1, float y1, float z1, int flag, int *totbins)
{

    float xdiff = x0-x1;
    float ydiff = y0-y1;
    float zdiff = z0-z1;
    float dist0 = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

    //int totbin = -999;

        int i0 = distance_to_bin(dist0,S_LO,S_HI,S_NBINS,flag); //Mpc/h, delta s=0.2

        //printf("%d %d %d\n",i0,i1,i2);
        if (i0<0)
        {
            *totbins = -999;
        } else {
            // Combine for big 1d rep of 3d histogram;
            //int nhistbins = nbins;
            //int nhistbins2 = nhistbins*nhistbins;
            *totbins = i0;
        }

    return 1;
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
    sprintf(defaultoutfilename,"default_out_2pcf.dat");

    char histout[1024];
    sprintf(histout,"0");

    int log_binning_flag = 0; // False

    int flag = 0;
    int voxel_division = -999;
    int voxel_index[2] = {0,0};
    int vtemp = -999;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    //while ((c = getopt(argc, argv, "ao:L:l:w:smx:X:")) != -1) {
    while ((c = getopt(argc, argv, "ao:l:x:X:")) != -1) {
        switch(c) {
            //case 'L':
                //printf("L is set\n");
                //hist_lower_range = atof(optarg);
                //break;
            //case 'w':
                //hist_bin_width = atof(optarg);
                //printf("Histogram bin width: %f\n",hist_bin_width);
                //break;
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

    if (argc < 2)
    {
        printf("\nMust pass in at least two input files on command line!\n");
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
        printf("vtemp: %d\n",vtemp);
        printf("Voxel indices are %03d %03d\n",voxel_index[0],voxel_index[1]);
    }

    printf("Log binning flag: %d\n",log_binning_flag);

    float *htemp_x[2], *htemp_y[2], *htemp_z[2];
    float *h_x[2], *h_y[2], *h_z[2];

    // Open the input files and the output file.
    FILE *infile[2], *outfile;
    for (int i=0;i<2;i++)
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
    if (strcmp(argv[optind+0],argv[optind+1])==0)
    {
        which_three_input_files = 0;
        printf("Using the same file! (DD or RR)\n");
    }
    else if (strcmp(argv[optind+0],argv[optind+1])!=0)
    {
        which_three_input_files = 1;
        printf("Not the same file! (DR or RD)\n");
    }

    ////////////////////////////////////////////////////////////////////////////
    // Read in the files.
    ////////////////////////////////////////////////////////////////////////////

    int NUM_GALAXIES[2] = {0,0};
    int size_of_galaxy_array[2];
    int idummy;
    float temp0, temp1, temp2, dummy;

    int max_ngals = 1000000;

    for (int i=0;i<2;i++)
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
            ///*
            //if (j<10)
            //{
                //printf("%f %f %f\n", htemp_x[i][j],htemp_y[i][j],htemp_z[i][j]);
            //}
            //*/
            NUM_GALAXIES[i] += 1;
            j += 1;
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Allocation the arrays of galaxies that we actually want to run over.
    ///////////////////////////////////////////////////////////////////////////

    //int min_gals[2] = {0,0,0};
    //int ngals[2] = {500,500,500};
    //int max_gals[2] = {min_gals[0]+ngals[0],min_gals[1]+ngals[1],min_gals[2]+ngals[2]};
    int min_gals[2] = {0,0};
    int max_gals[2] = {NUM_GALAXIES[0],NUM_GALAXIES[1]};
    int ngals[2] = {NUM_GALAXIES[0],NUM_GALAXIES[1]};
    int lohi[2] = {0,0};

    if (voxel_division>0 && vtemp>=0)
    {
        for (int i=0;i<2;i++)
        {
            vox2gal(voxel_division,voxel_index[i],NUM_GALAXIES[i],lohi);
            min_gals[i] = lohi[0];
            max_gals[i] = lohi[1];
            ngals[i] = max_gals[i]-min_gals[i];
        }
    }

    for (int i=0;i<2;i++)
    {
        printf("Galaxy indices: %d %d %d\n",min_gals[i],max_gals[i],ngals[i]);
        h_x[i] = (float*)malloc(ngals[i]*sizeof(float));
        h_y[i] = (float*)malloc(ngals[i]*sizeof(float));
        h_z[i] = (float*)malloc(ngals[i]*sizeof(float));

        //printf("ngals[%d] %d\n",i,ngals[i]);
        //h_x[i] = (float*)malloc(1000);
        //h_y[i] = (float*)malloc(1000);
        //h_z[i] = (float*)malloc(1000);

        int index = 0;
        for(int j=min_gals[i];j<max_gals[i];j++)
        {
            h_x[i][index] = htemp_x[i][j];
            h_y[i][index] = htemp_y[i][j];
            h_z[i][index] = htemp_z[i][j];

            //printf("indices: %d %d\n",i,index);
            //h_x[i][index] = 0.0;
            //h_y[i][index] = 0.0;
            //h_z[i][index] = 0.0;

            //if (index<10 || index>570)
            //{
                //printf("%d %f %f %f\n",index,h_x[i][index],h_y[i][index],h_z[i][index]);
            //}
            index++;
        }
    }

    printf("Finished filling the real galaxy arrays....\n");

    ////////////////////////////////////////////////////////////////////////////
    // Allocation of histograms.
    ///////////////////////////////////////////////////////////////////////////

    unsigned int *hist;
    unsigned int *temp_hist;
    //int nbins;
    int log_binning=flag;

    int size_hist = (S_NBINS);
    int size_hist_bytes = size_hist*sizeof(unsigned int);

    hist = (unsigned int*)malloc(size_hist_bytes);
    memset(hist, 0, size_hist_bytes);

    temp_hist = (unsigned int*)malloc(size_hist_bytes);
    memset(temp_hist, 0, size_hist_bytes);

    int x, y;
    float dist = 0;

    bool locked = false;
    int bin_index = 0;
    int calc_count = 0;
    int calc_count_max = 100;
    int bins = 0;

    int min_index[2] = {0,0};
    int max_index[2] = {ngals[0],ngals[1]};

    printf("About to enter the loops...\n");

    for(int i=min_index[0];i<max_index[0]; i++)
    {
        if (i%100==0)
        {
            printf("%d\n",i);
            fflush(stdout); 
        }
        int jmin = min_index[1];
        if (which_three_input_files==0) // DDD or RRR
            //jmin = i+1;
            jmin = i;
        else if (which_three_input_files==1) // DRR or RDD
            jmin = 0;
        for(int j=jmin;j<max_index[1];j++)
            //for(int j = 0; j < NUM_GALAXIES[1]; j++)
        {
                    bin_index = distance(h_x[0][i],h_y[0][i],h_z[0][i], \
                            h_x[1][j],h_y[1][j],h_z[1][j], \
                            flag, &bins);

                    //printf("bin_index: %d\n",bin_index);
                        bin_index = bins;
                        if (bin_index>=0)
                        {
                            hist[bin_index]++;
                            temp_hist[bin_index]++;
                            calc_count += 1;
                        }
        }
    }  

    //exit(0);

    int index = 0;
    unsigned long long total = 0;
    outfile = fopen(outfilename, "w");
    fprintf(outfile,"%d\n",NUM_GALAXIES[0]);
    fprintf(outfile,"%d\n",NUM_GALAXIES[1]);
    fprintf(outfile,"%d\n",S_NBINS);
    fprintf(outfile,"%-6.3f %-6.3f %-6.3f\n",S_LO,S_HI,(S_HI-S_LO)/S_NBINS);
    for(int i = 0; i < S_NBINS; i++)
    {
                index = i;
                printf("%lu\n",hist[index]);
                fprintf(outfile,"%lu\n",hist[index]);
                total += hist[index];
    }

    printf("Total: %lu\n",total);

    fclose(outfile);
    return 0;
}  
//////////////////////////////////////////////////////////////////////
