#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<cmath>
#include <unistd.h>

using namespace std;

#define HIST_MIN 0.0 // for degrees
#define HIST_MAX 100.0 // for degrees

//#define DEFAULT_NBINS 20 // for log binning
#define DEFAULT_NBINS 64 // for log binning
//#define DEFAULT_NBINS 16 // for log binning
//#define DEFAULT_NBINS 4 // for log binning
//#define DEFAULT_NBINS 126 // for log binning
//#define DEFAULT_NBINS 62 // for log binning

#define CONV_FACTOR 57.2957795 // 180/pi

float max_dist = 0.0;

int distance_to_bin(float dist, float hist_min, float hist_max, int nbins, float bin_width, int flag)
{
    int bin_index = 0;

    if(dist < hist_min)
        bin_index = 0;
    else if(dist >= hist_max)
        bin_index = nbins + 1;
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
            bin_index = int((log(dist)-log(hist_min))/bin_width) + 1;
        }
        else if (flag==2)// log 10 binning
        {
            bin_index = int((log10(dist)-log10(hist_min))/bin_width) + 1;
        }
    }

    return bin_index;
}

////////////////////////////////////////////////////////////////////////
int distance(float x0, float y0, float z0, float x1, float y1, float z1,float x2, float y2, float z2, float hist_min, float hist_max, int nbins, float bin_width, int flag)
{

    float xdiff = x0-x1;
    float ydiff = y0-y1;
    float zdiff = z0-z1;
    float dist0 = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

    xdiff = x0-x2;
    ydiff = y0-y2;
    zdiff = z0-z2;
    float dist1 = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

    xdiff = x1-x2;
    ydiff = y1-y2;
    zdiff = z1-z2;
    float dist2 = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

    if (dist0>max_dist) { max_dist = dist0; printf("max_dist: %f\n",max_dist); }
    if (dist1>max_dist) { max_dist = dist1; printf("max_dist: %f\n",max_dist); }
    if (dist2>max_dist) { max_dist = dist2; printf("max_dist: %f\n",max_dist); }

    //printf("---------\n");
    //printf("%f %f %f\n",x0,x1,x2);
    //printf("%f %f %f\n",y0,y1,y2);
    //printf("%f %f %f\n",z0,z1,z2);
    //printf("%f %f %f\n",dist0,dist1,dist2);

    // Sort the distances
    bool b0 = dist0<dist1;
    bool b1 = dist1<dist2;
    bool b2 = dist0<dist2;

    int i0=-1; // shortest
    int i1=-1; // middle
    int i2=-1; // longest

    i0 = distance_to_bin(dist0,hist_min,hist_max,nbins,bin_width,flag);
    i1 = distance_to_bin(dist1,hist_min,hist_max,nbins,bin_width,flag);
    i2 = distance_to_bin(dist2,hist_min,hist_max,nbins,bin_width,flag);

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
    
    //printf("%d %d %d\n",i0,i1,i2);

    // Combine for big 1d rep of 3d histogram;
    int nhistbins = nbins;
    int nhistbins2 = nhistbins*nhistbins;
    int totbin = nhistbins2*i2 + nhistbins*i1 + i0;

    return totbin;
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

    float hist_lower_range = 0.0000001;
    float hist_upper_range = 0;
    float hist_bin_width = 0.05;
    int log_binning_flag = 0; // False

    float scale_factor = 1.0; // For if we need to convert input to arcsec or arcmin
    float conv_factor_angle = 57.2957795; // 180/pi // For if we need to convert arcdistance to arcsec or arcmin

    float hist_min = 0;
    //float hist_max = 1.8;
    //float hist_max = 7000.0;
    float hist_max = sqrt(3*24*24);
    int nbins = DEFAULT_NBINS;
    float bin_width = (hist_max-hist_min)/nbins;
    int flag = 0;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////


    while ((c = getopt(argc, argv, "ao:L:l:w:sm")) != -1) {
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
                printf("Will use log binning.\n");
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

    float temp_lo = hist_lower_range;
    if (hist_upper_range == 0)
    {
        if (log_binning_flag==0)
        {
            for (int i=0;i<nbins;i++)
            {
                hist_upper_range = temp_lo + hist_bin_width;
                temp_lo = hist_upper_range;
            }
        }
        else if (log_binning_flag==1)
        {
            for (int i=0;i<nbins;i++)
            {
                hist_upper_range = exp(log(temp_lo) + hist_bin_width);
                temp_lo = hist_upper_range;
            }
        }
        else if (log_binning_flag==2)
        {
            for (int i=0;i<nbins;i++)
            {
                hist_upper_range = pow(10,(log10(temp_lo) + hist_bin_width));
                temp_lo = hist_upper_range;
            }
        }
    }
    printf("hist_upper_range: %f\n",hist_upper_range);

    float *h_x[3], *h_y[3], *h_z[3];

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

    ////////////////////////////////////////////////////////////////////////////
    // Read in the files.
    ////////////////////////////////////////////////////////////////////////////

    int NUM_GALAXIES[3];
    int size_of_galaxy_array[3];
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
            /*
            if (i<10)
                printf("%e %e %e\n", h_x[i][j],h_y[i][j],h_z[i][j]);
            */
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Allocation of histograms.
    ///////////////////////////////////////////////////////////////////////////

    unsigned long long *hist;
    unsigned long long *temp_hist;
    //int nbins;
    int log_binning=flag;

    int size_hist = (nbins)*(nbins)*(nbins);
    int size_hist_bytes = size_hist*sizeof(unsigned long long);

    hist = (unsigned long long*)malloc(size_hist_bytes);
    memset(hist, 0, size_hist_bytes);

    temp_hist = (unsigned long long*)malloc(size_hist_bytes);
    memset(temp_hist, 0, size_hist_bytes);

    int x, y;
    float dist = 0;

    unsigned long long int fake_tot = 0;

    bool locked = false;
    int bin_index = 0;
    int calc_count = 0;
    int calc_count_max = 100;
    int num_locked = 0;
    int num_not_locked = 0;
    for(int i = 0; i < NUM_GALAXIES[0]; i++)
    {
        if (i%100==0)
        {
            printf("%d\n",i);
            fflush(stdout); 
        }
        int jmin = 0;
        if (which_three_input_files==0) // DDD or RRR
            jmin = i+1;
        else if (which_three_input_files==1) // DRR or RDD
            jmin = 0;
        else if (which_three_input_files==2) // DRD or RDR
            jmin = 0;
        else if (which_three_input_files==3) // DDR or RRD
            jmin = i+1;
        //for(int j = jmin; j < NUM_GALAXIES[1]; j++)
        for(int j = 0; j < NUM_GALAXIES[1]; j++)
        {
            int kmin = 0;
            if (which_three_input_files==0)
                kmin = j+1;
            else if (which_three_input_files==1)
                kmin = j+1;
            else if (which_three_input_files==2)
                kmin = i+1;
            else if (which_three_input_files==3)
                kmin = 0;
            //for(int k =kmin; k < NUM_GALAXIES[2]; k++)
            for(int k =0; k < NUM_GALAXIES[2]; k++)
            {
                //bool do_calc = 1;
                //if (do_calc)
                {
                    bin_index = distance(h_x[0][i],h_y[0][i],h_z[0][i], \
                                         h_x[1][j],h_y[1][j],h_z[1][j], \
                                         h_x[2][k],h_y[2][k],h_z[2][k], \
                                         hist_min, hist_max, nbins, bin_width, flag);

                    //printf("%d\n",bin_index);
                    hist[bin_index]++;
                    temp_hist[bin_index]++;
                    calc_count += 1;
                    if (calc_count>calc_count_max)
                    {
                        for(int i=0;i<size_hist;i++)
                        {
                           if(temp_hist[i]>1)
                           {
                               printf("%d %d %d   ",calc_count,i,temp_hist[i]);
                               locked = true;
                           }
                        }
                       printf("--------\n");

                        if (locked)
                            num_locked++;
                        else
                            num_not_locked++;

                        memset(temp_hist, 0, size_hist_bytes);
                        calc_count = 0;
                        locked = false;
                    }
                }
                //*/
                fake_tot += 1;
            }
        }
    }  

    printf("Num locked: %llu\n",num_locked);
    printf("Num not locked: %llu\n",num_not_locked);
    printf("Fake tot: %llu\n",fake_tot);

    //exit(0);

    int index = 0;
    unsigned long long total = 0;
    fprintf(outfile,"%d %d %d\n",nbins,nbins,nbins);
    for(int i = 0; i < nbins; i++)
    {
        printf("%d --------------\n",i);
        fprintf(outfile,"%d\n",i);
        for(int j = 0; j < nbins; j++)
        {
            for(int k = 0; k < nbins; k++)
            {

                index = (nbins)*(nbins)*k + (nbins)*j + i; 
                printf("%lu ",hist[index]);
                fprintf(outfile,"%lu ",hist[index]);
                total += hist[index];
            }
            fprintf(outfile,"\n");
            printf("\n");
        }
    }

    printf("Total: %lu\n",total);

    fclose(outfile);
    return 0;
}  
//////////////////////////////////////////////////////////////////////
