#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<cmath>
#include <unistd.h>

using namespace std;

#define PI 3.14159

///////////////////////////////////////////////////////////////////////////////
// Histogram information

/*
#define S_NBINS 60
#define S_LO 0.
#define S_HI 120.

#define QS_NBINS 16
#define QS_LO 0.9
#define QS_HI 4.1

#define THETA_NBINS 25
#define THETA_LO 0.
#define THETA_HI 1.
*/

//////////////////////////
#define S_NBINS 50
#define S_LO 2.0
#define S_HI 12.0

#define QS_NBINS 16
#define QS_LO 0.9
#define QS_HI 4.1

#define THETA_NBINS 25
#define THETA_LO 0.
#define THETA_HI 1.

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
int distance(float x0, float y0, float z0, float x1, float y1, float z1,float x2, float y2, float z2, int flag, int *totbins)
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
    
    //if (dist0>max_dist) { max_dist = dist0; printf("max_dist: %f\n",max_dist); }
    //if (dist1>max_dist) { max_dist = dist1; printf("max_dist: %f\n",max_dist); }
    //if (dist2>max_dist) { max_dist = dist2; printf("max_dist: %f\n",max_dist); }

    //if (dist0<min_dist) { min_dist = dist0; printf("min_dist: %f\n",min_dist); }
    //if (dist1<min_dist) { min_dist = dist1; printf("min_dist: %f\n",min_dist); }
    //if (dist2<min_dist) { min_dist = dist2; printf("min_dist: %f\n",min_dist); }

    /*
    printf("---------\n");
    printf("%f %f %f\n",x0,x1,x2);
    printf("%f %f %f\n",y0,y1,y2);
    printf("%f %f %f\n",z0,z1,z2);
    printf("%f %f %f\n",dist0,dist1,dist2);
    */

    float s,qs,theta0,theta1,theta2;

    // Sort the distances
    //bool b0 = dist0<dist1;
    //bool b1 = dist1<dist2;
    //bool b2 = dist0<dist2;

    int i0=-1; // shortest
    int i1=-1; // middle
    int i2=-1; // longest

    float shortest,middle,longest;


    // From Stackoverflow 
    // http://stackoverflow.com/questions/13040240/the-fastest-way-to-sort-3-values-java
    ///*
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
    //*/

    /*
       if(dist0<dist1)
       {
       shortest=dist0;
       middle=dist1;
       longest=dist2; // Not really longest, because we never tested for it. 
       }
       if(dist0>=dist1)
       {
       shortest=dist1;
       middle=dist0;
       longest=dist2; // Not really longest, because we never tested for it. 
       }
     */

    int totbin = -999;

    float shortest2 = shortest*shortest;
    float middle2 = middle*middle;
    float longest2 = longest*longest;
    
    for (int k=0;k<3;k++)
    {
        if (k==0) {
            s = shortest;
            qs = middle/shortest;
            theta0 = (acosf((shortest2 + middle2 - longest2)/(2*shortest*middle)))/PI;
            i2 = distance_to_bin(theta0,THETA_LO,THETA_HI,THETA_NBINS,flag);
        } else if (k==1){
            s = middle;
            qs = longest/middle;
            theta1 = (acosf((middle2 + longest2 - shortest2)/(2*middle*longest)))/PI;
            i2 = distance_to_bin(theta1,THETA_LO,THETA_HI,THETA_NBINS,flag);
        } else if (k==2){
            s = shortest;
            qs = longest/shortest;
            //theta2 = (acosf((shortest2 + longest2 - middle2)/(2*shortest*longest)))/PI;
            theta2 = 1.0 - theta0 - theta1;
            i2 = distance_to_bin(theta2,THETA_LO,THETA_HI,THETA_NBINS,flag);
        }

        //printf("%f %f %f\n",s,qs,theta);

        i0 = distance_to_bin(s,S_LO,S_HI,S_NBINS,flag); //Mpc/h, delta s=0.2
        i1 = distance_to_bin(qs,QS_LO,QS_HI,QS_NBINS,flag); // delta qs = 0.2
        //i2 = distance_to_bin(theta,0,1.0,25,flag);

        //printf("%d %d %d\n",i0,i1,i2);
        if (i0<0 || i1<0 || i2<0)
        {
            totbin = -999;
        } else {
            // Combine for big 1d rep of 3d histogram;
            //int nhistbins = nbins;
            //int nhistbins2 = nhistbins*nhistbins;
            totbin = QS_NBINS*THETA_NBINS*i0 + THETA_NBINS*i1 + i2;
        }

        totbins[k] = totbin;
        if (dist0==0 || dist1==0 || dist2==0)
        {
            //printf("dists: %f %f %f\n",dist0,dist1,dist2);
            totbins[k] = -999;
        }

    }

    //return totbin;
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
    sprintf(defaultoutfilename,"default_out.dat");

    char histout[1024];
    sprintf(histout,"0");

    int log_binning_flag = 0; // False

    int flag = 0;
    int voxel_division = -999;
    int voxel_index[3] = {0,0,0};
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

    int size_hist = (S_NBINS)*(QS_NBINS)*(THETA_NBINS);
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
    int bins[3] = {0,0,0};

    int min_index[3] = {0,0,0};
    //int max_index[3] = {NUM_GALAXIES[0],NUM_GALAXIES[1],NUM_GALAXIES[2]};
    int max_index[3] = {ngals[0],ngals[1],ngals[2]};

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
        else if (which_three_input_files==2) // DRD or RDR
            jmin = 0;
        else if (which_three_input_files==3) // DDR or RRD
            jmin = i;
        for(int j=jmin;j<max_index[1];j++)
            //for(int j = 0; j < NUM_GALAXIES[1]; j++)
        {
            int kmin = min_index[2];
            if (which_three_input_files==0)
                //kmin = j+1;
                kmin = j;
            else if (which_three_input_files==1)
                //kmin = j+1;
                kmin = j;
            else if (which_three_input_files==2)
                //kmin = i+1;
                kmin = i;
            else if (which_three_input_files==3)
                kmin = 0;
            for(int k=kmin;k<max_index[2];k++)
                //for(int k =0; k < NUM_GALAXIES[2]; k++)
            {

                //printf("%d ",k);
                    bin_index = distance(h_x[0][i],h_y[0][i],h_z[0][i], \
                            h_x[1][j],h_y[1][j],h_z[1][j], \
                            h_x[2][k],h_y[2][k],h_z[2][k], \
                            flag, bins);

                    //bin_index = distance(htemp_x[0][i],htemp_y[0][i],htemp_z[0][i], \
                            //htemp_x[1][j],htemp_y[1][j],htemp_z[1][j], \
                            //htemp_x[2][k],htemp_y[2][k],htemp_z[2][k], \
                            //hist_min, hist_max, nbins, bin_width, flag, bins);

                    //printf("bin_index: %d\n",bin_index);
                    for (int b=0;b<3;b++)
                    {
                        bin_index = bins[b];
                        if (bin_index>=0)
                        {
                            hist[bin_index]++;
                            temp_hist[bin_index]++;
                            calc_count += 1;
                        }
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
    fprintf(outfile,"%d %d %d\n",S_NBINS,QS_NBINS,THETA_NBINS);
    fprintf(outfile,"%-6.3f %-6.3f %-6.3f\n",S_LO,S_HI,(S_HI-S_LO)/S_NBINS);
    fprintf(outfile,"%-6.3f %-6.3f %-6.3f\n",QS_LO,QS_HI,(QS_HI-QS_LO)/QS_NBINS);
    fprintf(outfile,"%-6.3f %-6.3f %-6.3f\n",THETA_LO,THETA_HI,(THETA_HI-THETA_LO)/THETA_NBINS);
    for(int i = 0; i < S_NBINS; i++)
    {
        printf("%d --------------\n",i);
        //fprintf(outfile,"%d\n",i);
        for(int j = 0; j < QS_NBINS; j++)
        {
            for(int k = 0; k < THETA_NBINS; k++)
            {

                index = (QS_NBINS)*(THETA_NBINS)*i + (THETA_NBINS)*j + k; 
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
