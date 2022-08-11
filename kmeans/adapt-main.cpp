/*****************************************************************************/
/*IMPORTANT:  READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.         */
/*By downloading, copying, installing or using the software you agree        */
/*to this license.  If you do not agree to this license, do not download,    */
/*install, copy or use the software.                                         */
/*                                                                           */
/*                                                                           */
/*Copyright (c) 2005 Northwestern University                                 */
/*All rights reserved.                                                       */

/*Redistribution of the software in source and binary forms,                 */
/*with or without modification, is permitted provided that the               */
/*following conditions are met:                                              */
/*                                                                           */
/*1       Redistributions of source code must retain the above copyright     */
/*        notice, this list of conditions and the following disclaimer.      */
/*                                                                           */
/*2       Redistributions in binary form must reproduce the above copyright   */
/*        notice, this list of conditions and the following disclaimer in the */
/*        documentation and/or other materials provided with the distribution.*/
/*                                                                            */
/*3       Neither the name of Northwestern University nor the names of its    */
/*        contributors may be used to endorse or promote products derived     */
/*        from this software without specific prior written permission.       */
/*                                                                            */
/*THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS    */
/*IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED      */
/*TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT AND         */
/*FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL          */
/*NORTHWESTERN UNIVERSITY OR ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT,       */
/*INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES          */
/*(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR          */
/*SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)          */
/*HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,         */
/*STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN    */
/*ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE             */
/*POSSIBILITY OF SUCH DAMAGE.                                                 */
/******************************************************************************/
/*************************************************************************/
/**   File:         example.c                                           **/
/**   Description:  Takes as input a file:                              **/
/**                 ascii  file: containing 1 data point per line       **/
/**                 binary file: first int is the number of objects     **/
/**                              2nd int is the no. of features of each **/
/**                              object                                 **/
/**                 This example performs a fuzzy c-means clustering    **/
/**                 on the data. Fuzzy clustering is performed using    **/
/**                 min to max clusters and the clustering that gets    **/
/**                 the best score according to a compactness and       **/
/**                 separation criterion are returned.                  **/
/**   Author:  Wei-keng Liao                                            **/
/**            ECE Department Northwestern University                   **/
/**            email: wkliao@ece.northwestern.edu                       **/
/**                                                                     **/
/**   Edited by: Jay Pisharath                                          **/
/**              Northwestern University.                               **/
/**                                                                     **/
/**   ================================================================  **/
/**																		**/
/**   Edited by: Sang-Ha  Lee											**/
/**				 University of Virginia									**/
/**																		**/
/**   Description:	No longer supports fuzzy c-means clustering;	 	**/
/**					only regular k-means clustering.					**/
/**					Simplified for main functionality: regular k-means	**/
/**					clustering.											**/
/**                                                                     **/
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <sys/types.h>
#include <fcntl.h>
#include "getopt.h"

#include "kmeans.h"
#include "kmeans-adapt.h"

#include "adapt.h"


/*---< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int opt;
    extern char *optarg;
    extern int optind;
    int nclusters = 5;
    char *filename = 0;
    float *buf;
    AD_real **attributes;
    AD_real **cluster_centres = NULL;

    int numAttributes;
    int numObjects;
    char line[1024];
    int isBinaryFile = 0;
    int nloops;
    float threshold = 0.001;
    double timing;

    while ((opt = getopt(argc, argv, "i:k:t:b")) != EOF)
    {
        switch (opt)
        {
        case 'i':
            filename = optarg;
            break;
        case 'b':
            isBinaryFile = 1;
            break;
        case 't':
            threshold = atof(optarg);
            break;
        case 'k':
            nclusters = atoi(optarg);
            break;
        case '?':
            usage(argv[0]);
            break;
        default:
            usage(argv[0]);
            break;
        }
    }

    if (filename == 0)
        usage(argv[0]);

    numAttributes = numObjects = 0;

    /* from the input file, get the numAttributes and numObjects ------------*/

    if (isBinaryFile)
    {
        FILE *infile;
        if ((infile = fopen(filename, "r")) == NULL)
        {
            fprintf(stderr, "Error: no such file (%s)\n", filename);
            exit(1);
        }
        fread(&numObjects, sizeof(int), 1, infile);
        fread(&numAttributes, sizeof(int), 1, infile);

        /* allocate space for attributes[] and read attributes of all objects */
        buf = (float *)malloc(numObjects * numAttributes * sizeof(float));

        fread(buf, sizeof(float), numObjects * numAttributes, infile);

        fclose(infile);
    }
    else
    {
        FILE *infile;
        if ((infile = fopen(filename, "r")) == NULL)
        {
            fprintf(stderr, "Error: no such file (%s)\n", filename);
            exit(1);
        }
        while (fgets(line, 1024, infile) != NULL)
            if (strtok(line, " \t\n") != 0)
                numObjects++;
        rewind(infile);
        while (fgets(line, 1024, infile) != NULL)
        {
            if (strtok(line, " \t\n") != 0)
            {
                /* ignore the id (first attribute): numAttributes = 1; */
                while (strtok(NULL, " ,\t\n") != NULL)
                    numAttributes++;
                break;
            }
        }

        /* allocate space for attributes[] and read attributes of all objects */
        buf = (float *)malloc(numObjects * numAttributes * sizeof(float));

        rewind(infile);

        int buf_ind = 0;

        while (fgets(line, 1024, infile) != NULL)
        {
            if (strtok(line, " \t\n") == NULL)
                continue;
            for (int j = 0; j < numAttributes; j++)
            {
                buf[buf_ind] = atof(strtok(NULL, " ,\t\n"));
                buf_ind++;
            }
        }
        fclose(infile);
    }

    nloops = 1;
    printf("I/O completed\n");

    attributes = (AD_real **)malloc(numObjects * sizeof(AD_real *));
    attributes[0] = (AD_real *)malloc(numObjects * numAttributes * sizeof(AD_real));

    for (int i = 1; i < numObjects; i++)
        attributes[i] = attributes[i - 1] + numAttributes;

    // memcpy(attributes[0], buf, numObjects * numAttributes * sizeof(float));
    for (int i = 0; i < numObjects * numAttributes; i++)
        attributes[0][i] = buf[i];

    int *membership = (int *)malloc(numObjects * sizeof(int));

    AD_real **clusters; /* out: [nclusters][numAttributes] */
    clusters = (AD_real **)malloc(nclusters * sizeof(AD_real *));
    clusters[0] = (AD_real *)malloc(nclusters * numAttributes * sizeof(AD_real));

    AD_real **new_centers; /* [nclusters][numAttributes] */
    new_centers = (AD_real **)malloc(nclusters * sizeof(AD_real *));
    new_centers[0] = (AD_real *)malloc(nclusters * numAttributes * sizeof(AD_real));

    int *new_centers_len; /* [nclusters]: no. of points in each cluster */
    new_centers_len = (int *)malloc(nclusters * sizeof(int));

    cluster_centres = NULL;

    //-------- Allocate memory before this line -------------//

    AD_begin();

    for (int loop = 0; loop < nloops; loop++)
    {

        // cluster(numObjects,
        //         numAttributes,
        //         attributes,           /* [numObjects][numAttributes] */
        //         nclusters,
        //         threshold,
        //         &cluster_centres
        //        );

        srand(7);

        // tmp_cluster_centres = kmeans_clustering(attributes,
        //                                         numAttributes,
        //                                         numObjects,
        //                                         nclusters,
        //                                         threshold,
        //                                         membership);
        int n = 0, index;
        float delta;

        /* allocate space for returning variable clusters[] */
        for (int i = 1; i < nclusters; i++)
            clusters[i] = clusters[i - 1] + numAttributes;

        /* randomly pick cluster centers */
        for (int i = 0; i < nclusters; i++)
        {
            // n = (int)rand() % numObjects;
            for (int j = 0; j < numAttributes; j++)
                clusters[i][j] = attributes[n][j];
            n++;
        }

        for (int i = 0; i < numObjects; i++)
            membership[i] = -1;

        /* need to initialize new_centers_len and new_centers[0] to all 0 */
        *new_centers_len = 0;

        for (int i = 0; i < nclusters * numAttributes; i++)
            new_centers[0][i] = 0.0;

        for (int i = 1; i < nclusters; i++)
            new_centers[i] = new_centers[i - 1] + numAttributes;

        do
        {
            delta = 0.0;

            for (int i = 0; i < numObjects; i++)
            {
                /* find the index of nestest cluster centers */

                // index = find_nearest_point(attributes[i], numAttributes, clusters, nclusters);
                AD_real min_dist = FLT_MAX;

                /* find the cluster center id with min distance to pt */
                for (int k = 0; k < nclusters; k++)
                {
                    AD_real dist = 0.0;

                    for (int l = 0; l < numAttributes; l++) {
                        AD_INDEPENDENT(attributes[i][l], "attributes");
                        AD_INDEPENDENT(clusters[k][l], "clusters");
                    }

                    dist = adapt::euclid_dist_2(attributes[i], clusters[k], numAttributes); /* no need square root */

                    AD_DEPENDENT(dist, "dist", 0.0);

                    // Error Estimation End
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        index = k;
                    }
                }

                /* if membership changes, increase delta by 1 */
                if (membership[i] != index)
                    delta += 1.0;

                /* assign the membership to object i */
                membership[i] = index;

                /* update new cluster centers : sum of objects located within */
                new_centers_len[index]++;
                for (int j = 0; j < numAttributes; j++)
                    new_centers[index][j] += attributes[i][j];
            }

            /* replace old cluster centers with new_centers */
            for (int i = 0; i < nclusters; i++)
            {
                for (int j = 0; j < numAttributes; j++)
                {
                    if (new_centers_len[i] > 0)
                        clusters[i][j] = new_centers[i][j] / new_centers_len[i];
                    new_centers[i][j] = 0.0; /* set back to 0 */
                }
                new_centers_len[i] = 0; /* set back to 0 */
            }

            // delta /= numObjects;
        } while (delta > threshold);

        cluster_centres = clusters;

        // kmeans_clustering end
        // cluster end
    }

    AD_report();
    AD_end();

    printf("number of Clusters %d\n", nclusters);
    printf("number of Attributes %d\n\n", numAttributes);
    /*printf("Cluster Centers Output\n");
    printf("The first number is cluster number and the following data is arribute value\n");
    printf("=============================================================================\n\n");

    for (i=0; i<nclusters; i++) {
        printf("%d: ", i);
        for (j=0; j<numAttributes; j++)
            printf("%f ", cluster_centres[i][j]);
        printf("\n\n");
    }*/
    printf("Time for process: %f\n", timing);

    //------------------ Deallocate memory after this line -------------//
    free(new_centers[0]);
    free(new_centers);
    free(new_centers_len);
    free(membership);
    free(attributes);
    free(cluster_centres[0]);
    free(cluster_centres);
    free(buf);
    return (0);
}
