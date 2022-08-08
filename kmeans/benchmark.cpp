#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <sys/types.h>
#include <fcntl.h>

#include "benchmark/benchmark.h"

#include "kmeans.h"
#include "kmeans-adapt.h"

#include "clad/Differentiator/Differentiator.h"
#include "../PrintModel/ErrorFunc.h"

#include "Derivative.h"

#include "adapt.h"

#define INPUT_FILE_NAME "data/3000.txt"
#define IS_BINARY_FILE 0

struct cout_redirect
{
    cout_redirect(std::streambuf *new_buffer)
        : old(std::cout.rdbuf(new_buffer))
    {
    }

    ~cout_redirect()
    {
        std::cout.rdbuf(old);
    }

private:
    std::streambuf *old;
};

static void ErrorEstimateKMeansAdapt(benchmark::State &state)
{
    int opt;
    extern char *optarg;
    extern int optind;
    int nclusters = 5;
    const char *filename = INPUT_FILE_NAME;
    float *buf;
    AD_real **attributes;
    AD_real **cluster_centres = NULL;

    int numAttributes;
    int numObjects;
    char line[1024];
    int isBinaryFile = IS_BINARY_FILE;
    int nloops;
    float threshold = 0.001;
    double timing;

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

    std::stringstream suppress;
    cout_redirect output(suppress.rdbuf());

    for (auto _ : state)
    {
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
    }

    //------------------ Deallocate memory after this line -------------//
    free(new_centers[0]);
    free(new_centers);
    free(new_centers_len);
    free(membership);
    free(attributes);
    free(cluster_centres[0]);
    free(cluster_centres);
    free(buf);
}

static void ErrorEstimateKMeansClad(benchmark::State &state)
{
    // auto df = clad::estimate_error(euclid_dist_2<double, double>);
int opt;
    extern char *optarg;
    extern int optind;
    int nclusters = 5;
    const char *filename = INPUT_FILE_NAME;
    float *buf;
    double **attributes;
    double **cluster_centres = NULL;

    int numAttributes;
    int numObjects;
    char line[1024];
    int isBinaryFile = IS_BINARY_FILE;
    int nloops;
    float threshold = 0.001;
    double timing;

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

    attributes = (double **)malloc(numObjects * sizeof(double *));
    attributes[0] = (double *)malloc(numObjects * numAttributes * sizeof(double));

    for (int i = 1; i < numObjects; i++)
        attributes[i] = attributes[i - 1] + numAttributes;

    // memcpy(attributes[0], buf, numObjects * numAttributes * sizeof(float));
    for (int i = 0; i < numObjects * numAttributes; i++)
        attributes[0][i] = buf[i];

    int *membership = (int *)malloc(numObjects * sizeof(int));

    double **clusters; /* out: [nclusters][numAttributes] */
    clusters = (double **)malloc(nclusters * sizeof(double *));
    clusters[0] = (double *)malloc(nclusters * numAttributes * sizeof(double));

    double **new_centers; /* [nclusters][numAttributes] */
    new_centers = (double **)malloc(nclusters * sizeof(double *));
    new_centers[0] = (double *)malloc(nclusters * numAttributes * sizeof(double));

    int *new_centers_len; /* [nclusters]: no. of points in each cluster */
    new_centers_len = (int *)malloc(nclusters * sizeof(int));

    cluster_centres = NULL;

    //-------- Allocate memory before this line -------------//
    std::stringstream suppress;
    cout_redirect output(suppress.rdbuf());

for (auto _ : state)
    {
    AD_begin();
    clad::resetErrors();

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
                double min_dist = FLT_MAX;

                /* find the cluster center id with min distance to pt */
                for (int k = 0; k < nclusters; k++)
                {
                    double dist = 0.0;

                    dist = euclid_dist_2<double, double>(attributes[i], clusters[k], numAttributes); /* no need square root */

                    // Error Estimation Begin
                    clad::array<double> d_clusters(numAttributes);
                    clad::array<double> d_attributes(numAttributes);
                    int d_numAttributes = 0;
                    double final_error = 0;

                    euclid_dist_2_grad(
                        attributes[i], clusters[k], numAttributes,
                        d_attributes, d_clusters, &d_numAttributes, final_error);

                    // printf("Final error: %f of object %d and cluster %d\n", final_error, i, k);
                    // printf("Actual error: %f\n", std::fabs(euclid_dist_2<float, double>(attributes[i], clusters[k], numAttributes) - dist));
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
    benchmark::DoNotOptimize(cluster_centres);

    clad::printErrorReport();
}

    //------------------ Deallocate memory after this line -------------//
    free(new_centers[0]);
    free(new_centers);
    free(new_centers_len);
    free(membership);
    free(attributes);
    free(cluster_centres[0]);
    free(cluster_centres);
    free(buf);
}

BENCHMARK(ErrorEstimateKMeansClad)->Unit(benchmark::kSecond)->Iterations(1);
BENCHMARK(ErrorEstimateKMeansAdapt)->Unit(benchmark::kSecond)->Iterations(1);

BENCHMARK_MAIN();