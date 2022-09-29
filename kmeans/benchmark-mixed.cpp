#include "benchmark/benchmark.h"
#include "../benchmark-utils/memory-manager.hpp"

#include <iostream>
#include <sstream>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <fcntl.h>

#include "kmeans.h"


#define IS_BINARY_FILE 0
#define NLOOPS 1

const char * input_files[] = {"data/100.txt", "data/1000.txt", "data/10000.txt", "data/100000.txt", "data/1000000.txt"};

struct cout_suppressor
{
  cout_suppressor()
      : buffer(), old(std::cout.rdbuf(buffer.rdbuf()))
  {
  }

  ~cout_suppressor()
  {
    std::cout.rdbuf(old);
  }

private:
  std::stringstream buffer;
  std::streambuf *old;
};

void EuclidDist_MixedPrecision(benchmark::State &state)
{
    const int numdims = state.range(1);
    const int numpoints = state.range(0);
    srand(time(NULL));
    float **pt1 = new float*[numpoints];
    double **pt2 = new double*[numpoints];

    for (int i = 0; i < numpoints; i++) {
        pt1[i] = new float[numdims];
        pt2[i] = new double[numdims];
        for (int j = 0; j < numdims; j++) {
            pt1[i][j] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            pt2[i][j] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        }
    }

    for (auto _ : state)
    {
        for (int i = 0; i < numpoints; i++) {
            double dist = euclid_dist_2<double, float>(pt1[i], pt2[i], numdims);
            benchmark::DoNotOptimize(dist);
        }
    }
    
    for (int i = 0; i < numpoints; i++) {
        delete[] pt1[i];
        delete[] pt2[i];
    }

    delete[] pt1;
    delete[] pt2;
}

void EuclidDist_HighPrecision(benchmark::State &state)
{
    const int numdims = state.range(1);
    const int numpoints = state.range(0);
    srand(time(NULL));
    double **pt1 = new double*[numpoints];
    double **pt2 = new double*[numpoints];

    for (int i = 0; i < numpoints; i++) {
        pt1[i] = new double[numdims];
        pt2[i] = new double[numdims];
        for (int j = 0; j < numdims; j++) {
            pt1[i][j] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            pt2[i][j] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        }
    }

    for (auto _ : state)
    {
        for (int i = 0; i < numpoints; i++) {
            double dist = euclid_dist_2<double, double>(pt1[i], pt2[i], numdims);
            benchmark::DoNotOptimize(dist);
        }
    }
    
    for (int i = 0; i < numpoints; i++) {
        delete[] pt1[i];
        delete[] pt2[i];
    }

    delete[] pt1;
    delete[] pt2;
}

void KMeans_MixedPrecision(benchmark::State &state)
{
    int opt;
    extern char *optarg;
    extern int optind;
    int nclusters = 5;
    const char *filename = input_files[state.range(0)];
    float *buf;
    float **attributes;
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

    nloops = NLOOPS;

    attributes = (float **)malloc(numObjects * sizeof(double *));
    attributes[0] = (float *)malloc(numObjects * numAttributes * sizeof(double));

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
    cout_suppressor suppressor;

    for (auto _ : state)
    {
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

                        dist = euclid_dist_2<double, float>(attributes[i], clusters[k], numAttributes); /* no need square root */

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

void KMeans_HighPrecision(benchmark::State &state)
{
    int opt;
    extern char *optarg;
    extern int optind;
    int nclusters = 5;
    const char *filename = input_files[state.range(0)];
    double *buf;
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
        buf = (double *)malloc(numObjects * numAttributes * sizeof(double));

        fread(buf, sizeof(double), numObjects * numAttributes, infile);

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
        buf = (double *)malloc(numObjects * numAttributes * sizeof(double));

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

    nloops = NLOOPS;

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
    cout_suppressor suppressor;

    for (auto _ : state)
    {
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

BENCHMARK(EuclidDist_MixedPrecision)->Args({10, 3})->Args({10, 10})->Args({100, 10})->Args({100, 100})->Args({1000, 100});
BENCHMARK(EuclidDist_HighPrecision)->Args({10, 3})->Args({10, 10})->Args({100, 10})->Args({100, 100})->Args({1000, 100});
BENCHMARK(KMeans_MixedPrecision)->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4);
BENCHMARK(KMeans_HighPrecision)->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4);

// BENCHMARK_MAIN();
int main(int argc, char** argv)
{
    ::benchmark::RegisterMemoryManager(mm.get());
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::RegisterMemoryManager(nullptr);
}
