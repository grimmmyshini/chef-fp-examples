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
#ifndef _H_KMEANS
#define _H_KMEANS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#ifndef FLT_MAX
#define FLT_MAX 3.40282347e+38
#endif

#define RANDOM_MAX 2147483647


// /* kmeans_clustering.c */
template <typename T, typename P>
__inline
T euclid_dist_2(P *pt1, P *pt2, int numdims)
{
    int i;
    T ans=0.0;

    for (i=0; i<numdims; i++)
        ans += (pt1[i]-pt2[i]) * (pt1[i]-pt2[i]);

    return(ans);
}

// int find_nearest_point(float  *pt,          /* [nfeatures] */
//                        int     nfeatures,
//                        float **pts,         /* [npts][nfeatures] */
//                        int     npts)
// {
//     int index, i;
//     float min_dist=FLT_MAX;

//     /* find the cluster center id with min distance to pt */
//     for (i=0; i<npts; i++) {
//         float dist;
//         dist = euclid_dist_2(pt, pts[i], nfeatures);  /* no need square root */
//         if (dist < min_dist) {
//             min_dist = dist;
//             index    = i;
//         }
//     }
//     return(index);
// }

// /*----< euclid_dist_2() >----------------------------------------------------*/
// /* multi-dimensional spatial Euclid distance square */



// /*----< kmeans_clustering() >---------------------------------------------*/
// float** kmeans_clustering(float **feature,    /* in: [npoints][nfeatures] */
//                           int     nfeatures,
//                           int     npoints,
//                           int     nclusters,
//                           float   threshold,
//                           int    *membership) /* out: [npoints] */
// {

//     int      i, j, n=0, index, loop=0;
//     int     *new_centers_len; /* [nclusters]: no. of points in each cluster */
//     float    delta;
//     float  **clusters;   /* out: [nclusters][nfeatures] */
//     float  **new_centers;     /* [nclusters][nfeatures] */
  

//     /* allocate space for returning variable clusters[] */
//     clusters    = (float**) malloc(nclusters *             sizeof(float*));
//     clusters[0] = (float*)  malloc(nclusters * nfeatures * sizeof(float));
//     for (i=1; i<nclusters; i++)
//         clusters[i] = clusters[i-1] + nfeatures;

//     /* randomly pick cluster centers */
//     for (i=0; i<nclusters; i++) {
//         //n = (int)rand() % npoints;
//         for (j=0; j<nfeatures; j++)
//             clusters[i][j] = feature[n][j];
// 		n++;
//     }

//     for (i=0; i<npoints; i++)
// 		membership[i] = -1;

//     /* need to initialize new_centers_len and new_centers[0] to all 0 */
//     new_centers_len = (int*) calloc(nclusters, sizeof(int));

//     new_centers    = (float**) malloc(nclusters *            sizeof(float*));
//     new_centers[0] = (float*)  calloc(nclusters * nfeatures, sizeof(float));
//     for (i=1; i<nclusters; i++)
//         new_centers[i] = new_centers[i-1] + nfeatures;
 
  
//     do {
		
//         delta = 0.0;

//         for (i=0; i<npoints; i++) {
// 	        /* find the index of nestest cluster centers */
// 	        index = find_nearest_point(feature[i], nfeatures, clusters, nclusters);
// 	        /* if membership changes, increase delta by 1 */
// 	        if (membership[i] != index) delta += 1.0;

// 	        /* assign the membership to object i */
// 	        membership[i] = index;

// 	        /* update new cluster centers : sum of objects located within */
// 	        new_centers_len[index]++;
// 	        for (j=0; j<nfeatures; j++)          
// 				new_centers[index][j] += feature[i][j];
//         }
      

// 	/* replace old cluster centers with new_centers */
//         for (i=0; i<nclusters; i++) {
//             for (j=0; j<nfeatures; j++) {
//                 if (new_centers_len[i] > 0)
// 					clusters[i][j] = new_centers[i][j] / new_centers_len[i];
// 				new_centers[i][j] = 0.0;   /* set back to 0 */
// 			}
// 			new_centers_len[i] = 0;   /* set back to 0 */
// 		}
            
//         //delta /= npoints;
//     } while (delta > threshold);

  
//     free(new_centers[0]);
//     free(new_centers);
//     free(new_centers_len);

//     return clusters;
// }

// /* cluster.c */
// /*---< cluster() >-----------------------------------------------------------*/
// int cluster(int      numObjects,      /* number of input objects */
//             int      numAttributes,   /* size of attribute of each object */
//             float  **attributes,      /* [numObjects][numAttributes] */
//             int      num_nclusters,
//             float    threshold,       /* in:   */
//             float ***cluster_centres /* out: [best_nclusters][numAttributes] */
    
//             )
// {
//     int     nclusters;
//     int    *membership;
//     float **tmp_cluster_centres;

//     membership = (int*) malloc(numObjects * sizeof(int));

//     nclusters=num_nclusters;

//     srand(7);
	
//     tmp_cluster_centres = kmeans_clustering(attributes,
//                                             numAttributes,
//                                             numObjects,
//                                             nclusters,
//                                             threshold,
//                                             membership);

//     if (*cluster_centres) {
// 		free((*cluster_centres)[0]);
//         free(*cluster_centres);
// 	}
// 	*cluster_centres = tmp_cluster_centres;

   
//     free(membership);

//     return 0;
// }

#endif
