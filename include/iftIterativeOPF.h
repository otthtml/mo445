#ifndef IFT_ITERATIVE_OPF_H_
#define IFT_ITERATIVE_OPF_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ift.h"
#include "iftDataSet.h"

typedef struct ift_graphnode {
    /**  Maximum arc weight from the node to its neighbors */
    float     maxarcw;
    /** Corresponding root node in the graph */
    int       root;
    /** Corresponding training sample in the original dataset */
    int       sample;
    /** List of adjacent nodes */
    iftAdjSet *adj;
    /** List of adjacent nodes on plateaus of density */
    iftSet    *adjplat;
    /** Predecessor node */
    int       pred;     // predecessor node
} iftGraphNode;

typedef struct ift_graph {
    /** Is the graph complete?  */
    bool         is_complete;
    /** List of nodes in the graph */
    iftGraphNode *node;
    /** List of path value of the nodes */
    float        *pathval;
    /** List of nodes ordered by its path value */
    int          *ordered_nodes;
    /** Number of nodes of the graph */
    int          nnodes;
    /** Priority queue */
    iftFHeap     *Q;
    /** Corresponding dataset */
    iftDataSet   *Z;
    /** Centroids */
    int          *centroids;
} iftGraph;

typedef enum ift_connectivity_function_type { 
    IFT_MAX, IFT_SUM } iftConnFuncType;

typedef enum ift_centers_update_type { 
    IFT_DATA_FEAT, IFT_IMG_COORD } iftCentUpdateType;

/**
 * @brief Updates the centroids using the feature vectors
 * @param graph --- input graph
 * @param S --- input dynamic sets
 * @param centroids --- centroids to be updated
 * @param n_clusters --- number of dynamic sets (clusters)
 * @author David Aparco Cardenas 
 * @date Jan 10th, 2020
 */

void iftUpdateCentroids(iftGraph *graph, iftDynamicSet **S, int *centroids, int n_clusters);

/**
 * @brief Updates the centroids using the image coordinates X, Y
 * @param graph --- input graph
 * @param S --- input dynamic sets
 * @param centroids --- centroids to be updated
 * @param n_clusters --- number of dynamic sets (clusters)
 * @author David Aparco Cardenas 
 * @date Jan 10th, 2020
 */

void iftMUpdateCentroids(iftGraph *graph, iftDynamicSet **S, int *centroids, int n_clusters);

/**
 * @brief Updates a single centroid using the nearest node to the mean
 * @param graph --- input graph
 * @param S --- dynamic set
 * @return Updated centroid
 * @author David Aparco Cardenas 
 * @date Jan 10th, 2020
 */

int iftUpdateSingleCentroid(iftGraph *graph, iftDynamicSet *S);

/**
 * @brief Inserts a new element in a dynamic set
 * @param graph --- input graph
 * @param S --- input dynamic set
 * @param p --- new element to be inserted
 * @author David Aparco Cardenas 
 * @date Jan 10th, 2020
 */

void iftInsertSetDynamicSet(iftGraph *graph, iftDynamicSet *S, int p);

/**
 * @brief Clusterizes a dataset using the Iterative OPF algorithm.
 * @param graph --- input graph
 * @param n_clusters --- number of clusters
 * @param max_iter --- maximum number of iterations
 * @param conn_function_type --- connectivity function type: IFT_MAX, IFT_SUM 
 * @param cent_update_type --- center update policy: centroids (IFT_DATA_FEAT), x-y coordinates (IFT_IMG_COORD)
 * @author David Aparco Cardenas 
 * @date Jan 10th, 2020
 */
  
void iftIterativeOPF(iftGraph *graph, int n_clusters, int max_iter, iftConnFuncType conn_function_type, iftCentUpdateType cent_update_type);

/**
 * @brief Clusterizes a graph using the Iterated Watersheds algorithm.
 * @param graph --- input graph
 * @param n_clusters --- desired number of clusters
 * @param max_iter --- maximum number of iterations
 * @param conn_function_type --- connectivity function type: IFT_MAX, IFT_SUM
 * @param cent_update_type --- center update policy: centroids (IFT_DATA_FEAT), x-y coordinates (IFT_IMG_COORD)
 * @author David Aparco Cardenas 
 * @date Jan 10th, 2020
 */
  
void iftIteratedWatersheds(iftGraph *graph, int n_clusters, int max_iter, iftConnFuncType conn_function_type, iftCentUpdateType cent_update_type);

/**
 * @brief Creates a graph.
 * @param Z --- input dataset
 * @return Created graph iftGraph
 * @author David Aparco Cardenas 
 * @date Jan 30th, 2020
 */

iftGraph *iftCreateGraph(iftDataSet *Z);

/**
 * @brief Destroys a graph.
 * @param graph --- input graph
 * @author David Aparco Cardenas 
 * @date Jan 30th, 2019
 */

void iftDestroyGraph(iftGraph **graph);

/**
 * @brief Sets the graph topology from a file.
 * @param pathname --- delaunay triangulation file
 * @param graph --- input graph
 * @author David Aparco Cardenas 
 * @date Jan 30th, 2020
 */

void iftReadDelaunayTriangulation(const char *pathname, iftGraph *graph);

/**
 * @brief Sets the adjacency sets for each node in the graph from MImage.
 * @param graph --- input graph
 * @param mimg --- input MImage
 * @param A --- input adjacency relation 
 * @author David Aparco Cardenas 
 * @date Feb 18th, 2020
 */

void iftSetMGraphAdjacencySets(iftGraph *graph, iftMImage *mimg, iftAdjRel *A);

/**
 * @brief Update the optimum path cost based on the global path cost.
 * @param graph --- input graph
 * @param opt_labels --- optimum group labels
 * @param opt_centroids --- optimum centroid
 * @param opt_cost --- optimum global path cost
 * @author David Aparco Cardenas 
 * @date Feb 26th, 2020
 */

void iftUpdateOptimumPathCost(iftGraph *graph, iftIntArray *opt_labels, int *opt_centroids, float *opt_cost);

/**
 * @brief Update the optimum path cost based on the global path cost.
 * @param graph --- input graph
 * @param opt_labels --- optimum group labels
 * @param mean_centroids --- centroids for k-means
 * @param opt_centroids --- optimum centroids
 * @param opt_cost --- optimum group labels
 * @author David Aparco Cardenas 
 * @date Feb 26th, 2020
 */

void iftUpdateOptimumPathCostKMeans(iftGraph *graph, iftIntArray *opt_labels, float **mean_centroids, float **opt_centroids, float *opt_cost);

/**
 * @brief Iterative optimum-path forest algorithm.
 * @param graph --- input graph
 * @param conn_function_type --- connectivity function type: IFT_FMAX, IFT_FSUM 
 * @author David Aparco Cardenas 
 * @date Feb 26th, 2020
 */

void iftDynamicOptimumPathForest(iftGraph *graph, iftDynamicSet **S, iftConnFuncType conn_function_type);

/**
 * @brief Optimum-path forest algorithm.
 * @param graph --- input graph
 * @param S --- input dynamic sets
 * @param conn_function_type --- connectivity function type: IFT_FMAX, IFT_FSUM
 * @author David Aparco Cardenas 
 * @date Feb 26th, 2020
 */

void iftOptimumPathForest(iftGraph *graph, iftDynamicSet **S, iftConnFuncType conn_function_type);

/**
 * @brief Creates a connected graph with a knn adjacency relation
 * @param Z --- input dataset
 * @param n_clusters --- k nearest neighbors
 * @return Connected graph iftGraph
 * @author David Aparco Cardenas 
 * @date Mar 20th, 2020
 */

iftGraph *iftCreateConnectedKnnGraph(iftDataSet *Z, int k);

/**
 * @brief Label connected components of a graph
 * @param graph --- input graph
 * @author David Aparco Cardenas 
 * @date Mar 20th, 2020
 */

void iftLabelConnectedComponents(iftGraph *graph);

/**
 * @brief Check the graph symmetry, i.e. if node u is adjacent to node v, then node v is adjacent to node u
 * @param graph --- input graph
 * @author David Aparco Cardenas 
 * @date Mar 20th, 2020
 */

int iftCheckGraphSymmetry(iftGraph *graph);

/**
 * @brief selects initial centroids using the k++ algorithm
 * @param graph --- input graph
 * @param n_clusters --- number of centroids
 * @author David Aparco Cardenas 
 * @date Apr 9th, 2020
 */

void iftSetInitialCenters(iftGraph *graph, int n_clusters);

#ifdef __cplusplus
}
#endif

#endif
