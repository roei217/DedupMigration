/**
Created by Roei Kisous
*/
#pragma once

#include "Node.hpp"
#include "Utility.hpp"
#include "AlgorithmDSManager.hpp"

class HierarchicalClustering final {

private:
    struct ClustersMergeOffer;
    struct ClusteringParams;

public:
    /**
     * creates a new instance of Hierarchical Clustering
     * @param ds - a AlgorithmDSManager's object which contains all matrices and data structures for the algorithm
     * @param lb_sizes - sizes to load balance to, in case you're not using load balance, this argument can be anything
     */
    explicit HierarchicalClustering(std::unique_ptr<AlgorithmDSManager>&  ds, const std::vector<double>& lb_sizes);

    HierarchicalClustering(const HierarchicalClustering&) = delete;
    HierarchicalClustering& operator=(const HierarchicalClustering&) = delete;
    ~HierarchicalClustering() = default;

    /**
     * runs the current instance of Hierarchical Clustering
     * @param load_balance - whether to use load balance
     * @param eps - increment the system size by this every iteration, relevant only with load balance constrain
     * @param number_of_clusters - number of cluster to create (number of output volumes)
     * @param traffics - list of traffics to examine
     * @param seeds - list of seeds to examine
     * @param gaps - list of gaps to examine, in % (0.5, 1, 3 for instance)
     * @param output_path_prefix - output file path prefix
     * Note that in case you want to add volumes to the system, for example from 5 to 7, use your 5 volumes as input
     * and add to this input 2 empty volumes. A total of 7 volumes where 2 are empty
     */
    void run(const bool load_balance,const int eps, const int number_of_clusters, const std::vector<double> &traffics,
             const std::vector<int> &seeds, const std::vector<double> &gaps,
             const std::string& output_path_prefix);

private:
    /**
     * perform the clustering process
     * @param clustering_params - clustering parameters
     * @return - whether the process was successful or not
     */
    bool performClustering(const ClusteringParams& clustering_params);

    /**
     * reset our data structures for the next execution of the algorithm
     */
    void resetDataStructures();

    /**
     * initialize the clusters vector
     */
    void initClusters();

    /**
     * @param initial_clusters - initial clusters as vector of sets
     * @return the calculated clustering result of our execution in a format of
     * mapping between workload name to its final files
     *
     * Should be called only after a successful performClustering
    */
    std::map<std::string, std::set<int>> getClusteringResult(const std::vector<std::set<int>>& initial_clusters) const;

    /**
     * @param final_clusters - the final clustering
     * @param workload_to_cluster_map - the mapping between the initial state clustering to the final clustering
     * @return arranged clustering result by convert the given parameters to a map where every key is a workload name
     * and the value is the set of its files in the clustering result
     */
    std::map<std::string, std::set<int>> getArrangedResult(const std::vector<std::set<int>>& final_clusters,
                                                           const std::vector<int>& workload_to_cluster_map) const;

    /**
     * activates only in case load balancing is in use
     * @param clustering_params - clustering parameters
     * @param cluster1 - a cluster 1 index
     * @param cluster2 - a cluster 2 index
     * @return whether a merge between the clusters 1+2 is valid, size wise
     */
    bool isValidClustersMerge(const ClusteringParams& clustering_params, const int cluster1, const int cluster2) const;

    /**
     * @param cluster1 - a cluster 1 index
     * @param cluster2 - a cluster 2 index
     * @return the size of the cluster in case clusters 1 +2 are merged
     */
    double getMergedClusterSize(const int cluster1, const int cluster2) const;

    /**
     * @param files_indices - set of cluster's files indexes
     * @return the size of the given cluster
     */
    double getClusterSize(const std::set<int>& files_indices) const;

    /**
     * finds a random merge offer within the merge offers with the smallest dissimilarity value
     * (in the lower triangular of dissimilarity_matrix)
     *
     * @param clustering_params - clustering parameters
     * @return a random merge offer from the smallest merge offers
     */
    ClustersMergeOffer findSmallestDissimilarity(const ClusteringParams& clustering_params);

    /**
     * performs complete linkage hierarchical clustering
     * @param merge_offer - the chosen merge offer
     */
    void completeLinkage(const ClustersMergeOffer& merge_offer);

    /**
     * merges the clusters in the given merge_offer (random merge offer within the merge offers with the smallest
     * dissimilarity value)
     *
     * @param merge_offer - the chosen merge offer
     */
    void mergeClusters(const ClustersMergeOffer& merge_offer);

    /**
     * inner function for considering a new merge offer between cluster 1 and cluster 2 when we already have merge
     * offers in the given sorted_merge_offers
     * @param cluster1 - a cluster index
     * @param cluster2 - a cluster index
     * @param sorted_merge_offers - vector of the 'best' merge offers
     * @param clustering_params - clustering parameters
     */
    void updateSortedMergeOffers(const ClusteringParams& clustering_params,
                                 std::vector<ClustersMergeOffer>&  sorted_merge_offers,
                                 const int cluster1,
                                 const int cluster2);

    /**
     * @return the current clustering result where every set in the given vector is a cluster
     */
    std::vector<std::set<int>> getCurrentClustering() const;

    /**
     * ClustersMergeOffer's ascending sort function
     * @param clusters - a clustering where every cluster's set contains its files indices
     * @return returns a clustering where every cluster's set contains all the blocks of its files
     */
    std::vector<std::set<int>> getBlocksInClusters(const std::vector<std::set<int>>& clusters) const;

    /**
     * returns the weighted dissimilarity value between two clusters
     * @param cluster1 - a cluster index
     * @param cluster2 - a cluster index
     * @param clustering_params - clustering parameters
     */
    double getWeightedValueBetweenClusters(const ClusteringParams& clustering_params, const int cluster1,
                                           const int cluster2);

    /**
     * @param initial_clusters - the initial clustering
     * @param final_clusters - the final clustering
     * @return a greedy mapping between the initial clustering to the final clustering based on blocks intersection
     * between the clustering
     */
    std::vector<int> getGreedyWorkloadToClusterMapping(const std::vector<std::set<int>>& initial_clusters,
                                                       const std::vector<std::set<int>>& final_clusters) const;

    /**
     * @param clustering - an arranged clustering mapping as returned from getArrangedResults
     * @return convert the algo based index of the file to the actual file index (SN) of the workload
     */
    std::map<std::string, std::set<int>> getConvertedResults(const std::map<std::string, std::set<int>>& clustering) const;

private:
    /**
     * ClustersMergeOffer's descending sort function
     * @param a - first merge offer
     * @param b - second merge offer
     * @return whether a > b by (weighted_dissimilarity, cluster1, cluster2)
     */
    static bool sortDesc(const ClustersMergeOffer &a, const ClustersMergeOffer &b);

    /**
     * ClustersMergeOffer's ascending sort function
     * @param a - first merge offer
     * @param b - second merge offer
     * @return whether a < b by (weighted_dissimilarity, cluster1, cluster2)
     */
    static bool sortAsc(const ClustersMergeOffer &a, const ClustersMergeOffer &b);

    /**
     * ClustersMergeOffer's ascending sort function
     * @param clusters1_blocks - a clustering where every cluster's ser contains all the blocks of its files
     * @param clusters2_blocks - a clustering where every cluster's set contains all the blocks of its files
     * @return a intersection matrix where every cell [i,j] holds the number of shared blocks between
     * cluster i and cluster j
     */
    static std::vector<std::vector<int>> getBlocksIntersectionOfClusters(
            const std::vector<std::set<int>>& clusters1_blocks,
            const std::vector<std::set<int>>& clusters2_blocks);

    /**
     * finds the best suited cluster and original volume to match in our assignment of cluster to volumes
     * @param mat - the matrix of blocks' number intersection between clusters and original volumes
     * @return the best suited cluster and original volume to match
     */
    static std::pair<int, int> findMax(std::vector<std::vector<int>> mat);

    /**
     * @param sorted_merge_offers - list of sorted merge offers
     * @param gap - gap param of the clustering
     * @return a randomized offer from all offers in sorted_merge_offers that within the gap range
     */
    static ClustersMergeOffer randClustersMergeOfferFromVector(
            const std::vector<ClustersMergeOffer>& sorted_merge_offers,
            const double gap);

    /**
     * @param clustering_params - clustering params
     * @param result_mapping - the final result mapping of the clustering
     * @param initial_mapping - the initial mapping of the system
     * @param output_path_prefix - the prefix for the output file path
     * @return write clustering result to csv which its name derived from the clustering params
     */
    static void outputConvertedResultsToCsv(const ClusteringParams& clustering_params,
                                            const std::map<std::string,std::set<int>>& result_mapping,
                                            const std::map<std::string,std::set<int>>& initial_mapping,
                                            const std::string& output_path_prefix);

    static std::string getResultFileName(const double traffic, const int seed, const double gap,
                                    const bool load_balance,const int eps, const std::string& output_path_prefix);

// m_clusters - list of nodes which represent the cluster. each cluster contains set of all the files in it
// m_lb_sizes - sizes to load balance to, relevant for load balance only
// m_ds - a AlgorithmDSManager's object which contains all the data structures for the algorithm
private:
    std::vector<std::unique_ptr<Node>> m_clusters;
    const std::vector<double> m_lb_sizes;
    const std::unique_ptr<AlgorithmDSManager> m_ds;
};
