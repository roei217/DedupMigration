/**
Created by Roei Kisous
*/
#pragma once

#include "Utility.hpp"

#include <fstream>
#include <utility>
#include <vector>
#include <queue>
#include <unordered_set>
#include <unordered_map>

class AlgorithmDSManager final {
public:
    /**
     * DissimilarityCell struct - a cell value for each cell in dissimilarity matrix
     */
    struct DissimilarityCell {
        DissimilarityCell() : jaccard_distance(0), origin_clusters()
        {
        }

        DissimilarityCell(const double jaccard_distance, std::unordered_set<int> origin_clusters) :
        jaccard_distance(jaccard_distance),
        origin_clusters(std::move(origin_clusters))
        {
        }

        double jaccard_distance;
        std::unordered_set<int> origin_clusters;
    };

public:
    /**
     * @param workloads_paths - paths' vector to the workloads
     * @param requested_number_of_fingerprints - how many fps to use
     * @param load_balance - whether to use load balance or not
     */
    explicit AlgorithmDSManager(const std::vector<std::string>& workloads_paths,
                                const int requested_number_of_fingerprints,
                                const bool load_balance);

    AlgorithmDSManager(const AlgorithmDSManager&) = delete;
    AlgorithmDSManager& operator=(const AlgorithmDSManager&) = delete;
    ~AlgorithmDSManager() = default;
    /**
     *
     * @return - the initial state files' algorithm sn to cluster mapping
     */
    std::unordered_map<int,int> getInitialFileToClusterMapping();

    /**
     *
     * resets the dissimilarities matrix's bottom triangle by copying corresponding values from its upper triangle
     */
    void resetDMBottomTriangle();

    /**
     *
     * @param cluster1 - a cluster index
     * @param cluster2 - a cluster index
     * @return the dissimilarity cell between cluster1 and cluster2 (as written in the dissimilarities matrix)
     */
    DissimilarityCell& getDissimilarityCell(const int cluster1, const int cluster2);

    /**
     *
     * @param cluster1 - a cluster index
     * @param cluster2 - a cluster index
     * @param dissimilarity_cell - a DissimilarityCell object
     * set the relevant cluster1,cluster2 dissimilarities matrix's cell to be the given dissimilarity_cell
     */
    void setDissimilarityCell(const int cluster1, const int cluster2, const DissimilarityCell& dissimilarity_cell);

    /**
     * deactivate the given cluster in the dissimilarity matrix
     * @param cluster_index - a cluster index
     */
    void deactivateClusterInDissimilarityMat(const int cluster_index);

    /**
     *
     * @return number of fingerprints for clustering
     */
    int getNumberOfFingerprintsForClustering() const;

    /**
     *
     * @return number of files for clustering
     */
    int getNumberOfFilesForClustering() const;

    /**
     * @param fp_index - a fingerprint index
     * @return the size of the corresponding fingerprint as given by fp_index
     */
    int getFingerprintSize(const int fp_index) const;

    /**
     * @param file_index - a file index
     * @param fp_index - a fingerprint index
     * @return whether file(file_index) contains fingerprint (fp_index). use the appearances' matrix
     */
    bool isFileHasFingerprint(const int file_index, const int fp_index) const;

    /**
     * @param file_index - a file index
     * @return the SN of the give file index
     */
    int getFileSN(const int file_index) const;

    /**
     * @return the initial clusters
     */
    std::vector<std::set<int>> getInitialClusters() const;

    /**
     * @param clusters_as_vector - clusters as vector. for example, as getInitialClusters results
     * @return the clusters as map where the key is the workload name and the value is the set of its files
     */
    std::map<std::string, std::set<int>> getClustersAsMap(const std::vector<std::set<int>>& clusters_as_vector) const;

    /**
     * @param workload_index - a workload index
     * @return the full path of the corresponding workload
     */
    std::string getWorkloadFullPath(const int workload_index) const;

    /**
     * @return num of workloads
     */
    int getNumOfWorkloads() const;

    /**
     * @return the initial system size
     */
    int getInitialSystemSize() const;

    /**
     * @return the optimal system size
     */
    int getOptimalSystemSize() const;

private:
    /**
     * initialize the Fingerprint related Data structures and returns the selected fps ordered by their SN
     */
    std::vector<int> initializeAndSelectFingerprints();

    /**
     * initialize the appearances matrix with the initial+optimal system size with deduplication fields
     * @param fingerprints_for_clustering_ordered_by_SN - fingerprints for clustering ordered by their SN
     * @param load_balance - whether to calculate load balance stats
     * @return - pair of initial system size and optimal system size
     */
    void initializeAppearancesMatrix(const std::vector<int>& fingerprints_for_clustering_ordered_by_SN,
                                     bool load_balance);

    /**
     * initialize the final dissimilarity matrix
     */
    void initializeDissimilaritiesMatrix();

    /**
     * clears the appearances matrix
     */
    void clearAppearancesMatrix();

    /**
     * returns vector of workloads' streams (files located at m_workloads_paths)
     */
    std::vector<std::ifstream> getWorkloadsStreams() const;

    /**
     * inner function of initializeAndSelectFingerprints
     * updates the given current_blocks,current_algo_file_index  and minhash_fingerprints
     * with the given workload file
     * @param workload_file - a workload's file
     * @param workload_index - index of the workload's file
     * @param current_algo_file_index - the current algo file index
     * @param current_blocks - set of blocks we already handled
     * @param requested_num_fps - threshold of number of fingerprints we want to collect
     * @param minhash_fingerprints - priority queue of fingerprints
     */
    void updateFPWithWorkload(std::ifstream& workload_file, const int workload_index,
                              int& current_algo_file_index,
                              std::unordered_set<int>& current_blocks,
                              const int requested_num_fps,
                              std::priority_queue<std::pair<std::string, int>>& minhash_fingerprints);

    /**
     * inner function of initializeAppearancesMatrix
     * updates the given initial_volumes_size and optimal_volumes_size
     * also updating the class fields: m_appearances_matrix, m_fingerprint_to_size
     * and m_initial_system_size_with_deduplication
     * @param fingerprints_for_clustering_ordered_by_SN - selected fingerprints ordered by the SN (ascending)
     * @param workload_file_line -  a workload line representing a file
     * @param file_index - the current algo file index
     * @param current_volume - set of blocks we already handled in current volume
     * @param optimal_volumes - set of blocks we already handled in all volumes
     * @param load_balance - is load balancing enabled
     */
    void updateAppearancesMatWithWorkloadFileLine(const std::vector<int>& fingerprints_for_clustering_ordered_by_SN,
                                                  const std::string& workload_file_line,
                                                  const int file_index,
                                                  std::unordered_set<int>& current_volume,
                                                  std::unordered_set<int>& optimal_volumes,
                                                  const bool load_balance);
private:
    /**
     * close the workloads' streams as given in workloads_streams
     * @param workloads_streams - workloads' streams vector
     */
    static void closeWorkloadsStreams(std::vector<std::ifstream>& workloads_streams);

    /**
     * inner function of initializeAndSelectFingerprints
     * updates the given minhash_fingerprints and current_blocks with the given workload_line
     * @param workload_line - a workload line representing a block
     * @param current_blocks - set of blocks we already handled
     * @param requested_num_fps - threshold of number of fingerprints we want to collect
     * @param minhash_fingerprints - priority queue of fingerprints
     */
    static void updateFPWithWorkloadBlockLine(const std::string& workload_line, std::unordered_set<int>& current_blocks,
                                              const int requested_num_fps,
                                              std::priority_queue<std::pair<std::string, int>>& minhash_fingerprints);

    /**
     * inner function of initializeAndSelectFingerprints
     * @param minhash_fingerprints - priority queue of selected fingerprints
     * @return selected fingerprints ordered by the SN (ascending)
     */
    static std::vector<int> getFPForClusteringOrderedBySN(
            std::priority_queue<std::pair<std::string, int>>& minhash_fingerprints);

// m_dissimilarities_matrix - the dissimilarities' matrix
// m_fingerprint_to_size - fingerprint to size mapping
// m_appearances_matrix - the appearances' matrix -> is file x contains fp y
// m_algo_file_index_to_workload_file_index - each spot on vector is a mapping of index x in the clustering to file number in workload
// m_workloads_paths - workload paths vector
// m_requested_number_of_fingerprints - number of fingerprints to use, -1 for 'all'
// m_number_of_files_for_clustering - number of fingerprints to use, -1 for 'all'
// m_number_of_fingerprints_for_clustering - actual number of fingerprints
// m_initial_system_size_with_deduplication - system initial size with deduplication
// m_optimal_system_size_with_deduplication - system optimal size with deduplication
private:
    std::vector<std::vector<DissimilarityCell>> m_dissimilarities_matrix;
    std::vector<int> m_fingerprint_to_size;
    std::vector<std::vector<bool>> m_appearances_matrix;
    std::vector<std::map<int, int>> m_algo_file_index_to_workload_file_index;
    std::vector<std::string> m_workloads_paths;
    const int m_requested_number_of_fingerprints;
    int m_number_of_files_for_clustering;
    int m_number_of_fingerprints_for_clustering;
    int m_initial_system_size_with_deduplication;
    int m_optimal_system_size_with_deduplication;
};
