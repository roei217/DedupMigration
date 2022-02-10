/**
Created by Roei Kisous
*/
#include "AlgorithmDSManager.hpp"

#include <iostream>
#include <algorithm>
#include <string>
#include <cfloat>

using namespace std;

AlgorithmDSManager::AlgorithmDSManager(
        const vector<string>& workloads_paths,
        const int requested_number_of_fingerprints,
        const bool load_balance) :
        m_requested_number_of_fingerprints(requested_number_of_fingerprints),
        m_workloads_paths(workloads_paths),
        m_number_of_files_for_clustering(0),
        m_number_of_fingerprints_for_clustering(0),
        m_initial_system_size_with_deduplication(0),
        m_optimal_system_size_with_deduplication(0)
{
    const vector<int> fingerprints_for_clustering_ordered_by_SN = initializeAndSelectFingerprints();
    initializeAppearancesMatrix(fingerprints_for_clustering_ordered_by_SN, load_balance);
    initializeDissimilaritiesMatrix();
}

std::vector<std::ifstream> AlgorithmDSManager::getWorkloadsStreams() const{

    std::vector<std::ifstream> workloads_ifstream;
    //read files
    for (int i = 0; i < m_workloads_paths.size(); ++i) {
        workloads_ifstream.emplace_back(m_workloads_paths[i]);
        if (!workloads_ifstream[i].is_open()) {
            workloads_ifstream[i].close();
            throw std::runtime_error("error opening workload with index " + to_string(i));
        }
    }

    return std::move(workloads_ifstream);
}

void AlgorithmDSManager::closeWorkloadsStreams(std::vector<std::ifstream>& workloads_streams){

    //close all workloads' streams we opened before if still opened
    for(auto& f : workloads_streams){
        if (!f.is_open())
            f.close();
    }
}

void AlgorithmDSManager::updateFPWithWorkloadBlockLine(const std::string& workload_line,
                                                       unordered_set<int>& current_blocks,
                                                       const int requested_num_fps,
                                                       std::priority_queue<std::pair<std::string, int>>& minhash_fingerprints){

    const int block_sn = Utility::getWorkloadBlockSN(workload_line);
    if (!current_blocks.insert(block_sn).second)//already used this fp
        return;

    const string block_fingerprint = Utility::getWorkloadBlockFingerprint(workload_line);
    static const int NO_LIMIT = -1;

    //if we yet to collect m_requested_number_of_fingerprints fps, or we need all the fps, add it
    if (requested_num_fps == NO_LIMIT || minhash_fingerprints.size() < requested_num_fps){
        minhash_fingerprints.push(make_pair(block_fingerprint, block_sn));
        return;
    }

    //if fp is lexi prior to the max value we currently have, replace them
    if (block_fingerprint < minhash_fingerprints.top().first) {
        minhash_fingerprints.pop();
        minhash_fingerprints.push(make_pair(block_fingerprint, block_sn));
    }
}

void AlgorithmDSManager::updateFPWithWorkload(std::ifstream& workload_file, const int workload_index,
                                              int& current_algo_file_index,unordered_set<int>& current_blocks,
                                              const int requested_num_fps,
                                              std::priority_queue<std::pair<std::string, int>>& minhash_fingerprints){
    string line_content;
    while (getline(workload_file, line_content)) {
        if(Utility::isWorkloadBlockLine(line_content)){
            updateFPWithWorkloadBlockLine(line_content, current_blocks, requested_num_fps, minhash_fingerprints);
        }
        else if (Utility::isWorkloadFileLine(line_content)) {
            //count files and match workload i and file m_number_of_files_for_clustering (algorithm's sn) to original file sn
            m_algo_file_index_to_workload_file_index[workload_index][current_algo_file_index++] =
                    Utility::getWorkloadFileSN(line_content);
        }
    }
}

vector<int> AlgorithmDSManager::getFPForClusteringOrderedBySN(priority_queue<pair<string, int>>& minhash_fingerprints){
    vector<int> fingerprints_for_clustering_ordered_by_SN;
    fingerprints_for_clustering_ordered_by_SN.reserve(minhash_fingerprints.size());

    while (!minhash_fingerprints.empty()){
        fingerprints_for_clustering_ordered_by_SN.emplace_back(minhash_fingerprints.top().second);
        minhash_fingerprints.pop();
    }

    //sort by SN
    sort(fingerprints_for_clustering_ordered_by_SN.begin(), fingerprints_for_clustering_ordered_by_SN.end());

    return std::move(fingerprints_for_clustering_ordered_by_SN);
}

vector<int> AlgorithmDSManager::initializeAndSelectFingerprints() {
    // current_blocks for not using same fp twice
    unordered_set<int> current_blocks;
    std::vector<std::ifstream> workloads_streams = getWorkloadsStreams();
    priority_queue<pair<string, int>> minhash_fingerprints; // priority queue of fps' <fp, sn>

    for (int workload_index = 0; workload_index < workloads_streams.size(); ++workload_index) {
        m_algo_file_index_to_workload_file_index.emplace_back();
        updateFPWithWorkload(workloads_streams[workload_index], workload_index, m_number_of_files_for_clustering,
                             current_blocks, m_requested_number_of_fingerprints, minhash_fingerprints);
    }

    vector<int> fingerprints_for_clustering_ordered_by_SN = getFPForClusteringOrderedBySN(minhash_fingerprints);
    m_number_of_fingerprints_for_clustering = fingerprints_for_clustering_ordered_by_SN.size();

    closeWorkloadsStreams(workloads_streams);

    return fingerprints_for_clustering_ordered_by_SN;
}

void AlgorithmDSManager::clearAppearancesMatrix(){
    //create the appearances_matrix with initial values false
    m_appearances_matrix = vector<vector<bool>>();
    m_appearances_matrix.reserve(m_number_of_files_for_clustering);
    for (int i = 0; i < m_number_of_files_for_clustering; ++i)
        m_appearances_matrix.emplace_back(std::vector<bool>(m_number_of_fingerprints_for_clustering, false));
}

void AlgorithmDSManager::updateAppearancesMatWithWorkloadFileLine(
        const vector<int>& fingerprints_for_clustering_ordered_by_SN,
        const string& workload_file_line,
        const int file_index,
        unordered_set<int>& current_volume,
        unordered_set<int>& optimal_volumes,
        const bool load_balance){
    const int number_of_blocks_in_file_line = Utility::getNumOfBlockInFile(workload_file_line);

    //for each fp for this file if the fp was selected, add to relevant data structures
    for (int block_index = 0; block_index < number_of_blocks_in_file_line; ++block_index) {
        const int fingerprint_index = Utility::binarySearch(
                fingerprints_for_clustering_ordered_by_SN,
                0,
                fingerprints_for_clustering_ordered_by_SN.size() - 1,
                Utility::getBlockSNInFileByIndex(workload_file_line, block_index));

        static constexpr int FINGERPRINT_NOT_FOUND = -1;
        if (fingerprint_index != FINGERPRINT_NOT_FOUND) {
            const int block_size= Utility::getBlockSizeInFileByIndex(workload_file_line, block_index);
            m_appearances_matrix[file_index][fingerprint_index] = true;
            m_fingerprint_to_size[fingerprint_index] = block_size;
            if (load_balance) {
                //new fp for specific volume
                if (current_volume.insert(fingerprint_index).second)
                    m_initial_system_size_with_deduplication += block_size;

                //new fp for entire system
                if (optimal_volumes.insert(fingerprint_index).second)
                    m_optimal_system_size_with_deduplication += block_size;
            }
        }
    }
}

void AlgorithmDSManager::initializeAppearancesMatrix(
        const vector<int>& fingerprints_for_clustering_ordered_by_SN,
        bool load_balance) {

    clearAppearancesMatrix();

    //initialize fields
    m_fingerprint_to_size = vector<int>(m_number_of_fingerprints_for_clustering, 0);
    m_initial_system_size_with_deduplication = 0;
    m_optimal_system_size_with_deduplication = 0;

    unordered_set<int> optimal_volumes_size;//fps when all fps in same volume
    int file_index = 0;
    std::vector<std::ifstream> workloads_streams = getWorkloadsStreams();
    for (int i = 0; i < workloads_streams.size(); ++i) {
        //fps for specific initial volume
        unordered_set<int> current_volume;

        string workload_file_line;
        while (getline(workloads_streams[i], workload_file_line)) {
            if(!Utility::isWorkloadFileLine(workload_file_line))
                continue;

            //if we got here it is file line
            updateAppearancesMatWithWorkloadFileLine(fingerprints_for_clustering_ordered_by_SN, workload_file_line,
                                                     file_index, current_volume, optimal_volumes_size,
                                                     load_balance);

            file_index++;
        }
    }

    closeWorkloadsStreams(workloads_streams);
}

void AlgorithmDSManager::initializeDissimilaritiesMatrix() {
    unordered_map<int, int> initial_clusters = getInitialFileToClusterMapping();

    //create the dissimilarity_matrix, do it row by row
    m_dissimilarities_matrix = vector<vector<DissimilarityCell>>();
    m_dissimilarities_matrix.reserve(m_number_of_files_for_clustering);

    for (int i = 0; i < m_number_of_files_for_clustering; ++i) {
        m_dissimilarities_matrix.emplace_back(m_number_of_files_for_clustering);

        // fill the diagonal - irrelevant
        m_dissimilarities_matrix[i][i] = {0, unordered_set<int>({initial_clusters[i]})};

        for (int j = i + 1; j < m_number_of_files_for_clustering; ++j) {
            const float distance = Utility::getJaccardDistance(m_appearances_matrix[i],m_appearances_matrix[j]);
            const unordered_set<int> clusters_of_files = {initial_clusters[i], initial_clusters[j]};
            m_dissimilarities_matrix[i][j] = {distance, clusters_of_files};
        }
    }

    resetDMBottomTriangle();
}

unordered_map<int,int> AlgorithmDSManager::getInitialFileToClusterMapping(){
    //files' algorithm sn to cluster mapping
    unordered_map<int, int> initial_clusters;

    for(int cluster_index=0 ; cluster_index < m_algo_file_index_to_workload_file_index.size(); ++cluster_index){
        for(const auto& file_mapping : m_algo_file_index_to_workload_file_index[cluster_index]){
            initial_clusters[file_mapping.first] = cluster_index;
        }
    }

    return std::move(initial_clusters);
}

void AlgorithmDSManager::resetDMBottomTriangle(){
    // fill the bottom half
    for (int i = 0; i < m_dissimilarities_matrix.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            m_dissimilarities_matrix[i][j] = m_dissimilarities_matrix[j][i];
        }
    }
}

AlgorithmDSManager::DissimilarityCell& AlgorithmDSManager::getDissimilarityCell(const int cluster1, const int cluster2){

    const int column=min(cluster1,cluster2);
    const int row=max(cluster1,cluster2);

    return m_dissimilarities_matrix[row][column];
}

void AlgorithmDSManager::setDissimilarityCell(const int cluster1, const int cluster2,
                                              const AlgorithmDSManager::DissimilarityCell& dissimilarity_cell){
    const int column=min(cluster1,cluster2);
    const int row=max(cluster1,cluster2);

    m_dissimilarities_matrix[row][column] = dissimilarity_cell;
}

void AlgorithmDSManager::deactivateClusterInDissimilarityMat(const int cluster_index){
    for (int i = 0; i < m_dissimilarities_matrix.size(); ++i)
        setDissimilarityCell(cluster_index, i, {DBL_MAX, unordered_set<int>()});
}

int AlgorithmDSManager::getNumberOfFilesForClustering() const{
    return m_number_of_files_for_clustering;
}

int AlgorithmDSManager::getNumberOfFingerprintsForClustering() const{
    return m_number_of_fingerprints_for_clustering;
}

int AlgorithmDSManager::getFingerprintSize(const int fp_index) const{
    return m_fingerprint_to_size[fp_index];
}

bool AlgorithmDSManager::isFileHasFingerprint(const int file_index, const int fp_index) const{
    return m_appearances_matrix[file_index][fp_index];
}

int AlgorithmDSManager::getFileSN(const int file_index) const{
    for(const map<int, int>& cluster_mapping : m_algo_file_index_to_workload_file_index){
        if(cluster_mapping.find(file_index) != cluster_mapping.cend())
            return cluster_mapping.at(file_index);
    }

    // should not arrive here since every file must exist exactly in one cluster
    exit(EXIT_FAILURE);
}

vector<set<int>> AlgorithmDSManager::getInitialClusters() const{
    // retrieve the initial clusters' files
    vector<set<int>> initial_clusters;
    for(const auto& cluster_index_mapping: m_algo_file_index_to_workload_file_index){
        set<int> current_cluster;
        for (const auto& file_number_mapping : cluster_index_mapping)
            current_cluster.insert(file_number_mapping.first);

        initial_clusters.push_back(current_cluster);
    }

    return std::move(initial_clusters);
}

map<string, set<int>> AlgorithmDSManager::getClustersAsMap(const vector<set<int>>& clusters_as_vector) const{
    map<string, set<int>> initial_clusters_as_map;
    for(int workload_index = 0; workload_index < clusters_as_vector.size(); ++workload_index)
        initial_clusters_as_map[getWorkloadFullPath(workload_index)] = clusters_as_vector[workload_index];

    return initial_clusters_as_map;
}

string AlgorithmDSManager::getWorkloadFullPath(const int workload_index) const{


    return Utility::getFullPath(m_workloads_paths[workload_index]);
}

int AlgorithmDSManager::getNumOfWorkloads() const{
    return m_workloads_paths.size();
}

int AlgorithmDSManager::getInitialSystemSize() const{
    return m_initial_system_size_with_deduplication;
}

int AlgorithmDSManager::getOptimalSystemSize() const{
    return m_optimal_system_size_with_deduplication;
}