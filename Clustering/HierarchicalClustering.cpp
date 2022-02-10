/**
Created by Roei Kisous
*/
#include "HierarchicalClustering.hpp"
#include <cmath>
#include <algorithm>
#include <cfloat>
#include <set>
#include <unordered_set>
#include <random>
#include <fstream>
#include <iostream>
#include <utility>
using namespace std;

struct HierarchicalClustering::ClustersMergeOffer final{
public:
    ClustersMergeOffer(double weighted_dissimilarity, int cluster1, int cluster2):
    weighted_dissimilarity(weighted_dissimilarity), cluster1(cluster1), cluster2(cluster2)
    {
    }

public:
    double weighted_dissimilarity;
    int cluster1;
    int cluster2;
};

struct HierarchicalClustering::ClusteringParams final{
public:
    ClusteringParams(
            const double traffic, const int seed, const double gap,
            const int initial_system_size, const int optimal_system_size,
            const bool load_balance, const int number_of_clusters, const int max_results_size,  const int eps) :
            seed(seed),
            gap(gap),
            traffic(traffic),
            load_balance(load_balance),
            number_of_clusters(number_of_clusters),
            max_results_size(max_results_size),
            eps(eps)
    {
        srand(seed);

        const double W_T = traffic / 100.0;
        approx_system_size = initial_system_size - W_T * (initial_system_size - optimal_system_size);
    }
public:
    double traffic;
    double gap;
    double approx_system_size;
    int seed;
    int eps;
    int number_of_clusters;
    int max_results_size;
    bool load_balance;
};

HierarchicalClustering::HierarchicalClustering(
        std::unique_ptr<AlgorithmDSManager>&  ds,
        const vector<double> &lb_sizes):
        m_lb_sizes(lb_sizes),
        m_ds(ds.release())
{
    initClusters();
}

void HierarchicalClustering::initClusters(){
    // if only need to reset
    if(m_clusters.size() == m_ds->getNumberOfFilesForClustering())
    {
        for(const auto& node : m_clusters)
            node->reset();

        return;
    }

    m_clusters = std::vector<std::unique_ptr<Node>>(m_ds->getNumberOfFilesForClustering());

    // create a node for each file with its size
    for (int i = 0; i < m_clusters.size(); ++i) {
        double size = 0;
        for (int j = 0; j < m_ds->getNumberOfFingerprintsForClustering(); ++j) {
            if (m_ds->isFileHasFingerprint(i, j))
                size += m_ds->getFingerprintSize(j);
        }

        m_clusters[i] = std::make_unique<Node>(i, size);
    }
}

HierarchicalClustering::ClustersMergeOffer HierarchicalClustering::findSmallestDissimilarity(
        const ClusteringParams& clustering_params){

    vector<ClustersMergeOffer> sorted_merge_offers;
    sorted_merge_offers.reserve(clustering_params.max_results_size);

    // iterate the lower triangular of the matrix and find our minimum dissimilarity values
    for (int i = 0; i < m_ds->getNumberOfFilesForClustering(); i++) {
        for (int j = 0; j < i; j++) {
            updateSortedMergeOffers(clustering_params, sorted_merge_offers, i, j);
        }
    }

    // in case we did not find anything, sanity check, not suppose to happen unless we use load balancing
    if (sorted_merge_offers.size() == 0)
        return {DBL_MAX, -1, -1};

    // sort them in *ascending* order
    sort(sorted_merge_offers.begin(), sorted_merge_offers.end(), HierarchicalClustering::sortAsc);

    return randClustersMergeOfferFromVector(sorted_merge_offers, clustering_params.gap);
}

void HierarchicalClustering::updateSortedMergeOffers(
        const ClusteringParams& clustering_params,
        vector<ClustersMergeOffer>&  sorted_merge_offers,
        const int cluster1,
        const int cluster2
        )
        {
    // if one of the given cluster is invalid there is no need to update the sorted dissimilarities values
    if(m_ds->getDissimilarityCell(cluster1, cluster2).jaccard_distance == DBL_MAX)
        return;

    const double current_weighted_dissimilarity = getWeightedValueBetweenClusters(clustering_params,cluster1, cluster2);
    // in case we already collect max_results_size values and this value is greater than our max value, we will not
    // consider inserting it to our "best" dissimilarities solution cache
    if (sorted_merge_offers.size() >= clustering_params.max_results_size &&
    current_weighted_dissimilarity >= sorted_merge_offers.front().weighted_dissimilarity)
        return;

    // in case we are using load balance and the new sizes in case of a merge are invalid
    if (clustering_params.load_balance && !isValidClustersMerge(clustering_params, cluster1, cluster2))
        return;


    if (sorted_merge_offers.size() < clustering_params.max_results_size){
        // we are yet to collect max_results_size values
        sorted_merge_offers.emplace_back(current_weighted_dissimilarity, cluster1, cluster2);
    }
    else if (current_weighted_dissimilarity < sorted_merge_offers.front().weighted_dissimilarity){
        // this value is smaller than our max value, replace
        sorted_merge_offers.front() = {current_weighted_dissimilarity, cluster1, cluster2};
    }

    // sort them in descending order for keeping the max dissimilarities at front
    sort(sorted_merge_offers.begin(), sorted_merge_offers.end(), HierarchicalClustering::sortDesc);
}

bool HierarchicalClustering::sortDesc(const ClustersMergeOffer &a, const ClustersMergeOffer &b) {
    if(a.weighted_dissimilarity > b.weighted_dissimilarity)
        return true;

    if(a.weighted_dissimilarity == b.weighted_dissimilarity && a.cluster1 > b.cluster1)
        return true;

    if(a.weighted_dissimilarity == b.weighted_dissimilarity && a.cluster1 == b.cluster1 && a.cluster2 > b.cluster2)
        return true;

    return false;
}

bool HierarchicalClustering::sortAsc(const ClustersMergeOffer &a, const ClustersMergeOffer &b) {
    return sortDesc(b,a);
}

double HierarchicalClustering::getWeightedValueBetweenClusters(const ClusteringParams& clustering_params,
                                                               const int cluster1, const int cluster2){
    const int max_number_of_clusters = m_ds->getNumOfWorkloads();
    const double W_T = clustering_params.traffic / 100.0;
    const AlgorithmDSManager::DissimilarityCell& cluster1_2_dissimilarity = m_ds->getDissimilarityCell(cluster1, cluster2);

    const double physical_distance =
            static_cast<double>(cluster1_2_dissimilarity.origin_clusters.size()) / max_number_of_clusters;

    return (1 - W_T) * physical_distance + W_T * cluster1_2_dissimilarity.jaccard_distance;
}

HierarchicalClustering::ClustersMergeOffer HierarchicalClustering::randClustersMergeOfferFromVector(
        const vector<ClustersMergeOffer>& sorted_merge_offers,
        const double gap){
    // get the last index which is in gap range
    int max_index_in_gap = sorted_merge_offers.size() - 1;
    while (max_index_in_gap >= 1){
        if(sorted_merge_offers[max_index_in_gap].weighted_dissimilarity <=
            sorted_merge_offers.front().weighted_dissimilarity * (1 + gap / 100.0))
            break;

        max_index_in_gap--;
    }

    // fetch random value from the vector
    const int random_index = rand() % (max_index_in_gap + 1);

    return sorted_merge_offers[random_index];
}

bool HierarchicalClustering::isValidClustersMerge(const ClusteringParams& clustering_params, const int cluster1,
                                                  const int cluster2) const{

    // compile a list of all the current clusters' sizes in the system beside those we merge
    vector<double> curr_sizes;
    curr_sizes.reserve(m_ds->getNumberOfFilesForClustering());

    for (int i = 0; i < m_ds->getNumberOfFilesForClustering(); ++i) {
        if (i != cluster1 && i != cluster2 && m_clusters[i]->isActivated())
            curr_sizes.push_back(m_clusters[i]->getSize());
    }

    // add the size of the clusters we merge into one
    curr_sizes.push_back(getMergedClusterSize(cluster1, cluster2));

    // sort in descending order
    std::sort(curr_sizes.begin(), curr_sizes.end(), std::greater<double>());

    // assert that our constraints are satisfied
    for (int i = 0; i < m_lb_sizes.size(); ++i) {
        if (curr_sizes[i] > (m_lb_sizes[i] / 100) * clustering_params.approx_system_size)
            return false;
    }

    return true;
}

double HierarchicalClustering::getClusterSize(const set<int>& files_indices) const{
    // find the files of the merged cluster
    unordered_set<int> cluster_blocks; // cluster's blocks
    int cluster_size = 0; // cluster's size

    for (int fp_index = 0; fp_index < m_ds->getNumberOfFingerprintsForClustering(); ++fp_index) {
        for (const auto& file_index : files_indices) {
            if (m_ds->isFileHasFingerprint(file_index, fp_index) && cluster_blocks.find(fp_index) == cluster_blocks.cend()){
                cluster_size += m_ds->getFingerprintSize(fp_index);
                cluster_blocks.insert(fp_index);
            }
        }
    }

    return cluster_size;
}

double HierarchicalClustering::getMergedClusterSize(const int cluster1, const int cluster2) const{
    // find the files of the merged cluster
    unordered_set<int> cluster_blocks; // cluster's blocks
    int cluster_size = 0; // cluster's size

    for (int fp_index = 0; fp_index < m_ds->getNumberOfFingerprintsForClustering(); ++fp_index) {
        for (const auto& file_index : m_clusters[cluster1]->getCurrentFiles()) {

            if (m_ds->isFileHasFingerprint(file_index, fp_index) && cluster_blocks.find(fp_index) == cluster_blocks.cend()){
                cluster_size += m_ds->getFingerprintSize(fp_index);
                cluster_blocks.insert(fp_index);
            }
        }

        for (const auto& file_index : m_clusters[cluster2]->getCurrentFiles()) {
            if (m_ds->isFileHasFingerprint(file_index, fp_index) && cluster_blocks.find(fp_index) == cluster_blocks.cend()){
                cluster_size += m_ds->getFingerprintSize(fp_index);
                cluster_blocks.insert(fp_index);
            }
        }
    }

    return cluster_size;
}

void HierarchicalClustering::mergeClusters(const ClustersMergeOffer& merge_offer) {
    completeLinkage(merge_offer);

    // disable cluster 2 files to cluster 1 since we gonna it as the merged cluster
    m_clusters[merge_offer.cluster1]->addFiles( m_clusters[merge_offer.cluster2]->getCurrentFiles());

    //calculate cluster 1 size after merging cluster 2 to it
    m_clusters[merge_offer.cluster1]->setSize(getClusterSize(m_clusters[merge_offer.cluster1]->getCurrentFiles()));

    // disable nodes since we gonna user cluster1
    m_clusters[merge_offer.cluster2]->disableNode();

    //deactivate merge_offer.cluster2 since we are going to use merge_offer.cluster1 index only
    m_ds->deactivateClusterInDissimilarityMat(merge_offer.cluster2);
}

void HierarchicalClustering::completeLinkage(const ClustersMergeOffer& merge_offer) {
    // complete linkage before performing the merge_ofer (update dissimilarities matrix)

    // iterate each file
    for (int i = 0; i < m_ds->getNumberOfFilesForClustering(); ++i) {
        const AlgorithmDSManager::DissimilarityCell cluster2_and_i_dissimilarity = m_ds->getDissimilarityCell(merge_offer.cluster2, i);

        if ((cluster2_and_i_dissimilarity.jaccard_distance != DBL_MAX)) {
            const AlgorithmDSManager::DissimilarityCell cluster1_and_i_dissimilarity = m_ds->getDissimilarityCell(merge_offer.cluster1, i);

            // if this file is still active, create the set of origin clusters
            const unordered_set<int> file1_origin_clusters = cluster2_and_i_dissimilarity.origin_clusters;
            const unordered_set<int> file2_origin_clusters = cluster1_and_i_dissimilarity.origin_clusters;

            const AlgorithmDSManager::DissimilarityCell new_dissimilarity_cell = {
                    max(cluster2_and_i_dissimilarity.jaccard_distance, cluster1_and_i_dissimilarity.jaccard_distance),
                    Utility::getUnorderedSetsUnion(file1_origin_clusters, file2_origin_clusters)
            };

            // set the new value ( we only set cluster1 and i since cluster2 will be soon deactivated)
            m_ds->setDissimilarityCell(merge_offer.cluster1, i, new_dissimilarity_cell);
        }
    }
}

string HierarchicalClustering::getResultFileName(const double traffic, const int seed, const double gap,
                                                 const bool load_balance,const int eps, const string& output_path_prefix){
    return output_path_prefix + "_WT" + Utility::getString(traffic/100) +"_S"+ Utility::getString(seed) + "_G" +
        Utility::getString(gap) + "_E" + to_string(eps) + (load_balance? "_lb" : "") + ".csv";
}

void HierarchicalClustering:: run(const bool load_balance,const int eps, const int number_of_clusters,
                                 const vector<double> &traffics, const vector<int> &seeds,
                                 const vector<double> &gaps, const string& output_path_prefix) {

    // for each combo of traffic, seed, gap execute the algorithm
    const vector<set<int>> initial_clusters = m_ds->getInitialClusters();
    const map<string, set<int>> init_cluster_map_sn = getConvertedResults(m_ds->getClustersAsMap(initial_clusters));
    for (const double traffic :traffics){
        for (const int seed : seeds){
            for (const double gap :gaps){
                // create our results csv
                static constexpr int MAX_RESULTS_SIZE = 10;
                ClusteringParams clustering_params(traffic, seed, gap, m_ds->getInitialSystemSize(),
                                                   m_ds->getOptimalSystemSize(), load_balance,
                                                   number_of_clusters, MAX_RESULTS_SIZE, eps);

                const double original_approx_system_size = clustering_params.approx_system_size;

                // loop until we succeed to build a valid dendrogram
                for(int attempt=1; !performClustering(clustering_params); ++attempt){
                    clustering_params.approx_system_size = original_approx_system_size *
                            pow(1 + eps / 100.0, attempt);

                }

                const map<string, set<int>> clustering_result = getClusteringResult(initial_clusters);

                // output results to csv
                outputConvertedResultsToCsv(clustering_params, getConvertedResults(clustering_result),
                                            init_cluster_map_sn, output_path_prefix);

            }
        }
    }
}


map<string, set<int>> HierarchicalClustering::getClusteringResult(const vector<set<int>>& initial_clusters) const {
    const vector<set<int>> final_clusters = getCurrentClustering();

    const vector<int> workload_to_cluster_map = getGreedyWorkloadToClusterMapping(initial_clusters, final_clusters);

    return std::move(getArrangedResult(final_clusters, workload_to_cluster_map));
}

map<string, set<int>> HierarchicalClustering::getArrangedResult(const vector<set<int>>& final_clusters,
                                                                          const vector<int>& workload_to_cluster_map) const{
    map<string, set<int>> arranged_results;
    for (int i = 0; i < workload_to_cluster_map.size(); ++i)
        arranged_results.insert({ m_ds->getWorkloadFullPath(i), final_clusters[workload_to_cluster_map[i]]});

    return std::move(arranged_results);
}

map<string, set<int>> HierarchicalClustering::getConvertedResults(const map<string, set<int>>& clustering) const{
    map<string, set<int>> converted_results;
    for(const auto& workload_to_files_mapping: clustering){
        const string workload_name = workload_to_files_mapping.first;
        set<int> cluster_file_sn;
        for (const int file_index : workload_to_files_mapping.second) {
            cluster_file_sn.insert(m_ds->getFileSN(file_index));
        }
        converted_results.insert({ workload_name, cluster_file_sn});
    }

    return converted_results;
}

vector<int> HierarchicalClustering::getGreedyWorkloadToClusterMapping(const vector<set<int>> &initial_clusters,
                                                                      const vector<set<int>> &final_clusters) const{

    const vector<set<int>> initial_clusters_blocks = getBlocksInClusters(initial_clusters);
    const vector<set<int>> final_clusters_blocks = getBlocksInClusters(final_clusters);

    vector<vector<int>> intersections = getBlocksIntersectionOfClusters(initial_clusters_blocks, final_clusters_blocks);

    // index i is the i'th original workload
    vector<int> workload_to_cluster_map(initial_clusters_blocks.size(), 0);

    // pair the volumes and clusters based on that matrix in a greedy manner
    for (int i = 0; i < initial_clusters_blocks.size(); ++i) {
        const pair<int, int> max_indices = findMax(intersections);

        for (int j = 0; j < final_clusters_blocks.size(); j++)
            intersections[max_indices.first][j] = -1;

        for(int k=0;k<initial_clusters_blocks.size(); ++k)
            intersections[k][max_indices.second] = -1;

        workload_to_cluster_map[max_indices.first] = max_indices.second;
    }

    return workload_to_cluster_map;
}

vector<vector<int>> HierarchicalClustering::getBlocksIntersectionOfClusters(const vector<set<int>>& clusters1_blocks,
                                                                            const vector<set<int>>& clusters2_blocks){
    // init intersection matrix
    vector<vector<int>> intersections;
    intersections.reserve(clusters1_blocks.size());
    for (int i = 0; i < clusters1_blocks.size(); ++i)
        intersections.emplace_back(vector<int>(clusters2_blocks.size(), 0));

    // set each cell to the size of block intersection between the clusters1 to a clusters2
    for (int i = 0; i < clusters1_blocks.size(); ++i) {
        for (int j = 0; j < clusters2_blocks.size(); ++j) {
            intersections[i][j] =
                    count_if(clusters1_blocks[i].cbegin(), clusters1_blocks[i].cend(),
                             [&](int fp_index)
                             {
                        return clusters2_blocks[j].find(fp_index) != clusters2_blocks[j].cend();
                             });
        }
    }

    return intersections;
}

vector<set<int>> HierarchicalClustering::getBlocksInClusters(const vector<set<int>>& clusters) const{
    const int max_number_of_clusters = m_ds->getNumOfWorkloads();

    vector<set<int>> clusters_blocks;
    clusters_blocks.reserve(max_number_of_clusters);

    //find blocks for each cluster in the given clusters
    for(const auto& cluster : clusters){
        set<int> cluster_blocks;
        for(const int file_index : cluster){
            for (int fp_index = 0; fp_index < m_ds->getNumberOfFingerprintsForClustering(); ++fp_index) {
                if (m_ds->isFileHasFingerprint(file_index, fp_index)) // if block exists in file - add it
                    cluster_blocks.insert(fp_index);
            }
        }

        clusters_blocks.emplace_back(cluster_blocks);
    }

    return std::move(clusters_blocks);
}

pair<int, int> HierarchicalClustering::findMax(vector<vector<int>> mat) {
    static constexpr int UNSET_VALUE = -1;

    int maxX=UNSET_VALUE;
    int maxY=UNSET_VALUE;
    int maxElement=UNSET_VALUE;

    // look for the largest intersection
    for (int i = 0; i < mat.size(); ++i) {
        for (int j = 0; j < mat[i].size(); ++j) {
            if (mat[i][j] > maxElement) {
                maxElement = mat[i][j];
                maxX = i;
                maxY = j;
            }
        }
    }

    // return our result
    return {maxX, maxY};
}

vector<set<int>> HierarchicalClustering::getCurrentClustering() const{
    vector<set<int>> result;
    result.reserve(m_ds->getNumberOfFilesForClustering());

    for(const auto& cluster : m_clusters){
        if(cluster->isActivated())
            result.emplace_back(cluster->getCurrentFiles());
    }

    return std::move(result);
}

void HierarchicalClustering::outputConvertedResultsToCsv(const ClusteringParams& clustering_params,
                                                         const map<string, set<int>>& result_mapping,
                                                         const map<string, set<int>>& initial_mapping,
                                                         const string& output_path_prefix){

    ofstream result_csv_file = ofstream(getResultFileName(clustering_params.traffic, clustering_params.seed,
                                                          clustering_params.gap, clustering_params.load_balance,
                                                          clustering_params.eps, output_path_prefix));
    try {
        for(const auto& work_mapping : result_mapping){
            result_csv_file<< work_mapping.first +",";
            int i=1;
            for(const auto& file_sn : initial_mapping.at(work_mapping.first)){
                result_csv_file<< file_sn << (i++==initial_mapping.at(work_mapping.first).size() ? "" : "-");
            }
            result_csv_file<< ",";

            i=1;
            for(const auto& file_sn : work_mapping.second){
                result_csv_file<< file_sn << (i++==work_mapping.second.size() ? "" : "-");
            }
            result_csv_file << ",1" << std::endl;
        }
    }
    catch (...){
        cerr<<"error in outputConvertedResultsToCsv"<<endl;
        result_csv_file.close();
    }

}
bool HierarchicalClustering::performClustering(const ClusteringParams& clustering_params) {
    resetDataStructures();

    for (int iteration = m_ds->getNumberOfFilesForClustering(); iteration > clustering_params.number_of_clusters; iteration--) {
        const ClustersMergeOffer chosen_merge_offer = findSmallestDissimilarity(clustering_params);

        // check if we did not find any suitable clusters to merge, and we are yet to receive number_of_clusters clusters
        if (chosen_merge_offer.weighted_dissimilarity == DBL_MAX)
            return false;

        // merge the clusters we found and update our data structures
        mergeClusters(chosen_merge_offer);
    }

    return true;
}

void HierarchicalClustering::resetDataStructures() {
    m_ds->resetDMBottomTriangle();

    initClusters();
}
