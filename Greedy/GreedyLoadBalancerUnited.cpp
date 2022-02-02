#include <algorithm>
#include <numeric>
#include <cfloat>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <chrono>
#include <set>
#include <boost/algorithm/string.hpp>
#include <unordered_map> 
#include <cmath>  

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    Data structure explanation:
    _____________________________________
    std::vector<std::string> source_volume_list:
        A vector containing the paths of all source volumes.

    std::vector<std::string> target_volume_list:
        A vector containing the paths of all target volumes.

    std::vector<double> desired_balance:
        A vector containing the desired_balance of all volumes. given fractions can be between 0 to 1 and must sum up to a total of 1.

    std::vector<std::vector<short>> block_volume_sourceRefCount:
        The data structure saving for each source volume how many references it has for a block. IE how many files are using it.
        If there are 4 files using block 10 in volume 3. then block_volume_sourceRefCount[10][3] == 4

    std::vector<std::vector<short>> block_volume_targetRefCount:
        The data structure saving for each target volume how many references it has for a block. IE how many files are using it.
        If there are 4 files using block 10 in volume 3. then block_volume_targetRefCount[10][3] == 4
        
    std::vector<std::vector<std::string>> volume_file_fileList:
        The data structure saving for each file it's representation in the Michal and Paulina standard. 

    std::pair<int, int> lastSourceSn_block_file:
        A pair of the last SN in the source volume of blocks and of files

    std::pair<int, int> lastTargetSn_block_file:
        A pair of the last SN in the target volume of blocks and of files

    std::vector<double> blockSizes:
        The data structure saving every size of the blocks. 
        blockSizes[10] is the size in KB of the block with SN 10.

    std::vector<double> volumeSizes:
        The data structure saving the current size of all of the volumes.
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    Global variables and their default values
*/

bool fullMigration = true;
std::vector<int> loopEvade(2, -2);
double totalSourceSize = 0; 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    Utility general cpp functions
*/

void updateLoopEvade(int fileSn) {
    loopEvade.insert(loopEvade.begin(), fileSn);
    loopEvade.pop_back();
}

bool isLoop(int fileSn) {
    return std::find(loopEvade.begin(), loopEvade.end(), fileSn) != loopEvade.end();
}

/**
 * @brief Returns the indexes of the sorted array from the given array
 *        [1,2,3,4] -> [0,1,2,3]
 *        [4,2,3,1] -> [3,1,2,0]
 * @param v Is the array from which we want the sorted indexes
 */
template <typename T>
std::vector<std::size_t> sort_indexes(const std::vector<T> &v) {

  // Initialize original index locations
  std::vector<std::size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // Sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](std::size_t i1, std::size_t i2) {return v[i1] < v[i2];});

  return idx;
}

/**
 * @brief Converts a string to a bool value
 * @param v The string to convert
 */
bool string2bool (const std::string & v) {
    return !v.empty () &&
        (strcasecmp (v.c_str (), "true") == 0 ||
         atoi (v.c_str ()) != 0);
}

/**
 * @brief Splits a string to an array by a delimiter
 * @param str The string
 * @param delimiter The delimiter
 */
std::vector<std::string> split_string(std::string str, const std::string &delimiter)
{
    std::vector<std::string> result;
    boost::split(result, str, boost::is_any_of(delimiter));
    return result;
}

/**
 * @brief Replaces all occurrences of 'from' to 'to' in the string 
 * @param str The string
 * @param from The string to be replaced
 * @param to The string to be replaced with
 */
std::string replaceAll(std::string& str, const std::string& from, const std::string& to) {
    std::string returnString = str;
    if(from.empty())
        return returnString;
    size_t start_pos = 0;
    while((start_pos = returnString.find(from, start_pos)) != std::string::npos) {
        returnString.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
    return returnString;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    Utility problem specific functions
*/

/**
 * @brief Calculates the volume size
 * @param volumeSizes Reference where the result would be put
 * @param block_volume_targetRefCount The block reference count object
 * @param blockSizes An array with the sizes of all the blocks
 */
void getVolumeSizes(std::vector<double> &volumeSizes, std::vector<std::vector<short>> &block_volume_targetRefCount, std::vector<double> &blockSizes) {
    for(int volume = 0; volume < volumeSizes.size(); volume++) {
        volumeSizes[volume] = 0.0;
        for(int block = 0; block < blockSizes.size(); block++) {
            if(block_volume_targetRefCount[block][volume] > 0) {
                volumeSizes[volume] += blockSizes[block];
            }
        }
    }
}

/**
 * @brief Calculate each volume's percentage of the total size
 * @param volumeSizes Reference where the result would be put
 * @param block_volume_targetRefCount The block reference count object
 * @param blockSizes An array with the sizes of all the blocks
 */
void getVolumePercentages(std::vector<double> &volumeSizes, std::vector<std::vector<short>> &block_volume_targetRefCount, std::vector<double> &blockSizes) {
    double totalSize = 0;
    for(int volume = 0; volume < volumeSizes.size(); volume++) {
        volumeSizes[volume] = 0.0;
        for(int block = 0; block < blockSizes.size(); block++) {
            if(block_volume_targetRefCount[block][volume] > 0) {
                volumeSizes[volume] += blockSizes[block];
            }
        }
        totalSize += volumeSizes[volume];
    }

    for(int volume = 0; volume < volumeSizes.size(); volume++) {
        volumeSizes[volume] = volumeSizes[volume] / totalSize;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    Utility ingestion phase functions
*/

/**
 * @brief Calculate each volume's percentage of the total size
 * @param blockSizes An array where the block sizes would be put
 * @param block_volume_targetRefCount The block reference count object
 */
void calcBlockSizes(std::vector<double> &blockSizes, std::vector<std::string> &volumeList) {
    std::string content;
    std::vector<std::string> splitted_content;
    for (auto &volume : volumeList) {
        std::ifstream volume_stream(volume.c_str(), std::ifstream::in);
        if (!volume_stream.is_open())
        {
            std::cout << "error opening volume file - calcBlockSizes" << volume << std::endl;
            exit(1);
        }

        while (std::getline(volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "F")
            {
                int number_of_blocks_in_file_line = std::stoi(splitted_content[4]);
                for (register int i = 0; i < 2 * number_of_blocks_in_file_line; i += 2) //Read block_sn and block_size simultaneously and add constrains to the model.
                {
                    int blockSn = std::stoi(splitted_content[5 + i]);
                    int size_read = std::stoi(splitted_content[6 + i]);
                    if(blockSizes[blockSn]  == 0.0) {
                        double size = ((double)size_read) / 1024.0;
                        blockSizes[blockSn] = size;
                    }
                }
            }
        }    
    }
}

/**
 * @brief Find the last block index and last file index and return a pair of <lastBlock, lastFile>
 * @param volumeList A list of all the volume paths
 */
std::pair<int, int> getLastBlockAndFileSn(std::vector<std::string> &volumeList) 
{
    int lastBlockSn = 0;
    int lastFileSn = 0;
    for (auto &source_volume : volumeList) {
        std::ifstream source_volume_stream(source_volume.c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            std::cout << "error opening volume file - getLastBlockSn" << source_volume << std::endl;
            exit(1);
        }
        
        std::string content;
        std::vector<std::string> splitted_content;

        while (std::getline(source_volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "F")
            {
                lastFileSn = std::max(lastFileSn, std::stoi(splitted_content[1]));
            }
            if (splitted_content[0] == "B")
            {
                lastBlockSn = std::max(lastBlockSn, std::stoi(splitted_content[1]));
            }
        }
        source_volume_stream.close();
	}   
    return std::make_pair(lastBlockSn, lastFileSn);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    Utility optimization phase functions
*/

/**
 * @brief Check if the provided solution object is an empty solution object
 * solution is either marked by an illegal sourceVolume of -1 or is de facto empty by moving nothing
 * @param solution A tuple representing a transaction: (sourceVolume, targetVolume, fileSn, best replicate, best delete, best move, csvLineOfMovedFile)
 */

bool isEmptySolution(std::tuple<int, int, int, double, double, double, std::string> &solution) {
    
    if(std::get<0>(solution) == -1) {
        return true;
    }
    // return ((std::get<3>(solution) == 0) && (std::get<4>(solution) == 0) && (std::get<5>(solution) == 0));
}

/**
 * @brief Check if the given state is a legal state
 * a legal state is when all the volumes' relative size to the total size is within the given margin
 * @param volumeSizes The current size of the volumes
 * @param desired_balance The desired balance we want in the system
 * @param margin The error margin we allow ourselves from the desired balance for each volume
 */
bool isLegalState(std::vector<double> &volumeSizes, std::vector<double> &desired_balance, double margin) {
    double totalSize = 0.0;
    for(int volume = 0; volume < volumeSizes.size(); volume++) {
        totalSize += volumeSizes[volume];
    }
    for(int volume = 0; volume < volumeSizes.size(); volume++) {
        double load = volumeSizes[volume] / totalSize;
        if((load > desired_balance[volume] + margin) || (load < desired_balance[volume] - margin)) {
            return false;
        }
    }
    return true;
}

/*
    bestTransfer: (fileSn, replicated, deleted, movement, csv line of the best file)
*/
/**
 * @brief Check if a given transfer would result in a legal state
 * a legal state is when all the volumes' relative size to the total size is within the given margin
 * a given bestTransfer isBalanced if simulating the transfer would result in a legal state
 * 
 * @param volumeSizes The current size of the volumes
 * @param bestTransfer: the transfer given to ve checked (fileSn, replicated, deleted, movement, csv line of the best file)
 * @param source: Index of source volume
 * @param target: Index of target volume
 * @param desired_balance The desired balance we want in the system
 * @param margin The error margin we allow ourselves from the desired balance for each volume
 */
bool isBalanced(std::vector<double> &volumeSizes, std::tuple<int, double, double, double, std::string> bestTransfer, int source, int target, std::vector<double> &desired_balance, double margin) {
    std::vector<double> copy(volumeSizes);
    double replicated = std::get<1>(bestTransfer);
    double deleted   = std::get<2>(bestTransfer);
    double movement   = std::get<3>(bestTransfer);
    copy[source] -= deleted;
    copy[target] += movement + replicated;
    return isLegalState(copy, desired_balance, margin);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    Core logic function
*/

/**
 * @brief Parse the input and initialize all the needed databases with input
 * a seemingly superfluous distinction between source and target is made, this is due to the need to allow the code to regress back to previous models
 * in which not all volumes were both source and target volumes.
 * @param source_volume_list A list of all source volumes
 * @param target_volume_list A list of all target volumes
 * @param block_volume_sourceRefCount Reference to the block reference object for source volume
 * @param block_volume_targetRefCount Reference to the block reference object for target volume
 * @param lastSourceSn_block_file A tuple of the last index of the blocks and files in the source volumes
 * @param lastTargetSn_block_file A tuple of the last index of the blocks and files in the target volumes
 * @param volume_file_fileList A list saving the unparsed lines of the volumes as are saved in the input csv files
 */
void ingest(std::vector<std::string> &source_volume_list, std::vector<std::string> &target_volume_list, std::vector<std::vector<short>> &block_volume_sourceRefCount,
 std::vector<std::vector<short>> &block_volume_targetRefCount, std::pair<int, int> &lastSourceSn_block_file,
  std::pair<int, int> &lastTargetSn_block_file, std::vector<std::vector<std::string>> &volume_file_fileList) {
    
    int lastBlock = std::max(lastSourceSn_block_file.first, lastTargetSn_block_file.first);
    
    for(int i = 0; i < lastBlock + 1; i++) {
        std::vector<short> volumeList = std::vector<short> (source_volume_list.size(), 0);
        block_volume_sourceRefCount.push_back(volumeList);
    }
    for(int i = 0; i < lastBlock + 1; i++) {
        std::vector<short> volumeList = std::vector<short> (target_volume_list.size(), 0);
        block_volume_targetRefCount.push_back(volumeList);
    }   

    for (int source = 0; source < source_volume_list.size(); source++) {
        //File line only relevant for source volumes since target volumes cant move files
        std::vector<std::string> fileLine;

        std::ifstream source_volume_stream(source_volume_list[source].c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            std::cout << "error opening volume file - addConstraint_RemmapedFileHasAllItsBlocks" << source_volume_list[source] << std::endl;
            exit(1);
        }

        std::string content;
        std::vector<std::string> splitted_content;
        while (std::getline(source_volume_stream, content)) {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "F") {
                fileLine.push_back(content);
            }
            if (splitted_content[0] == "B") {
                int blockSn = std::stoi(splitted_content[1]);
                int refCount = std::stoi(splitted_content[3]);
                block_volume_sourceRefCount[blockSn][source] = refCount;
            }
        }
        source_volume_stream.close();
        volume_file_fileList.push_back(fileLine);
	}   

    for (int target = 0; target < target_volume_list.size(); target++) {
        std::ifstream target_volume_stream(target_volume_list[target].c_str(), std::ifstream::in);
        if (!target_volume_stream.is_open())
        {
            std::cout << "error opening volume file - addConstraint_RemmapedFileHasAllItsBlocks" << target_volume_list[target] << std::endl;
            exit(1);
        }

        std::string content;
        std::vector<std::string> splitted_content;
        while (std::getline(target_volume_stream, content)) {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "B") {
                int blockSn = std::stoi(splitted_content[1]);
                int refCount = std::stoi(splitted_content[3]);
                block_volume_targetRefCount[blockSn][target] = refCount;
            }
        }
        target_volume_stream.close();
	}       
}

/**
 * @brief Find the best file to transfer between a source and target volume
 * retVal: (fileSn, replicated, deleted, movement, csv line of the best file)
 * @param source Source volume index
 * @param target Target volume index
 * @param block_volume_sourceRefCount Reference to the block reference object for source volume
 * @param block_volume_targetRefCount Reference to the block reference object for target volume
 * @param currentTraffic How much traffic had allready been used
 * @param trafficSize How much traffic in total can be used
 * @param file_fileList All files
 * @param checkLegalState A configuration, if true we make sure the transfer is legal within margin, else we allow illegal transfers. Illegal transfers are necessary 
 * for the balancing phase since the system needs several transfers before ending in a legal state
 * @param volumeSizes A list of the current volumeSizes
 * @param desired_balancethe Desired relative sizes of the volumes
 * @param margin Error margin from the desired relative size
 */
std::tuple<int, double, double, double, std::string> getBestTransfer(int source, int target,
    std::vector<std::vector<short>> &block_volume_sourceRefCount,
    std::vector<std::vector<short>> &block_volume_targetRefCount,
    double currentTraffic, double trafficSize,
    std::vector<std::string> &file_fileList,
    bool checkLegalState,
    std::vector<double> &volumeSizes,
    std::vector<double> &desired_balance,
    double margin) {

    //Create empty transfer
    std::tuple<int, double, double, double, std::string> retVal;
    std::get<0>(retVal) = -1;       //Best fileSn
    std::get<1>(retVal) = DBL_MAX;  //Best totalReplication
    std::get<2>(retVal) = -DBL_MAX;  //Best totalDeletion
    std::get<3>(retVal) = DBL_MAX;  //Best movement
    std::get<4>(retVal) = "";       //Csv line of the best file
    if(fullMigration && source == target) {
        return retVal;
    }

    std::vector<std::string> splitted_content;
    
    for (auto &fileLine : file_fileList) {
        splitted_content = split_string(fileLine, ",");
        if (splitted_content[0] == "F")
        {
            int fileSn = std::stoi(splitted_content[1]);
            if(isLoop(fileSn)) {
                continue;
            }
            int number_of_blocks_in_file_line = std::stoi(splitted_content[4]);
            double replication = 0;
            double deletion = 0;
            double movement = 0;
            for (register int i = 0; i < 2 * number_of_blocks_in_file_line; i += 2) //Read block_sn and block_size simultaneously and add constrains to the model.
            {
                int block_sn = std::stoi(splitted_content[5 + i]);
                int size_read = std::stoi(splitted_content[6 + i]);
                double blockSize = ((double)size_read) / 1024.0;
                int refInSource = block_volume_sourceRefCount[block_sn][source];
                int refInTarget = block_volume_targetRefCount[block_sn][target];
                if(refInSource == 1) {
                    if(refInTarget > 0) {
                        deletion += blockSize; // Can delete
                    } else {
                        movement += blockSize; // Needs to move
                    }
                } else {
                    if(refInTarget == 0) {
                        replication += blockSize; //Must replicate
                        movement += blockSize; // Needs to move
                    }
                }
            }
            // if((deletion == 0 && movement == 0) || (deletion == 0 && replication == 0)) {
            if((checkLegalState && deletion == 0) ||
                (!checkLegalState && (deletion == 0 && movement == 0 && replication == 0))) {
                continue;
            }
            double bestReclaim = std::get<1>(retVal) / std::max(1.0, std::get<2>(retVal));
            double currentReclaim = replication / std::max(1.0, deletion);
            if((currentTraffic + movement <= trafficSize) && (currentReclaim < bestReclaim)) {
                std::tuple<int, double, double, double, std::string> currentTransfer;            
                std::get<0>(currentTransfer) = fileSn;
                std::get<1>(currentTransfer) = replication;
                std::get<2>(currentTransfer) = deletion;
                std::get<3>(currentTransfer) = movement;
                std::get<4>(currentTransfer) = fileLine;
                if(!checkLegalState || isBalanced(volumeSizes, currentTransfer, source ,target, desired_balance, margin)) {
                    std::get<0>(retVal) = fileSn;
                    std::get<1>(retVal) = replication;
                    std::get<2>(retVal) = deletion;
                    std::get<3>(retVal) = movement;
                    std::get<4>(retVal) = fileLine;
                }
            }
        }
    }
    return retVal;
}

/**
 * @brief Perform a single greedy iteration and return the greediest iteration
 * retVal: (sourceVolume, targetVolume, fileSn, replicate, delete, move, csv line of the best file to move)
 * @param source_volume_list List of source volumes
 * @param target_volume_list List of target volumes
 * @param block_volume_sourceRefCount Reference to the block reference object for source volume
 * @param block_volume_targetRefCount Reference to the block reference object for target volume
 * @param currentTraffic How much traffic had allready been used
 * @param trafficSize How much traffic in total can be used
 * @param volume_file_fileList List saving the files of each volume as string lines
 * @param makeLegal A configuration, if true then the iteration would attempt to move the system to a legal state, if false then the iteration
 * would assume it's in a legal state and try to optimize the system without leaving the legal state
 * @param volumeSizes A list of the current volumeSizes
 * @param desired_balancethe Desired relative sizes of the volumes
 * @param margin Error margin from the desired relative size
 */
std::tuple<int, int, int, double, double, double, std::string> greedyIterate(std::vector<std::string> &source_volume_list,
    std::vector<std::string> &target_volume_list,
    std::vector<std::vector<short>> &block_volume_sourceRefCount,
    std::vector<std::vector<short>> &block_volume_targetRefCount,
    double currentTraffic, double trafficSize,
    std::vector<std::vector<std::string>> &volume_file_fileList,
    bool makeLegal,
    std::vector<double> &volumeSizes,
    std::vector<double> &desired_balance,
    double margin
) {
    // Create an empty solution
    std::tuple<int, int, int, double, double, double, std::string> retVal;
    std::get<0>(retVal) = -1;       //Best sourceVolume
    std::get<1>(retVal) = -1;       //Best targetVolume
    std::get<2>(retVal) = -1;       //Best fileSn
    std::get<3>(retVal) = DBL_MAX;  //Best replicate
    std::get<4>(retVal) = -DBL_MAX; //Best delete
    std::get<5>(retVal) = DBL_MAX;  //Best move
    std::get<6>(retVal) = "";       //Csv line of the best file to move

    //Capacity-reduction step
    if(!makeLegal) {
        for(int source = 0; source < source_volume_list.size(); source++) {
            for(int target = 0; target < target_volume_list.size(); target++) {
                std::tuple<int, double, double, double, std::string> bestTransfer;
                bestTransfer = getBestTransfer(source, target, block_volume_sourceRefCount, block_volume_targetRefCount,
                    currentTraffic, trafficSize,
                    volume_file_fileList[source],
                    true, volumeSizes, desired_balance, margin);
                double bestReclaim = std::get<3>(retVal) / std::max(1.0, std::get<4>(retVal));
                double currentReclaim = std::get<1>(bestTransfer) / std::max(1.0, std::get<2>(bestTransfer));
                if( currentTraffic + std::get<3>(bestTransfer) <= trafficSize && currentReclaim < bestReclaim ) {
                    if(isBalanced(volumeSizes, bestTransfer, source ,target, desired_balance, margin)) {
                        std::get<0>(retVal) = source;
                        std::get<1>(retVal) = target;
                        std::get<2>(retVal) = std::get<0>(bestTransfer);
                        std::get<3>(retVal) = std::get<1>(bestTransfer);
                        std::get<4>(retVal) = std::get<2>(bestTransfer);
                        std::get<5>(retVal) = std::get<3>(bestTransfer);
                        std::get<6>(retVal) = std::get<4>(bestTransfer);
                    }
                }
            }
        }
    }

    //Load-balancing step
    if(makeLegal) {
        // Heuristically to achieve balance we try to move from the biggest volume to the smallest volume
        // this heuristic would need to be changed should the assumption of a desired equal balance be broken
        // it would still work, just be extremely inefficient if not changed
        std::vector<std::size_t> sorted_indexes = sort_indexes(volumeSizes);

        int bigVolume_index = sorted_indexes.size() - 1;
        int smallVolume_index = 0;

        std::tuple<int, double, double, double, std::string> bestTransfer = getBestTransfer(sorted_indexes[bigVolume_index], sorted_indexes[smallVolume_index], block_volume_sourceRefCount, block_volume_targetRefCount,
            currentTraffic, trafficSize,
            volume_file_fileList[sorted_indexes[bigVolume_index]],
            false, volumeSizes, desired_balance, margin);

        while((std::get<0>(bestTransfer) == -1) && bigVolume_index > 0) {
            smallVolume_index++;
            if(smallVolume_index >= sorted_indexes.size()) {
                bigVolume_index--;
                smallVolume_index = 0;
            }
            if(bigVolume_index < 0) {
                break;
            }
            if(bigVolume_index == smallVolume_index) {
                continue;
            }
            bestTransfer = getBestTransfer(sorted_indexes[bigVolume_index], sorted_indexes[smallVolume_index], block_volume_sourceRefCount, block_volume_targetRefCount, 
                currentTraffic, trafficSize, 
                volume_file_fileList[sorted_indexes[bigVolume_index]],
                false, volumeSizes, desired_balance, margin);
        }
        if(std::get<0>(bestTransfer) == -1) {
            std::get<0>(retVal) = -1;
            std::get<1>(retVal) = -1;          
        } else {
            std::get<0>(retVal) = sorted_indexes[bigVolume_index];
            std::get<1>(retVal) = sorted_indexes[smallVolume_index];
        }
        std::get<2>(retVal) = std::get<0>(bestTransfer);
        std::get<3>(retVal) = std::get<1>(bestTransfer);
        std::get<4>(retVal) = std::get<2>(bestTransfer);
        std::get<5>(retVal) = std::get<3>(bestTransfer);
        std::get<6>(retVal) = std::get<4>(bestTransfer);
    }

    return retVal;
}

/**
 * @brief Receives a transfer object and applies the transfer to all of the objects saving the system state
 * @param bestIterationTransfer The transfer to apply
 * (sourceVolume, targetVolume, fileSn, replicate, delete, move, csv line of the best file to move)
 * @param currentTraffic How much traffic had allready been used
 * @param maxTraffic How much traffic in total can be used
 * @param block_volume_sourceRefCount Reference to the block reference object for source volume
 * @param block_volume_targetRefCount Reference to the block reference object for target volume
 * @param volume_file_fileList List saving the files of each volume as string lines
 */
void updateDbs(std::tuple<int, int, int, double, double, double, std::string> &bestIterationTransfer, double currentTraffic, double maxTraffic, 
    std::vector<std::vector<short>> &block_volume_sourceRefCount, std::vector<std::vector<short>> &block_volume_targetRefCount,
    std::vector<std::vector<std::string>> &volume_file_fileList) {
    std::vector<std::string> splitted_content = split_string(std::get<6>(bestIterationTransfer), ",");
    int fileSn = std::stoi(splitted_content[1]);
    updateLoopEvade(fileSn);
    int number_of_blocks_in_file_line = std::stoi(splitted_content[4]);

    int source = std::get<0>(bestIterationTransfer);
    int target = std::get<1>(bestIterationTransfer);

    for (register int i = 0; i < 2 * number_of_blocks_in_file_line; i += 2) //Read block_sn and block_size simultaneously and add constrains to the model.
    {
        int block_sn = std::stoi(splitted_content[5 + i]);
        block_volume_sourceRefCount[block_sn][source]--;
        block_volume_targetRefCount[block_sn][target]++;
        if(fullMigration) {
            block_volume_targetRefCount[block_sn][source]--;
            block_volume_sourceRefCount[block_sn][target]++;
        }
    }

    auto itr = std::find(volume_file_fileList[std::get<0>(bestIterationTransfer)].begin(), volume_file_fileList[std::get<0>(bestIterationTransfer)].end(), std::get<6>(bestIterationTransfer));
    if (itr != volume_file_fileList[std::get<0>(bestIterationTransfer)].end()) {
        volume_file_fileList[std::get<0>(bestIterationTransfer)].erase(itr);
    }
    //If we are under fullMigration then target volumes are also source volumes and so we need to give them the moved files as they can move again
    if(fullMigration) {
        volume_file_fileList[std::get<1>(bestIterationTransfer)].push_back(std::get<6>(bestIterationTransfer));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    Logging and output
*/

/**
 * @brief Print the chosen iteration into a csv solution file
 * @param outputCsvFile Solution file
 * @param iteration Iteration number
 * @param bestIterationTransfer The transfer chosen for the iteration
 * bestIterationTransfer: (sourceVolume, targetVolume, fileSn, replicate, delete, move, csv line of the best file to move)
 * @param currentTraffic Total traffic used
 * @param currentMigration Total data deleted
 * @param time Time that passed
 * @param iterationType Either balancing or optimization
 */
void printIteration(std::ofstream &outputCsvFile, int iteration, std::tuple<int, int, int, double, double, double, std::string> &bestIterationTransfer, double currentTraffic,
 double currentMigration, double time, std::string iterationType) {
    double trafficPercent = ((double)currentTraffic / totalSourceSize) * 100;
    double migrationPercent = ((double)currentMigration / totalSourceSize) * 100;

    outputCsvFile << iteration << "," << std::get<0>(bestIterationTransfer) << "," << std::get<1>(bestIterationTransfer) << "," << std::get<2>(bestIterationTransfer) << "," << std::get<3>(bestIterationTransfer) << "," << std::get<4>(bestIterationTransfer) << "," << std::get<5>(bestIterationTransfer) <<"," << currentTraffic << "," << currentMigration<< "," << trafficPercent << "," << migrationPercent << "," << time << "," << iterationType << "\n";
}

/**
 * @brief Saves input after the solution, for debugging purposes
 * @param write_solution Solution file
 * @param volume_list Input list of the volumes
 * @param T_percentage Input traffic percentage
 * @param Margin Input error margin
 * @param model_time_limit Input time limit
 */
void saveInput(std::string write_solution, std::string volume_list, double T_percentage, double Margin, int model_time_limit) {
    std::ofstream solution(write_solution, std::ios_base::app);
    if (!solution)
    {
        std::cout << "Cannot open output file" << write_solution << std::endl;
        exit(1);
    }
    solution << "volume_list: " << volume_list << std::endl;
    solution << "T_percentage: " << T_percentage<< std::endl;
    solution << "Margin: " << Margin << std::endl;
    solution << "model_time_limit: " << model_time_limit << std::endl;
}

/**
 * @brief Saves total run conclusion line, for statistic purposes
 * @param conclusionFilePath Path to the conclusion file
 * @param volume_list Input list of the volumes
 * @param maxTraffic Input traffic percentage
 * @param margin Input error margin
 * @param ingestionTime How long the ingestion took
 * @param totalTime Total runtime 
 * @param totalTraffic Total traffic used
 * @param totalMigration Total deletion 
 * @param initialVolumePercentage Volume state of input
 * @param finalVolumePercentage Volume state after migration plan
 */
void makeConclusionLine(std::string conclusionFilePath, std::string volume_list, int maxTraffic, double margin, double ingestionTime, double totalTime, 
    double totalTraffic, double totalMigration, std::vector<double> &initialVolumePercentage, std::vector<double> &finalVolumePercentage) {
    std::ofstream out(conclusionFilePath, std::ios_base::app);
    if (!out)
    {
        std::cout << "Cannot open output file\n";
    }

    auto splitSlashes = split_string(volume_list, "/");
    std::vector<std::string> splitted_content = split_string(replaceAll(splitSlashes[splitSlashes.size() - 1], "_ck8_", "_"), "_");
    
    std::string volumeName = splitted_content[0];
    std::string k = replaceAll(splitted_content[1], "k", "");
    std::string volumeNumber = "5";
    if(splitted_content.size() > 2) {
        volumeNumber = splitted_content[2];
    }
    double OptimizationTime = totalTime - ingestionTime;

    out << "GREEDY load balance" << ", "
        << volumeName << ", "
        << volumeNumber << ", "
        << k << ", "
        << maxTraffic << ", "
        << margin << ", "
        << ingestionTime << ", "
        << OptimizationTime << ","
        << totalTraffic << ", "
        << totalMigration << ", " ;
        
    for(int i = 0; i < initialVolumePercentage.size(); i++) {
        out << initialVolumePercentage[i] << ", "
            << finalVolumePercentage[i] << ", " ;
    }
    out << std::endl;
    out.close();
}

/**
 * @brief Saves the data required for creating the migration plan file and prints it into the given file
 * @param migrationPlanFilePath Path to where the migration plan should be saved
 * @param source_volume_list the data structures saving the paths of all the source volumes
 * @param volume_file_fileList the data structures of the current file location
 * @param isInitilization we want to save the file locations at the start and at the end. if false it means we are at the end and should print
 */            
std::vector<std::vector<std::string> > parsedCSV;
void saveMigrationPlanState(std::string migrationPlanFilePath, std::vector<std::string> source_volume_list, std::vector<std::vector<std::string>> volume_file_fileList, bool isInitilization) {
    std::cout << "saveMigrationPlanState " << isInitilization << std::endl;
    try 
    {
        if(isInitilization) 
        {
            for(int i = 0; i < source_volume_list.size(); i++) {
                parsedCSV.push_back(std::vector<std::string>());
                parsedCSV[i].push_back(source_volume_list[i]);
                if(volume_file_fileList[i].size()) {
                    std::string fileList = split_string(volume_file_fileList[i][0], ",")[1];
                    if(volume_file_fileList[i].size() > 1) {
                        for(int j = 1; j < volume_file_fileList[i].size(); j++) {
                            std::vector<std::string> splitted_content = split_string(volume_file_fileList[i][j], ",");
                            fileList += "-" + splitted_content[1];
                        }
                    }
                    parsedCSV[i].push_back(fileList);
                } else {
                    parsedCSV[i].push_back("");
                }
            } 
        } else {
            std::ofstream out(migrationPlanFilePath, std::ios_base::trunc);
            if (!out)
            {
                std::cout << "Cannot open migration plan output file\n";
            } 
            for(int i = 0; i < source_volume_list.size(); i++) {
                if(volume_file_fileList[i].size()) {
                    std::string fileList = split_string(volume_file_fileList[i][0], ",")[1];
                    if(volume_file_fileList[i].size() > 1) {
                        for(int j = 1; j < volume_file_fileList[i].size(); j++) {
                            std::vector<std::string> splitted_content = split_string(volume_file_fileList[i][j], ",");
                            fileList += "-" + splitted_content[1];
                        }
                    }
                    parsedCSV[i].push_back(fileList);
                } else {
                    parsedCSV[i].push_back("");
                }
                parsedCSV[i].push_back("1");
            }                           
            for(int i = 0; i < parsedCSV.size(); i++) {
                for(int j = 0; j < parsedCSV[i].size(); j++) {
                    out << parsedCSV[i][j] << ",";
                }
                out << std::endl;
            }            
        }
    } catch(const std::exception& e) {
        std::cout << e.what();
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    /*
        Read configurations
    */
    const auto begin = std::chrono::high_resolution_clock::now(); //Start the stopwatch for the total time.
    if (argc != 7)                       
    {
        std::cout
            << "arguments format is: {volumelist} {output} {conclusionFile} {timelimit} {Traffic} {margin}"
            << std::endl;
        return 0;
    }

    std::string volume_list = std::string(argv[1]);
    std::string output = std::string(argv[2]);
    std::string conclusionFile = std::string(argv[3]);
    int timelimit_inSeconds = std::stoi(std::string(argv[4]));
    int trafficPercent = std::stoi(std::string(argv[5]));
    double margin = std::stod(std::string(argv[6]));
    // The current model has no differentiation between source ant target volumes, but for regression we leave the configuration here
    // Making this false would make source volumes only capable of sending files to target volumes
    fullMigration = true; 

    /*
        Initialize data structures
    */
    std::string content;
    std::vector<std::string> source_volume_list;
    std::vector<std::string> target_volume_list;
    std::vector<double> desired_balance;

    std::ifstream volume_list_f(volume_list.c_str(), std::ifstream::in);
    if (!volume_list_f.is_open())
    {
        std::cout << "error opening volume list." << std::endl;
        exit(1);
    }

    /*
        Parse volume strings and desired balances
        The input files are made of lines made of 
        [path to volume file in Michal & Paulina standard], [a|s|t where s is source, t is target and a is all], [fraction of desired balance]
    */
    while (std::getline(volume_list_f, content)) 
    {
        std::size_t firstComma = content.find(", ");
        std::size_t secondComma = content.find(", ", firstComma + 1);

        std::string file_path = content.substr(0, firstComma);
        std::string volume_type = content.substr(firstComma + 1, secondComma - firstComma - 1);
        std::string::iterator end_pos = std::remove(volume_type.begin(), volume_type.end(), ' ');
        volume_type.erase(end_pos, volume_type.end());

        double desiredLoad = std::stod(content.substr(secondComma + 1));
        if(volume_type == "s" || volume_type == "a") {
            source_volume_list.push_back(file_path);
        }
        if(volume_type == "t" || volume_type == "a") {
            target_volume_list.push_back(file_path);            
            desired_balance.push_back(desiredLoad);
        }
    }

    /*
        Validate correct load distribution
    */
    double total = 0;
    for(int volume = 0; volume < desired_balance.size(); volume++) {
        total += desired_balance[volume];
    }

    // if(total != 1) {
    //     std::cout << "illegal load distribution! " << total << " != 1" << std::endl;
    // }

    /*
        Ingest data
    */
    std::vector<std::vector<short>> block_volume_sourceRefCount;
    std::vector<std::vector<short>> block_volume_targetRefCount;
    std::vector<std::vector<std::string>> volume_file_fileList;
    
    std::pair<int, int> lastSourceSn_block_file = getLastBlockAndFileSn(source_volume_list);
    std::pair<int, int> lastTargetSn_block_file = getLastBlockAndFileSn(target_volume_list);
    std::vector<double> blockSizes (lastSourceSn_block_file.first + 1, 0.0);
    std::vector<double> volumeSizes (target_volume_list.size(), 0.0);

    calcBlockSizes(blockSizes, source_volume_list);
    ingest(source_volume_list, target_volume_list, block_volume_sourceRefCount, block_volume_targetRefCount, lastSourceSn_block_file, lastTargetSn_block_file, volume_file_fileList);
    getVolumeSizes(volumeSizes, block_volume_targetRefCount, blockSizes);
    totalSourceSize = 0;
    for(int volume = 0; volume < volumeSizes.size(); volume++) {
        totalSourceSize += volumeSizes[volume];
    }
    std::vector<double> startingVolumePercentages (target_volume_list.size(), 0.0);
    getVolumePercentages(startingVolumePercentages, block_volume_targetRefCount, blockSizes);

    double trafficSize = ((double)trafficPercent * totalSourceSize) / 100;
    // Greedy seldom uses it's traffic due to superfluous movements, using only around 83.33 of it's given traffic
    // So we allow it a relative extra 20% as to allow it to use all it's traffic
    double handicappedTraffic = trafficSize + 0.2 * trafficSize;
    double marginSize = ((double)margin * totalSourceSize);
    
    std::cout << "totalSourceSize: " << totalSourceSize << std::endl;
    std::cout << "max Traffic: " << handicappedTraffic << std::endl;
    std::cout << "margin: " << marginSize << std::endl;

    /*
        Prepare first iteration
    */
    double previousTraffic = -1;
    double currentTraffic = 0;
    double currentMigration = 0;
    int iteration = 0;
    std::ofstream outputCsvFile;
    outputCsvFile.open (output);
    outputCsvFile << "Iteration,sourceVolume,targetVolume,fileSn,replicated KB,deleted KB,moved KB,totalTraffic KB,totalMigration KB,traffic% of Total source,deletion% of Total source,Time(s),iteration type\n";
    
    double elapsed_secs = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - begin).count();
    double ingestionTime = elapsed_secs;
    
    saveMigrationPlanState(output + "_migrationPlan", source_volume_list, volume_file_fileList, true);
    /*
        Empirically we found that optimizing to the final margin is worse than slowly refining the margin
        empirically we found the best results for 5 phases starting from 1.5 * margin and ending in 1 * margin
        we also found that an extension smaller than 0.05 results in worse results
        for every margin we first make sure we are in a legal state under the margin and then we optimize under the margin
        margin + marginExtension; margin + 0.75 * marginExtension; margin + 0.5 * marginExtension; margin + 0.25 * marginExtension; margin
    */
    bool skipBalancing = false;
    double marginExtension = std::max(0.05, 0.5 * margin);
    /*
        Receiving a margin of above 1 indicates we wish no load balancing and can skip the strenuous load balancing logic
    */
    if(margin >= 1) {
        marginExtension = 0;
        skipBalancing = true;
    }
    for(double currentMargin = margin + marginExtension; currentMargin >= margin; currentMargin -= (marginExtension / 4)) {
        double maxIterationTrafic = currentTraffic + (handicappedTraffic - currentTraffic) / 5;
        if(currentMargin == margin) {
            maxIterationTrafic = handicappedTraffic;
        }
        getVolumeSizes(volumeSizes, block_volume_targetRefCount, blockSizes);

        std::tuple<int, int, int, double, double, double, std::string> bestIterationTransfer = greedyIterate(source_volume_list, target_volume_list, block_volume_sourceRefCount, block_volume_targetRefCount, currentTraffic, handicappedTraffic, volume_file_fileList, true, volumeSizes, desired_balance, currentMargin);
        /*
            Load-balancing step
        */
        while(!skipBalancing && !isEmptySolution(bestIterationTransfer) && !isLegalState(volumeSizes, desired_balance, currentMargin)) {
            iteration++;
            currentTraffic += std::get<5>(bestIterationTransfer);
            currentMigration += std::get<4>(bestIterationTransfer);
        
            elapsed_secs = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - begin).count();
            printIteration(outputCsvFile, iteration, bestIterationTransfer, currentTraffic, currentMigration, elapsed_secs, "balance");
            updateDbs(bestIterationTransfer, currentTraffic, maxIterationTrafic, block_volume_sourceRefCount, block_volume_targetRefCount, volume_file_fileList);
            getVolumeSizes(volumeSizes, block_volume_targetRefCount, blockSizes);

            bestIterationTransfer = greedyIterate(source_volume_list, target_volume_list, block_volume_sourceRefCount, block_volume_targetRefCount, currentTraffic, maxIterationTrafic, volume_file_fileList, true, volumeSizes, desired_balance, currentMargin);
        }
        if(!skipBalancing) {
            if( !isLegalState(volumeSizes, desired_balance, currentMargin)) {
                std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!cannot balance to " << currentMargin << " !!!!!!!!!!!!!!!!!!" << std::endl;
            } else {
                std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!balanced to " << currentMargin << " !!!!!!!!!!!!!!!!!!" << std::endl;
            }
        }

        double totalSize = 0.0;
        for(int volume = 0; volume < volumeSizes.size(); volume++) {
            totalSize += volumeSizes[volume];
        }

        getVolumeSizes(volumeSizes, block_volume_targetRefCount, blockSizes);
        bestIterationTransfer = greedyIterate(source_volume_list, target_volume_list, block_volume_sourceRefCount, block_volume_targetRefCount, currentTraffic, maxIterationTrafic, volume_file_fileList, false, volumeSizes, desired_balance, currentMargin);
        /*
            Capacity-reduction step
        */
        while(!isEmptySolution(bestIterationTransfer)) {
            iteration++;
            currentTraffic += std::get<5>(bestIterationTransfer);
            currentMigration += std::get<4>(bestIterationTransfer);
            
            elapsed_secs = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - begin).count();
            printIteration(outputCsvFile, iteration, bestIterationTransfer, currentTraffic, currentMigration, elapsed_secs, "optimize");
            if(elapsed_secs > timelimit_inSeconds) {
                break;
            }
            updateDbs(bestIterationTransfer, currentTraffic, maxIterationTrafic, block_volume_sourceRefCount, block_volume_targetRefCount, volume_file_fileList);
            getVolumeSizes(volumeSizes, block_volume_targetRefCount, blockSizes);

            bestIterationTransfer = greedyIterate(source_volume_list, target_volume_list, block_volume_sourceRefCount, block_volume_targetRefCount, currentTraffic, maxIterationTrafic, volume_file_fileList, false, volumeSizes, desired_balance, currentMargin);
        }
        if(isLegalState(volumeSizes, desired_balance, currentMargin)) {
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!optimized " << currentMargin << " !!!!!!!!!!!!!!!!!!" << std::endl;
        } else {
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!unable to optimize " << currentMargin << " !!!!!!!!!!!!!!!!!!" << std::endl;
        }
        if(skipBalancing) {
            break;
        }
    }
    
    std::vector<double> finalVolumePercentages (target_volume_list.size(), 0.0);
    getVolumePercentages(finalVolumePercentages, block_volume_targetRefCount, blockSizes);
    saveMigrationPlanState(output + "_migrationPlan", source_volume_list, volume_file_fileList, false);
    /*
        Save statistics
    */
    makeConclusionLine(conclusionFile, volume_list, trafficPercent, margin, ingestionTime, elapsed_secs, currentTraffic / totalSourceSize, currentMigration / totalSourceSize, startingVolumePercentages, finalVolumePercentages);
    saveInput(output, volume_list, trafficPercent, margin, timelimit_inSeconds);

}
