/**
Created by Roei Kisous
*/
#pragma once

#include <string>
#include <map>
#include <vector>
#include <set>
#include <unordered_set>

namespace Utility{
    /**
     * simple binary search for a sorted array
     * @param array - array (as vector) to search
     * @param left_index - current left bound
     * @param right_index - current right bound
     * @param value - value to find
     * @return - index of value or -1 if does not exist
     */
    int binarySearch(const std::vector<int>& array, int left_index, int right_index, int value);

    /**
     * splits a string to tokens
     * @param line- string to split
     * @param delimiter - char to tokenize by
     * @return a vector with the tokens
     */
    std::vector<std::string> splitString(const std::string &line, char delimiter);

    /**
     * splits a string to tokens
     * @param path- absolute or relative path
     * @return the full path of the given path
     */
    std::string getFullPath(const std::string &path);

    /**
     *
     * @param path - path to be checked
     * @return true if file path exists, false otherwise
     */
    bool isFileExists(const std::string& path) ;

    /**
     *
     * @param val - value
     * @return string of the given val only with precision of 2 after dot
     */
    std::string getString(const double val) ;

    /**
     *
     * @param set1 - set 1 in the sets union
     * @param set2 - set 2 in the sets union
     * @return the union set of set1 + set2
     */
    template <class T>
    std::set<T> getSetsUnion(const std::set<T>& set1, const std::set<T>& set2)
    {
        std::set<T> union_set = set1;
        union_set.insert(set2.cbegin(), set2.cend());

        return std::move(union_set);
    }

    /**
     *
     * @param set1 - set 1 in the sets union
     * @param set2 - set 2 in the sets union
     * @return the union set of set1 + set2
     */
    template <class T>
    std::unordered_set<T> getUnorderedSetsUnion(const std::unordered_set<T>& set1, const std::unordered_set<T>& set2)
    {
        std::unordered_set<T> union_set = set1;
        union_set.insert(set2.cbegin(), set2.cend());

        return std::move(union_set);
    }

    /**
     * simple function to calculate the dissimilarity between two rows (with same lengths)
     * @param row1 - the first row to compare
     * @param row2 - the second row to compare
     * @return the Jaccard distance
     */
    float getJaccardDistance(const std::vector<bool>& row1,const std::vector<bool>& row2);

    /**
     * @param workload_line - line from a workload csv file
     * @return whether the line describes a block
     */
    bool isWorkloadBlockLine(const std::string& workload_line);

    /**
     * @param workload_line - line from a workload csv file
     * @return whether the line describes a file
     */
    bool isWorkloadFileLine(const std::string& workload_line);

    /**
     * @param workload_line - line from a workload csv file
     * @return return the workload's block sn as described in workload_line
     */
    int getWorkloadBlockSN(const std::string& workload_line);

    /**
     * @param workload_line - line from a workload csv file
     * @return return the workload's file sn as described in workload_line
     */
    int getWorkloadFileSN(const std::string& workload_line);

    /**
     * @param workload_line - a line from a workload csv file describing a file
     * @return return the num of blocks in file as described in workload_line
     */
    int getNumOfBlockInFile(const std::string &workload_line);

    /**
     * @param workload_line - line from a workload csv file  describing a file
     * @param block_index - index of the block in the file recipe
     * @return return the block's SN
     */
    int getBlockSNInFileByIndex(const std::string &workload_line, const int block_index);

    /**
    * @param workload_line - line from a workload csv file  describing a file
    * @param block_index - index of the block in the file recipe
    * @return return the block's size
    */
    int getBlockSizeInFileByIndex(const std::string &workload_line, const int block_index);

    /**
     * @param workload_line - line from a workload csv file
     * @return return the workload's block fingerprint as described in workload_line
     */
    std::string getWorkloadBlockFingerprint(const std::string &workload_line);
}
