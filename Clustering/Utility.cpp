#include "Utility.hpp"

#include <sstream>
#include <sys/stat.h>
#include <limits.h>
#include <stdlib.h>
#include <exception>
#include <iomanip>
#include <sstream>

using namespace std;

int Utility::binarySearch(const std::vector<int>& array, int left_index, int right_index, int value) {
    while (left_index <= right_index) {
        int middle_index = left_index + (right_index - left_index) / 2;
        // check if x is here
        if (array[middle_index] == value)
            return middle_index;
        // if x is greater, ignore left
        if (array[middle_index] < value)
            left_index = middle_index + 1;
        // if x is smaller, ignore right
        else
            right_index = middle_index - 1;
    }
    // if we reach here, we did not find anything
    return -1;
}

vector<string> Utility::splitString(const string &line, const char delimiter) {
    stringstream ss(line);
    string item;
    vector<string> elements;
    while (getline(ss, item, delimiter)) {
        elements.push_back(move(item));
    }

    return std::move(elements);
}

bool Utility::isFileExists(const string& path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

float Utility::getJaccardDistance(const vector<bool>& row1,const vector<bool>& row2) {
    int co_absence = 0;
    int xor_exists = 0;

    for (int k = 0; k < row1.size(); ++k) {
        co_absence += ((!row1[k]) && (!row2[k]))? 1 : 0;
        xor_exists += (((!row1[k]) && (row2[k])) || ((row1[k]) && (!row2[k]))) ? 1: 0;
    }

    const int diff = row1.size() - co_absence;
    // in case we divide by 0, set value to 0
    if(diff == 0)
        return 0;

    return static_cast<float>(xor_exists) / static_cast<float>(diff);
}

bool Utility::isWorkloadBlockLine(const std::string& workload_line){
    const vector<string> splitted_line_content = Utility::splitString(workload_line, ',');

    static const std::string BLOCK_LINE_START_STRING = "B";
    return splitted_line_content.front() == BLOCK_LINE_START_STRING;
}

bool Utility::isWorkloadFileLine(const std::string& workload_line){
    const vector<string> splitted_line_content = Utility::splitString(workload_line, ',');

    static const std::string FILE_LINE_START_STRING = "F";

    return splitted_line_content.front() == FILE_LINE_START_STRING;
}

int Utility::getNumOfBlockInFile(const std::string &workload_line) {
    const vector<string> splitted_line_content = Utility::splitString(workload_line, ',');

    static constexpr int FILE_NUM_OF_BLOCKS_INDEX = 4;
    return stoi(splitted_line_content[FILE_NUM_OF_BLOCKS_INDEX]);
}

int Utility::getBlockSNInFileByIndex(const std::string &workload_line, const int block_index) {
    const vector<string> splitted_line_content = Utility::splitString(workload_line, ',');

    static constexpr int FILE_RECIPE_BLOCK_SN_START_INDEX = 5;
    static constexpr int BLOCK_SN_JUMP = 2;
    return stoi(splitted_line_content[FILE_RECIPE_BLOCK_SN_START_INDEX + block_index * BLOCK_SN_JUMP]);
}

int Utility::getBlockSizeInFileByIndex(const std::string &workload_line, const int block_index) {
    const vector<string> splitted_line_content = Utility::splitString(workload_line, ',');

    static constexpr int FILE_RECIPE_BLOCK_SIZE_START_INDEX = 6;
    static constexpr int BLOCK_SIZE_JUMP = 2;
    return stoi(splitted_line_content[FILE_RECIPE_BLOCK_SIZE_START_INDEX + block_index * BLOCK_SIZE_JUMP]);
}

int Utility::getWorkloadBlockSN(const std::string &workload_line) {
    const vector<string> splitted_line_content = Utility::splitString(workload_line, ',');

    static constexpr int BLOCK_LINE_SN_INDEX = 1;
    return stoi(splitted_line_content[BLOCK_LINE_SN_INDEX]);
}

std::string Utility::getWorkloadBlockFingerprint(const std::string &workload_line) {
    vector<string> splitted_line_content = Utility::splitString(workload_line, ',');

    static constexpr int BLOCK_LINE_FP_INDEX = 2;
    return std::move(splitted_line_content[BLOCK_LINE_FP_INDEX]);
}

int Utility::getWorkloadFileSN(const std::string &workload_line) {
    const vector<string> splitted_line_content = Utility::splitString(workload_line, ',');

    static constexpr int FILE_LINE_SN_INDEX = 1;
    return stoi(splitted_line_content[FILE_LINE_SN_INDEX]);
}

std::string Utility::getFullPath(const string &path) {
    char resolved_path[PATH_MAX];
    if(!realpath(path.c_str(), resolved_path))
        throw std::runtime_error("getFullPath failed");
    return std::string(resolved_path);
}

std::string Utility::getString(const double val) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << val;
    return stream.str();
}
