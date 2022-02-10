/**
Created by Ariel Kolikant
*/
#include "gurobi_c++.h"
#include <fstream>
#include <chrono>
#include <set>
#include <boost/algorithm/string.hpp>
#include <unordered_map> 
#include <unistd.h>
#include <bitset>
#include <sstream>

#define UNDEFINED_STRING "-1"
#define UNDEFINED_DOUBLE -1.0
#define UNDEFINED_INT -1
#define UNDEFINED_STATUS "UNDEFINED_STATUS"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * Data structure explanation:
 * ______________________________________________
    std::vector<std::string> source_volume_list:
        A vector containing the paths of all source volumes.

    std::vector<std::string> target_volume_list:
        A vector containing the paths of all target volumes.

 *  std::vector<double> load_volume_list:
        A vector containing the desired_balance of all volumes. given fractions can be between 0 to 1 and must sum up to a total of 1.

 *  std::vector<std::vector<std::vector<int>>> intersects_source_target_blocksn;
        The data structure of all the intersects between source and target volumes. 
        intersects_source_target_blocksn[1][2] contains a vector with all block SN inside of the intersect between source volume 1 and target volume 2

 *  std::vector<std::pair<int, int>> num_of_blocks_and_files_sourceVolumes:
        A vector of pairs with the number of blocks and the number of files in each source volume

 *  std::vector<std::pair<int, int>> num_of_blocks_and_files_targetVolumes:
        A vector of pairs with the number of blocks and the number of files in each target volume

 *  std::pair<int, int> lastSourceSn_block_file:
        A pair of the last SN in the source volume of blocks and of files

 *  std::vector<double> block_sizes:
        The data structure saving every size of the blocks. 
        blockSizes[10] is the size in KB of the block with SN 10.
 
 * Gurobi data structure explanation:
 * ______________________________________________
 *  std::vector<std::vector<GRBVar*>> C_i_s_t:
 *      Saves the parameter for copying block i from volume s to volume t.
 *      C_i_s_t[0][1][2] is true if block 0 is copied from volume 1 to volume 2 and false if it isn't
 * 
 *  std::vector<GRBVar*> D_i_s:
 *      Saves the parameter for deleting block i from volume s
 *      D_i_s[0][1] is true if block 0 is deleted from volume 1 and false if it isn't
 * 
 *  std::vector<std::vector<GRBVar*>> X_l_s_t;
 *      Saves the parameter for moving file l from volume s to volume t.
 *      X_l_s_t[0][1][2] is true if file 0 is moved from volume 1 to volume 2 and false if it isn't
 * 
 **/

/**
 * Config params as global variables
 **/
int T_percentage = UNDEFINED_DOUBLE;              //The % we want to migrate to an empty destination.
double SOURCE_SIZE = UNDEFINED_DOUBLE;
double Margin = UNDEFINED_DOUBLE;
int InnerFilterSize = UNDEFINED_DOUBLE;
std::string conclusionFileName = UNDEFINED_STRING; //File we save our benchmarks, every line will be a different migration plan summary.
std::string input_file_name = UNDEFINED_STRING;      //Name of the input file, contains all the information needed.
std::string volume_list = UNDEFINED_STRING;      //Name of the input file, contains all the information needed.
double model_time_limit = UNDEFINED_DOUBLE;          //Time limit for the solver.
std::string seed = UNDEFINED_STRING;                 //Seed for the solver.
std::string number_of_threads = UNDEFINED_STRING;    //Number of threads we restrict our solver to run with.

/**
 * ---------------------------------------------------------------------------------------------------------------------------------------------
 * Statistic params for conclusion lines
 **/
double Kbytes_to_replicate = UNDEFINED_DOUBLE;       //KB to replicated as a result from the migration plan.
long int num_of_blocks = UNDEFINED_INT;              //Number of blocks in the input file
double actual_M_percentage = UNDEFINED_DOUBLE;         //The % of physical blocks we decided to migrate.
double actual_M_Kbytes = UNDEFINED_DOUBLE;           //Number of blocks we decided to migrate.
double actual_R_percentage = UNDEFINED_DOUBLE;         //The % of physical blocks we decided to replicate.
double actual_traffic_Percentage = UNDEFINED_DOUBLE;         //The % of physical blocks we decided to replicate.
double actual_deletion_Percentage = UNDEFINED_DOUBLE;         //The % of physical blocks we decided to replicate.
double actual_volume_clean_percentage = UNDEFINED_DOUBLE;         //The % of physical blocks we decided to replicate.
double actual_volume_add_percentage = UNDEFINED_DOUBLE;         //The % of physical blocks we decided to replicate.
double actual_volume_change_percentage = UNDEFINED_DOUBLE;         //The % of physical blocks we decided to replicate.
double total_traffic_and_clean_percentage = UNDEFINED_DOUBLE;         //The % of physical blocks we decided to replicate.
double actual_R_Kbytes = UNDEFINED_DOUBLE;           //Number of blocks we decided to replicate.
int num_of_files = UNDEFINED_INT;                    //Number of files in our input.
std::string solution_status = UNDEFINED_STATUS;      //Solver status at the end of the optimization.
double total_block_size_Kbytes = 0;                  //The total physical size of the system in KB units.
int variable_number = UNDEFINED_INT;
int constraint_number = UNDEFINED_INT;
int source_block_number = UNDEFINED_INT;
int target_block_number = UNDEFINED_INT;
/**
 * ---------------------------------------------------------------------------------------------------------------------------------------------
 * Util functions
 **/

/**
 * @brief Translates a char to its equivalent bit representation
 * @param c The char
 */
const char* hex_char_to_bin(char c)
{
    switch(toupper(c))
    {
        case '0': return "0000";
        case '1': return "0001";
        case '2': return "0010";
        case '3': return "0011";
        case '4': return "0100";
        case '5': return "0101";
        case '6': return "0110";
        case '7': return "0111";
        case '8': return "1000";
        case '9': return "1001";
        case 'A': return "1010";
        case 'B': return "1011";
        case 'C': return "1100";
        case 'D': return "1101";
        case 'E': return "1110";
        case 'F': return "1111";
        default : throw("bad char");
    }
}

/**
 * @brief Splits str according to delimiter, the result strings are stored in a std::vector.
 * 
 * @param str String to split.
 * @param delimiter Split according to delimiter.
 * @return std::vector<std::string> The split std::string in a std::vector.
 */
std::vector<std::string> split_string(std::string str, const std::string &delimiter)
{
    std::vector<std::string> result;
    boost::split(result, str, boost::is_any_of(delimiter));
    return result;
}

/**
 * @brief Replaces all instances of a substring with a different substring
 * 
 * @param str Base string
 * @param from Substring to find
 * @param to Substring to replace with
 * @return String after the replacement
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

/**
 * ---------------------------------------------------------------------------------------------------------------------------------------------
 * Problem specific util functions
 **/

/**
 * @brief Counts the number of metadata lines in the input file.
 * Metadata line is defined to have # symbol at the start of the line.
 * 
 * @param input_file_name The input file name.
 * @return Number of metadata lines in the input file
 */
int get_num_of_metadata_lines(std::string &input_file_name)
{
    std::ifstream f(input_file_name.c_str(), std::ifstream::in);
    int counter = 0;
    std::string content;
    if (!f.is_open())
    {
        std::cout << "error opening file: " << input_file_name << std::endl;
        exit(1);
    }
    std::getline(f, content);
    while (content[0] == '#')
    {
        counter++;
        std::getline(f, content);
    }
    f.close();
    return counter;
}

/**
 * @brief Retrieves the number of files and blocks in the input.
 * metadata lines are in format: "# <type_of_information>:<value>"
 * @param f Reference to the input stream
 * @param num_of_metadata_lines The number of metadata lines
 */
std::pair<int, int> get_num_of_blocks_and_files(std::string f, int num_of_metadata_lines)
{
    std::ifstream f_stream(f.c_str(), std::ifstream::in);
    if (!f_stream.is_open())
    {
        std::cout << "error opening volume file - get_num_of_blocks_and_files" << f << std::endl;
        exit(1);
    }

    const std::string type_of_info_file = "# Num files";
    const std::string type_of_info_block = "# Num Blocks";
    std::string content;
    std::string number_as_string;
    std::string type_of_info;
    bool set_num_files = false, set_num_blocks = false;

    for (int i = 0; i < num_of_metadata_lines; i++)
    {
        std::getline(f_stream, content);
        type_of_info = content.substr(0, content.find(": "));
        if (type_of_info == type_of_info_file)
        {
            num_of_files = std::stoi(content.substr(2 + content.find(": "))); //Sets global variable
            set_num_files = true;
        }
        if (type_of_info == type_of_info_block)
        {
            num_of_blocks = std::stol(content.substr(2 + content.find(": "))); //Sets global variable
            source_block_number = num_of_blocks;
            set_num_blocks = true;
        }
    }
    if (!set_num_blocks || !set_num_files)
    {
        std::cout << "cannot retrieve number of files or number of blocks from the input" << std::endl;
        exit(1);
    }
    f_stream.close();
    return std::make_pair(num_of_blocks, num_of_files);
}

/**
 * @brief Returns all intersecting blocks between two volumes
 * @param i_volume Path to first volume
 * @param j_volume Path to second volume
 */
std::vector<int> calcVolumeIntersect(std::string i_volume, std::string j_volume) 
{
    std::ifstream i_volume_stream(i_volume.c_str(), std::ifstream::in);
    if (!i_volume_stream.is_open())
    {
        std::cout << "error opening volume file - calcVolumeIntersect" << i_volume << std::endl;
        exit(1);
    }
    
    std::ifstream j_volume_stream(j_volume.c_str(), std::ifstream::in);
    if (!j_volume_stream.is_open())
    {
        std::cout << "error opening volume file - calcVolumeIntersect" << j_volume << std::endl;
        i_volume_stream.close();
        exit(1);
    }
    std::vector<int> blocks_i;
    std::vector<int> blocks_j;
    std::vector<int> intersect_i_j;
    std::string content;
    std::vector<std::string> splitted_content;


    while (std::getline(i_volume_stream, content))
    {
        splitted_content = split_string(content, ",");
        if (splitted_content[0] == "B")
        {
            blocks_i.push_back(std::stoi(splitted_content[1]));
        }
    }
    i_volume_stream.close();

    while (std::getline(j_volume_stream, content))
    {
        splitted_content = split_string(content, ",");
        if (splitted_content[0] == "B")
        {
            blocks_j.push_back(std::stoi(splitted_content[1]));
        }
    }
    j_volume_stream.close();

    std::set_intersection(blocks_i.begin(), blocks_i.end(), blocks_j.begin(), blocks_j.end(), back_inserter(intersect_i_j));
    return intersect_i_j;
}

/**
 * @brief Returns an array correlating each block index with a block size
 * @param source_volume_list A list of all source volume paths
 * @param lastSourceSn_block_file A pair including the last index if the block and of the files
 */
std::vector<double> getBlockSizes(std::vector<std::string> &source_volume_list, std::pair<int, int> &lastSourceSn_block_file) 
{
    std::vector<double> blockSizes(lastSourceSn_block_file.first + 1, 0);
    std::cout << blockSizes.size() << std::endl;
    for (auto &source_volume : source_volume_list) {
        std::ifstream source_volume_stream(source_volume.c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            std::cout << "error opening volume file - getBlockSizes" << source_volume << std::endl;
            exit(1);
        }
        
        std::string content;
        std::vector<std::string> splitted_content;

        while (std::getline(source_volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "F")
            {
                int number_of_blocks_in_file_line = std::stoi(splitted_content[4]);
                for (register int i = 0; i < 2 * number_of_blocks_in_file_line; i += 2) //Read block_sn and block_size simultaneously and add constrains to the model.
                {
                    int block_sn = std::stoi(splitted_content[5 + i]);
                    double size_read = std::stod(splitted_content[6 + i]); //Update block size histogram
                    double size = ((double)size_read) / 1024.0;
                    if (blockSizes[block_sn] == 0)
                    {
                        blockSizes[block_sn] = size;
                    }
                }
            }
        }
        source_volume_stream.close();
	}

    return blockSizes;
}

/**
 * @brief Find the index of the last block and the last file and pile them together
 * @param source_volume_list A list of all source volume paths
 */
std::pair<int, int> getLastBlockAndFileSn(std::vector<std::string> &source_volume_list) 
{
    int lastBlockSn = 0;
    int lastFileSn = 0;
    for (auto &source_volume : source_volume_list) {
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

/**
 * @brief Find which file indexes are used in a volume and store in an array
 * @param source_volume_list A list of all source volume paths
 */
std::vector<std::set<int>> getFileSnInVolumes(std::vector<std::string> &source_volume_list)
{
    std::vector<std::set<int>> fileSnInVolume;
    for (int source = 0; source < source_volume_list.size(); source++) {
        std::ifstream source_volume_stream(source_volume_list[source].c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            std::cout << "error opening volume file - getFileSnInVolumes" << source_volume_list[source] << std::endl;
            exit(1);
        }

        std::set<int> fileSnSet;
        fileSnInVolume.push_back(fileSnSet);
        std::string content;
        std::vector<std::string> splitted_content;
        while (std::getline(source_volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "F")
            {
                fileSnInVolume[source].insert(std::stoi(splitted_content[1]));
            }
        }
        source_volume_stream.close();
	}   
    return fileSnInVolume;
}

/**
 * ---------------------------------------------------------------------------------------------------------------------------------------------
 * Write run statistics
 **/

std::string volumeName;
std::string volumeNumber = "5"; //This was in our expirements the default size
std::string k;
double ingestionTime;
double optimizationTime;
std::vector<double> initialVolume;
std::vector<double> finalVolume;

/**
 * @brief Calculates the run statistics and puts them into global variables accessible to the statistics printer
 * @param name Name of the volume list used
 * @param C_i_s_t The array storing all the Copy type gurobi solution parameters 
 * @param D_i_s The array storing all the Delete type gurobi solution parameters
 * @param source_volume_list List of all source volumes
 * @param block_sizes The size of every block
 * @param intersects_source_target_blocksn A reference to the intersection array
 */
void computeConclusionFields(std::string name, std::vector<std::vector<GRBVar*>> &C_i_s_t, std::vector<GRBVar*> &D_i_s, std::vector<std::string> &source_volume_list, std::vector<double> &block_sizes, std::vector<std::vector<std::vector<int>>> &intersects_source_target_blocksn) {
    auto splitSlashes = split_string(name, "/");
    std::vector<std::string> splitted_content = split_string(replaceAll(splitSlashes[splitSlashes.size() - 1], "_ck8_", "_"), "_"); //ck8 was a common suffix for metadata purposes which was irrelevant to the algorithm so we had it removed in code.
    
    volumeName = splitted_content[0];
    k = replaceAll(splitted_content[1], "k", "");
    if(splitted_content.size() > 2) {
        volumeNumber = splitted_content[2];
    }

    //Compute initial system size
    double initialSystemSize = 0.0;
    for(int source = 0 ; source < intersects_source_target_blocksn.size(); source ++) {
        double size = 0;
        for(int i = 0; i < intersects_source_target_blocksn[source][source].size(); i++) {
            size += block_sizes[intersects_source_target_blocksn[source][source][i]];
        }
        initialSystemSize += size;
        initialVolume.push_back(size);
    }

    //Save the initial fractions of each volume
    for(int volume = 0; volume < initialVolume.size(); volume++) {
        initialVolume[volume] /= initialSystemSize;
    }
    double TotalChange = 0.0;
    double TotalTraffic = 0.0;

    //Create data_structure saving which block is in which volume
    std::vector<std::vector<bool>> exist_s_i;
    for (int source = 0; source < source_volume_list.size(); source++) {
        std::vector<bool> exist_s(C_i_s_t.size(), false);
        std::ifstream source_volume_stream(source_volume_list[source].c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            exit(1);
        }

        std::string content;
        std::vector<std::string> splitted_content;
        while (std::getline(source_volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "B")
            {
                exist_s[std::stoi(splitted_content[1])] = true;
            }
        }
        exist_s_i.push_back(exist_s);
        source_volume_stream.close();
	}

    //Use the gurobi solution saved in the gurobi parameters to track the changes to the system
    //copied blocks included in the intersection do not count for traffic
    for (int source = 0; source < intersects_source_target_blocksn.size(); source++) {
        for (int target = 0; target < intersects_source_target_blocksn[0].size(); target++) {
            std::vector<int> relevantIntersect = intersects_source_target_blocksn[source][target];
            if(source == target) {
                continue;
            }
            int previousIntersectIndex = -1;
            for(int j = 0; j < relevantIntersect.size(); j++) {
                int nextIntersectIndex = relevantIntersect[j];
                for(int i = previousIntersectIndex + 1; i < nextIntersectIndex; i++) {
                    if (C_i_s_t[i][source][target].get(GRB_DoubleAttr_X) != 0.0)
                    {
                        TotalChange += block_sizes[i];
                        TotalTraffic += block_sizes[i];
                    }
                }
                previousIntersectIndex = nextIntersectIndex;
            }
            for(int i = previousIntersectIndex + 1; i < C_i_s_t.size(); i++) {
                if (C_i_s_t[i][source][target].get(GRB_DoubleAttr_X) != 0.0)
                {
                    TotalChange += block_sizes[i];
                    TotalTraffic += block_sizes[i];
                }
            }
        }

        for(int i = 0; i < C_i_s_t.size(); i++) {
            if (D_i_s[i][source].get(GRB_DoubleAttr_X) != 0.0)
                {
                    TotalChange -= block_sizes[i];
                }
        }
    }

    //Find the new system size by finding how many blocks weren't deleted and how many new blocks were created
    double size = 0.0;
    double fraction = 0.0;
    for (int target = 0; target < intersects_source_target_blocksn[0].size(); target++) {
        double SumNewBlocks = 0.0;
        double SumOldBlocks = 0.0;
        for(int i = 0; i < block_sizes.size(); i++) {
            if(exist_s_i[target][i]) {
                if (D_i_s[i][target].get(GRB_DoubleAttr_X) == 0.0)
                {
                    SumOldBlocks += block_sizes[i];
                }
            }
        }
        for (int source = 0; source < intersects_source_target_blocksn.size(); source++) {
            std::vector<int> relevantIntersect = intersects_source_target_blocksn[source][target];
            if(source == target) {
                continue;
            }
            int previousIntersectIndex = -1;
            for(int j = 0; j < relevantIntersect.size(); j++) {
                int nextIntersectIndex = relevantIntersect[j];
                for(int i = previousIntersectIndex + 1; i < nextIntersectIndex; i++) {
                    if (C_i_s_t[i][source][target].get(GRB_DoubleAttr_X) != 0.0)
                    {
                        SumNewBlocks += block_sizes[i];
                    }
                }
                previousIntersectIndex = nextIntersectIndex;
            }
            for(int i = previousIntersectIndex + 1; i < C_i_s_t.size(); i++) {
                if (C_i_s_t[i][source][target].get(GRB_DoubleAttr_X) != 0.0)
                {
                    SumNewBlocks += block_sizes[i];
                }
            }
        }

        finalVolume.push_back((SumNewBlocks + SumOldBlocks) / (initialSystemSize + TotalChange));
    }   

    actual_deletion_Percentage = TotalChange/initialSystemSize;
    actual_traffic_Percentage = TotalTraffic/initialSystemSize;

}

/**
 * @brief Prints the conclusion line into the statistics file
 */
void save_runConclusionLine()
{
    std::ofstream out(conclusionFileName, std::ios_base::app);
    if (!out)
    {
        std::cout << "Cannot open output file\n";
    }   
    out << "ILP" << ", "
        << volumeName << ", "
        << volumeNumber << ", "
        << k << ", "
        << T_percentage << ", "
        << Margin << ", "
        << seed << ", "
        << model_time_limit << ","
        << ingestionTime << ", "
        << optimizationTime << ","
        << actual_traffic_Percentage << ", "
        << actual_deletion_Percentage << ", " ;
        
    for(int i = 0; i < initialVolume.size(); i++) {
        out << initialVolume[i] << ", "
            << finalVolume[i] << ", " ;
    }
    out << std::endl;
    out.close();
}

/**
 * @brief Prints the failure reason into where the solution should be written
 * 
 * @param print_to Where the solution should be written 
 * @param reason String of error
 */
void recordFail(std::string print_to, std::string reason) {
    std::ofstream solution(print_to, std::ios_base::app);
    if (!solution)
    {
        std::cout << "Cannot open output file" << print_to << std::endl;
        exit(1);
    }
    solution << "failed"; 
    solution << reason;
}

/**
 * @brief Records the input of the run, for debugging purposes
 * 
 * @param write_solution Where the solution should be written
 * @param volume_list The list of all volume files and their desired balance
 * @param T_percentage Input max traffic
 * @param Margin Input error margin
 * @param InnerFilterSize Input inner filter
 * @param model_time_limit Input time limit for optimization
 * @param seed The random seed for the gurobi solver
 */
void saveInput(std::string write_solution, std::string volume_list, double T_percentage, double Margin, int InnerFilterSize, int model_time_limit, std::string seed) {
    std::ofstream solution(write_solution, std::ios_base::app);
    if (!solution)
    {
        std::cout << "Cannot open output file" << write_solution << std::endl;
        exit(1);
    }
    solution << "volume_list: " << volume_list << std::endl;
    solution << "T_percentage: " << T_percentage<< std::endl;
    solution << "Margin: " << Margin << std::endl;
    solution << "InnerFilterSize: " << InnerFilterSize << std::endl;
    solution << "model_time_limit: " << model_time_limit << std::endl;
    solution << "seed: " << seed << std::endl;
}

/**
 * @brief Records information on solver performance
 * 
 * @param write_soprint_tolution Where the solution should be written
 * @param modelCreationTime How long it took to create the model itseld (s)
 * @param solverTime How long the gurobi solver took (s)
 * @param solution_status What status has the gurobi solver given the run
 */
void saveRunSolvingData(std::string print_to, int modelCreationTime, int solverTime, std::string solution_status) {
    std::ofstream solution(print_to, std::ios_base::app);
    if (!solution)
    {
        std::cout << "Cannot open output file" << print_to << std::endl;
        exit(1);
    }
    solution << "solver time(s): " << std::to_string(solverTime) << std::endl;
    solution << "model creation time(s): " << std::to_string(modelCreationTime) << std::endl;
    solution << "status: " << solution_status << std::endl;
}

/**
 * @brief Saves the gurobi solution as a migration plan (which file goes from where to where)
 * 
 * @param X_l_s_t The array saving the soltion file transfer gurobi parameter type
 * @param C_i_s_t The array saving the solution copy block gurobi parameter type
 * @param targetVolumeNum The number of target volumes in the system
 * @param print_to Where the solution should be written
 */
void saveSolution(std::vector<std::vector<GRBVar*>> &X_l_s_t, std::vector<std::vector<GRBVar*>> &C_i_s_t, int targetVolumeNum, std::string print_to) {
    std::ofstream solution(print_to, std::ios_base::app);
    if (!solution)
    {
        std::cout << "Cannot open output file" << print_to << std::endl;
        exit(1);
    }

    for (int i = 0; i < X_l_s_t.size(); i++) {
        for (int s = 0; s < X_l_s_t[i].size(); s++) {
            for (int t = 0; t < targetVolumeNum; t++) {
                if (X_l_s_t[i][s][t].get(GRB_DoubleAttr_X) != 0.0) //file is moved
                {
                    solution << i << ": " << s << " -> " << t <<  std::endl;
                }
            }
        }
    } 

    solution << "X_l_s_t size" <<X_l_s_t.size() << std::endl;
    if(X_l_s_t.size()) {
        solution <<"," << X_l_s_t[0].size() << "," << targetVolumeNum << std::endl;    
    }
}

/**
 * @brief Saves the assumed system balance (assumed since we work on a filtered portion of the block amount)
 * 
 * @param C_i_s_t The array saving the solution copy block gurobi parameter type
 * @param D_i_s The array saving the soltion delete block gurobi parameter type
 * @param source_volume_list List of all source volumes
 * @param block_sizes Array correlatin block index to its size
 * @param intersects_source_target_blocksn Array of blocks in intersect for each volume
 */
void saveFilterBalance(std::string print_to, std::vector<std::vector<GRBVar*>> &C_i_s_t, std::vector<GRBVar*> &D_i_s, std::vector<std::string> &source_volume_list, std::vector<double> &block_sizes, std::vector<std::vector<std::vector<int>>> &intersects_source_target_blocksn) {
    std::ofstream solution(print_to, std::ios_base::app);
    if (!solution)
    {
        std::cout << "Cannot open output file" << print_to << std::endl;
        exit(1);
    }
    double TotalChange = 0.0;
    double TotalDelete = 0.0;
    double TotalTraffic = 0.0;
    std::vector<double> totalDeletions(intersects_source_target_blocksn.size(), 0);
    std::vector<std::vector<double> > totalTraffics( intersects_source_target_blocksn.size(), std::vector<double>(intersects_source_target_blocksn[0].size(), 0));

    //Create data structure for which volumes contain which block
    std::vector<std::vector<bool>> exist_s_i;
    for (int source = 0; source < source_volume_list.size(); source++) {
        std::vector<bool> exist_s(C_i_s_t.size(), false);
        std::ifstream source_volume_stream(source_volume_list[source].c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            solution << "error opening volume file - saveFilterBalance" << source_volume_list[source] << std::endl;
            exit(1);
        }

        std::string content;
        std::vector<std::string> splitted_content;
        while (std::getline(source_volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "B")
            {
                exist_s[std::stoi(splitted_content[1])] = true;
            }
        }
        exist_s_i.push_back(exist_s);
        source_volume_stream.close();
	}

    //Use the solution by gurobi saved in the gurobi parameters to track of volume size changes
    //copied blocks inside the intersection do not count for traffic
    for (int source = 0; source < intersects_source_target_blocksn.size(); source++) {
        for (int target = 0; target < intersects_source_target_blocksn[0].size(); target++) {
            std::vector<int> relevantIntersect = intersects_source_target_blocksn[source][target];
            if(source == target) {
                continue;
            }
            for(int i = 0; i < C_i_s_t.size(); i++) {
                if(exist_s_i[target][i]) {
                    continue;
                } 
                if (C_i_s_t[i][source][target].get(GRB_DoubleAttr_X) != 0.0) {
                    TotalTraffic += block_sizes[i];
                    totalTraffics[source][target] += block_sizes[i];
                }
            }             
        }

        for(int i = 0; i < C_i_s_t.size(); i++) {
            if (D_i_s[i][source].get(GRB_DoubleAttr_X) != 0.0)
                {
                    TotalDelete -= block_sizes[i];
                    totalDeletions[source] -= block_sizes[i];            
                }
        }
    }

    //Track the current size by computing how many blocks were not deleted and how many new blocks were created
    double size = 0.0;
    double fraction = 0.0;
    for (int target = 0; target < intersects_source_target_blocksn[0].size(); target++) {
        double SumNewBlocks = 0.0;
        double SumOldBlocks = 0.0;
        for(int i = 0; i < block_sizes.size(); i++) {
            if(exist_s_i[target][i]) {
                if (D_i_s[i][target].get(GRB_DoubleAttr_X) == 0.0)
                {
                    SumOldBlocks += block_sizes[i];
                }
            }
        }
        for (int source = 0; source < intersects_source_target_blocksn.size(); source++) {
            std::vector<int> relevantIntersect = intersects_source_target_blocksn[source][target];
            if(source == target) {
                continue;
            }
            for(int i = 0; i < C_i_s_t.size(); i++) {
                if(exist_s_i[target][i]) {
                    continue;
                } 
                if (C_i_s_t[i][source][target].get(GRB_DoubleAttr_X) != 0.0) {
                    SumNewBlocks += block_sizes[i];
                }
            }            
        }
        TotalChange = TotalDelete + TotalTraffic;
        solution << std::endl << "volume " << target << ":" << SumNewBlocks + SumOldBlocks << std::endl;
        solution << "volume " << target << ":" << (SumNewBlocks + SumOldBlocks) / (SOURCE_SIZE + TotalChange) << std::endl;
        size += SumNewBlocks + SumOldBlocks;
        fraction += (SumNewBlocks + SumOldBlocks) / (SOURCE_SIZE + TotalChange);
    }  

        solution << std::endl << "original size: " << SOURCE_SIZE << std::endl;
        solution << "total delete " << TotalDelete << std::endl;
        solution << "total traffic " << TotalTraffic << std::endl;
        for (int source = 0; source < intersects_source_target_blocksn.size(); source++) {
            solution << "volume: " << source << " deletion: " << totalDeletions[source] << std::endl;           
            double traffic = 0.0;
            for (int target = 0; target < intersects_source_target_blocksn[0].size(); target++) {    
                traffic +=  totalTraffics[target][source];
            }
            solution << "volume: " << source << " traffic: " << traffic << std::endl;           
        }
        solution << "total change " << TotalChange << std::endl;        
        solution << "new size by change: " << SOURCE_SIZE + TotalChange << std::endl;
        solution << "new size by block counting: " << size << std::endl;
        solution << fraction << std::endl;
}

void saveMigrationPlan(std::vector<std::string> &source_volume_list, std::vector<std::set<int>> &fileSnInVolume, std::vector<std::vector<GRBVar*>> &X_l_s_t, std::string printTo, int targetVolumeNum) {
    std::vector<std::vector<std::string> > migrationPlanCSV;
    std::ofstream migrationPlan(printTo, std::ios_base::trunc);
    std::vector<std::set<int>> fileSnInVolumeFinalState(fileSnInVolume);

    for (int i = 0; i < X_l_s_t.size(); i++) {
        for (int s = 0; s < X_l_s_t[i].size(); s++) {
            for (int t = 0; t < targetVolumeNum; t++) {
                if (X_l_s_t[i][s][t].get(GRB_DoubleAttr_X) != 0.0) //file is moved
                {
                    fileSnInVolumeFinalState[s].erase(i);
                    fileSnInVolumeFinalState[t].insert(i);
                }
            }
        }
    } 

    if (!migrationPlan)
    {
        std::cout << "Cannot open output file" << printTo << std::endl;
        exit(1);
    }

    for(int i = 0; i < source_volume_list.size(); i++) {
        migrationPlanCSV.push_back(std::vector<std::string>());
        migrationPlanCSV[i].push_back(source_volume_list[i]);
        std::string fileList = "";
        for(auto elem : fileSnInVolume[i]) {
            fileList += std::to_string(elem);
            fileList += "-";
        }
        if(fileList != "") {
            fileList.pop_back();
        }
        migrationPlanCSV[i].push_back(fileList);
        
        fileList = "";
        for(auto elem : fileSnInVolumeFinalState[i]) {
            fileList += std::to_string(elem);
            fileList += "-";
        }
        if(fileList != "") {
            fileList.pop_back();
        }
        migrationPlanCSV[i].push_back(fileList);        
        migrationPlanCSV[i].push_back("1");        
    }

    for(int i = 0; i < migrationPlanCSV.size(); i++) {
        for(int j = 0; j < migrationPlanCSV[i].size(); j++) {
            migrationPlan << migrationPlanCSV[i][j] << ",";
        }
        migrationPlan << std::endl;
    }    


}


/**
 * ---------------------------------------------------------------------------------------------------------------------------------------------
 * Add gurobi constraints.
 * The contraint enumaration is similar to the article, with constraints 0a, 0b not mention in the article due to them being technical constraints.
 **/

// Contraint 0a. Blocks inside the intersect are considered to be copied
void addConstraint_allIntersectsAreCopied(GRBModel &model, std::vector<std::vector<std::vector<int>>> &intersects_source_target_blocksn, std::vector<std::vector<GRBVar*>> &C_i_s_t)
{
    int constraintAdded = 0;
    for (int source = 0; source < intersects_source_target_blocksn.size(); source++) {
        for (int target = 0; target < intersects_source_target_blocksn[source].size(); target++) {
            for (int i = 0; i < intersects_source_target_blocksn[source][target].size(); i++) {
                model.addConstr(C_i_s_t[intersects_source_target_blocksn[source][target][i]][source][target], GRB_EQUAL, 1) ;
            }      
        }
	}
}

// Contraint 0b. A file cannot be remapped from somewhere it doesn't exist
void addConstraint_DontRemapNonExistantFilesOfRemapToSelf(GRBModel &model, std::vector<std::vector<GRBVar*>> &X_l_s_t, std::vector<std::set<int>> &fileSnInVolumes, int numOfTargetVolumes)
{
    for (int l = 0; l < X_l_s_t.size(); l++) {
        for (int source = 0; source < X_l_s_t[l].size(); source ++) {
            model.addConstr(X_l_s_t[l][source][source], GRB_EQUAL, 0);
            if(fileSnInVolumes[source].count(l) == 1) {
                continue;
            } else {
                for (int target = 0; target < numOfTargetVolumes; target++) {
                    model.addConstr(X_l_s_t[l][source][target], GRB_EQUAL, 0);
                }   
            }
        }
    }
}

// Contraint 2. A file can be remapped to at most one volume.
void addConstraint_RemapFilesToOnlyOneVolume(GRBModel &model, std::vector<std::vector<GRBVar*>> &X_l_s_t, std::vector<std::pair<int, int>> &num_of_blocks_and_files_sourceVolumes, int numOfTargetVolumes)
{
    for (int l = 0; l < X_l_s_t.size(); l++) {
        for (int source = 0; source < X_l_s_t[l].size(); source ++) {
            GRBLinExpr Sum_X_l_s_t = 0.0;
            for (int target = 0; target < numOfTargetVolumes; target++) {
                Sum_X_l_s_t += X_l_s_t[l][source][target];
            }
            model.addConstr(Sum_X_l_s_t <= 1);  
        }
	}
}

// Contraint 3. A block can only be deleted or copied from a volume it was originally stored in.
void addConstraint_blockNeedsToExistToBeCloneOrDeleated(GRBModel &model, std::vector<std::vector<GRBVar*>> &C_i_s_t, std::vector<GRBVar*> &D_i_s, std::vector<std::string> &source_volume_list, std::vector<std::string> &target_volume_list)
{
    std::vector<std::vector<bool>> exist_s_i;
    for (int source = 0; source < source_volume_list.size(); source++) {
        std::vector<bool> exist_s(C_i_s_t.size(), false);
        std::ifstream source_volume_stream(source_volume_list[source].c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            std::cout << "error opening volume file - addConstraint_blockNeedsToExistToBeCloneOrDeleated" << source_volume_list[source] << std::endl;
            exit(1);
        }

        std::string content;
        std::vector<std::string> splitted_content;
        while (std::getline(source_volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "B")
            {
                exist_s[std::stoi(splitted_content[1])] = true;
            }
        }
        exist_s_i.push_back(exist_s);
        source_volume_stream.close();
	}    

    for (int block = 0; block < C_i_s_t.size(); block++) {
        for (int source = 0; source < source_volume_list.size(); source ++) {
            if(!exist_s_i[source][block]) {
                model.addConstr(D_i_s[block][source], GRB_EQUAL, 0);
                for (int target = 0; target < target_volume_list.size(); target ++) {
                    model.addConstr(C_i_s_t[block][source][target], GRB_EQUAL, 0);
                }
            }
        }
    }
}

// Contraint 4. A block can be deleted from a volume only if all the files containing it are remapped to other volumes.
// Contraint 5. A block can be deleted from a volume only if no file containing it is remapped to this volume.
// Contraint 6. When a file is remapped, all its blocks are either copied to the target volume, or are initially there.
void add_Three_constraints_BlockIsDeletedOnlyIfNoLocalFileUsingItRemains_BlockIsDeletedOnlyIfNoNewFileTransferredIsUsingIt_RemmapedFileHasAllItsBlocks(GRBModel &model, std::vector<std::vector<GRBVar*>> X_l_s_t, std::vector<std::vector<GRBVar*>> &C_i_s_t, std::vector<GRBVar*> &D_i_s, std::vector<std::string> &source_volume_list, int numOfTargetVolumes) {
    std::vector<std::vector<bool>> exist_s_i;

    for (int source = 0; source < source_volume_list.size(); source++) {
        std::ifstream source_volume_stream(source_volume_list[source].c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            std::cout << "error opening volume file - add_Three_constraints_BlockIsDeletedOnlyIfNoLocalFileUsingItRemains_BlockIsDeletedOnlyIfNoNewFileTransferredIsUsingIt_RemmapedFileHasAllItsBlocks" << source_volume_list[source] << std::endl;
            exit(1);
        }


        std::string content;
        std::vector<std::string> splitted_content;
        while (std::getline(source_volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "F")
            {
                int fileSn = std::stoi(splitted_content[1]);
                int number_of_blocks_in_file_line = std::stoi(splitted_content[4]);
                for (register int i = 0; i < 2 * number_of_blocks_in_file_line; i += 2) //read block_sn and block_size simultaneously and add constrains to the model.
                {
                    int block_sn = std::stoi(splitted_content[5 + i]);
                    GRBLinExpr BlockIsDeletedOnlyIfNoLocalFileUsingItRemains = 0.0;
                    for (int target = 0; target < numOfTargetVolumes; target++) {
                        GRBLinExpr RemmapedFileHasAllItsBlocks = 0.0;
                        model.addConstr(D_i_s[block_sn][target] <= 1 -  X_l_s_t[fileSn][source][target]); // BlockIsDeletedOnlyIfNoNewFileTransferredIsUsingIt
                        BlockIsDeletedOnlyIfNoLocalFileUsingItRemains += X_l_s_t[fileSn][source][target];
                        for (int arbitraryVolume = 0; arbitraryVolume < source_volume_list.size(); arbitraryVolume++) {
                            if(arbitraryVolume == target) {
                                continue;
                            }
                            RemmapedFileHasAllItsBlocks += C_i_s_t[block_sn][arbitraryVolume][target];
                        }
                        model.addConstr(X_l_s_t[fileSn][source][target] <= RemmapedFileHasAllItsBlocks);
                    }
                    model.addConstr(D_i_s[block_sn][source] <= BlockIsDeletedOnlyIfNoLocalFileUsingItRemains);
                }
            }
        }
        source_volume_stream.close();
	}    
}

// Contraint 7. A block can be copied to a target volume only from one source volume.
void addConstraint_doNotCopyFromMoreThanOneSource(GRBModel &model, std::vector<std::vector<GRBVar*>> &C_i_s_t, std::vector<std::string> &source_volume_list, int numOfTargetVolumes)
{
    std::vector<std::vector<bool>> exist_s_i;
    for (int source = 0; source < source_volume_list.size(); source++) {
        std::vector<bool> exist_s(C_i_s_t.size(), false);
        std::ifstream source_volume_stream(source_volume_list[source].c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            std::cout << "error opening volume file - addConstraint_doNotCopyFromMoreThanOneSource" << source_volume_list[source] << std::endl;
            exit(1);
        }

        std::string content;
        std::vector<std::string> splitted_content;
        while (std::getline(source_volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "B")
            {
                exist_s[std::stoi(splitted_content[1])] = true;
            }
        }
        exist_s_i.push_back(exist_s);
        source_volume_stream.close();
	}    


    for (int i = 0; i < C_i_s_t.size(); i++) {
        for(int target = 0; target < numOfTargetVolumes; target++) {
            GRBLinExpr sumS = 0.0;
            for(int source = 0; source < C_i_s_t[0].size(); source ++) {
                if(exist_s_i[target][i] ) {
                    continue;
                }
                sumS += C_i_s_t[i][source][target];
            }
            model.addConstr(sumS <= 1);
        }
    }
}

// Contraint 8. A block must be deleted if there are no files containing it on the volume.
void addConstraint_no_orphans(GRBModel &model, std::vector<std::vector<GRBVar*>> &X_l_s_t, std::vector<std::vector<GRBVar*>> &C_i_s_t, std::vector<GRBVar*> &D_i_s, std::vector<std::string> &source_volume_list, int numOfTargetVolumes, int numOfBlocks, int numOfFiles) {
    std::vector<std::vector<bool>> block_i_in_File_f(numOfBlocks + 1, std::vector<bool>(numOfFiles + 1, false));
    std::vector<std::vector<bool>> volume_v_has_file_f(source_volume_list.size(), std::vector<bool>(numOfFiles + 1, false));
    std::vector<std::vector<bool>> volume_v_has_block_i(source_volume_list.size(), std::vector<bool>(numOfBlocks + 1, false));;
    for (int source = 0; source < source_volume_list.size(); source++) {
        std::ifstream source_volume_stream(source_volume_list[source].c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            std::cout << "error opening volume file - addConstraint_no_orphans" << source_volume_list[source] << std::endl;
            exit(1);
        }

        std::string content;
        std::vector<std::string> splitted_content;

        while (std::getline(source_volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "F")
            {
                int fileSn = std::stoi(splitted_content[1]);
                volume_v_has_file_f[source][fileSn] = true;
                int number_of_blocks_in_file_line = std::stoi(splitted_content[4]);
                for (register int i = 0; i < 2 * number_of_blocks_in_file_line; i += 2) //read block_sn and block_size simultaneously and add constrains to the model.
                {
                    int block_sn = std::stoi(splitted_content[5 + i]);
                    
                    block_i_in_File_f[block_sn][fileSn] = true;
                    volume_v_has_block_i[source][block_sn] = true;
                }
            }
        }
        source_volume_stream.close();
	}    

    for (int block = 0; block <= numOfBlocks; block++) {
        for (int source = 0; source < source_volume_list.size(); source++) {
            if(!volume_v_has_block_i[source][block]){
                continue;
            }
            GRBLinExpr Sum_l_s_1MinusSum_t_Xl_s_st = 0.0;
            GRBLinExpr Sum_l_tXl_t_ts = 0.0;

            for(int target = 0; target < numOfTargetVolumes; target++) {
                GRBLinExpr Sum_l_t_Xlt_t_s = 0.0;
                for(int file = 0; file <= numOfFiles; file++) {
                    if(block_i_in_File_f[block][file]) {
                        if(volume_v_has_file_f[target][file]) {
                            Sum_l_tXl_t_ts += X_l_s_t[file][target][source];
                        }
                    }
                }
            }

            for(int file = 0; file <= numOfFiles; file++) {
                if(block_i_in_File_f[block][file]) {
                    if(volume_v_has_file_f[source][file]) {
                        GRBLinExpr Sum_l_sXl_s_st = 0.0;   
                        for(int target = 0; target < numOfTargetVolumes; target++) {
                            Sum_l_sXl_s_st += X_l_s_t[file][source][target];
                        }
                        Sum_l_s_1MinusSum_t_Xl_s_st += 1 - Sum_l_sXl_s_st;
                    }
                }
            }
            model.addConstr(D_i_s[block][source] >= 1 - (Sum_l_s_1MinusSum_t_Xl_s_st + Sum_l_tXl_t_ts));
        }
    }

}

// Contraint 9. A block cannot be copied to a target volume if no file will contain it there.
void addConstraint_no_extraBlocks(GRBModel &model, std::vector<std::vector<GRBVar*>> &X_l_s_t, std::vector<std::vector<GRBVar*>> &C_i_s_t, std::vector<GRBVar*> &D_i_s, std::vector<std::string> &source_volume_list, int numOfTargetVolumes, int numOfBlocks, int numOfFiles) {
    std::vector<std::vector<bool>> block_i_in_File_f(numOfBlocks + 1, std::vector<bool>(numOfFiles + 1, false));
    std::vector<std::vector<bool>> volume_v_has_file_f(source_volume_list.size(), std::vector<bool>(numOfFiles + 1, false));
    std::vector<std::vector<bool>> volume_v_has_block_i(source_volume_list.size(), std::vector<bool>(numOfBlocks + 1, false));;
    for (int source = 0; source < source_volume_list.size(); source++) {
        std::ifstream source_volume_stream(source_volume_list[source].c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            std::cout << "error opening volume file - addConstraint_no_extraBlocks" << source_volume_list[source] << std::endl;
            exit(1);
        }

        std::string content;
        std::vector<std::string> splitted_content;

        while (std::getline(source_volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "F")
            {
                int fileSn = std::stoi(splitted_content[1]);
                volume_v_has_file_f[source][fileSn] = true;
                int number_of_blocks_in_file_line = std::stoi(splitted_content[4]);
                for (register int i = 0; i < 2 * number_of_blocks_in_file_line; i += 2) //Read block_sn and block_size simultaneously and add constrains to the model.
                {
                    int block_sn = std::stoi(splitted_content[5 + i]);
                    
                    block_i_in_File_f[block_sn][fileSn] = true;
                    volume_v_has_block_i[source][block_sn] = true;
                }
            }
        }
        source_volume_stream.close();
	}     

    for(int target = 0; target < numOfTargetVolumes; target++) {
        for(int block = 0; block < C_i_s_t.size(); block++) {
            if(volume_v_has_block_i[target][block]) {
                continue;
            }
            GRBLinExpr sumCist = 0.0;
            GRBLinExpr sumXlst = 0.0;

            for(int source = 0; source < C_i_s_t[0].size(); source++) {
                if(!volume_v_has_block_i[source][block]) {
                    continue;
                }
                sumCist += C_i_s_t[block][source][target];
                for(int file = 0; file < X_l_s_t.size(); file++) {
                    if(!volume_v_has_file_f[source][file]) {
                        continue;
                    }
                    if(!block_i_in_File_f[block][file]) {
                        continue;
                    }                    
                    sumXlst += X_l_s_t[file][source][target];
                }
            }
            model.addConstr(sumCist <= sumXlst);
        }
    }
}

// Contraint 10. The size of all the copied blocks is not higher than themaximum allowed traffic.
// Set objective: Maximize the sum of sizes of all blocks that are deleted minus all the blocks that are copied. This is equivalent to minimizing the overall system size.
void addConstraint_traffic_and_setObjective(GRBModel &model, std::vector<std::vector<GRBVar*>> &C_i_s_t, std::vector<GRBVar*> &D_i_s, std::vector<std::vector<std::vector<int>>> &intersects_source_target_blocksn, std::vector<double> &block_sizes, int maximumTrafficPercentage) {
    double maxTrafficInKbytes = SOURCE_SIZE * maximumTrafficPercentage / 100;                //assign the number of bytes to migrate

    GRBLinExpr trafficConstraint = 0.0;
    GRBLinExpr objective = 0.0;

    for (int source = 0; source < intersects_source_target_blocksn.size(); source++) {
        bool setSourceDel = false;
        for (int target = 0; target < intersects_source_target_blocksn[0].size(); target++) {
            std::vector<int> relevantIntersect = intersects_source_target_blocksn[source][target];
            if(source == target) {
                continue;
            }
            int previousIntersectIndex = -1;
            for(int j = 0; j < relevantIntersect.size(); j++) {
                int nextIntersectIndex = relevantIntersect[j];
                for(int i = previousIntersectIndex + 1; i < nextIntersectIndex; i++) {
                    trafficConstraint += C_i_s_t[i][source][target] * block_sizes[i];
                    objective += C_i_s_t[i][source][target] * block_sizes[i];
                    if(!setSourceDel) {
                        objective -= D_i_s[i][source] * block_sizes[i];
                    }
                }
                if(!setSourceDel) {
                    objective -= D_i_s[nextIntersectIndex][source] * block_sizes[nextIntersectIndex];
                }
                previousIntersectIndex = nextIntersectIndex;
            }
            for(int i = previousIntersectIndex + 1; i < C_i_s_t.size(); i++) {
                trafficConstraint += C_i_s_t[i][source][target] * block_sizes[i];
                objective += C_i_s_t[i][source][target] * block_sizes[i];
                if(!setSourceDel) {
                    objective -= D_i_s[i][source] * block_sizes[i];
                }
            }
            setSourceDel = true;
        }
    }
    model.addConstr(trafficConstraint <= maxTrafficInKbytes);
    model.setObjective(objective, GRB_MINIMIZE); 
}

void addConstraint_traffic_and_setObjectiveFullMigration(GRBModel &model, std::vector<std::vector<GRBVar*>> &C_i_s_t, std::vector<GRBVar*> &D_i_s, std::vector<std::vector<std::vector<int>>> &intersects_source_target_blocksn, std::vector<double> &block_sizes, int maximumTrafficPercentage) {
    
    double SOURCE_SIZE = 0.0;
    for(int source = 0 ; source < intersects_source_target_blocksn.size(); source ++) {
        for(int i = 0; i < intersects_source_target_blocksn[source][source].size(); i++) {
            SOURCE_SIZE += block_sizes[intersects_source_target_blocksn[source][source][i]];
        }
    }
    double maxTrafficInKbytes = SOURCE_SIZE * maximumTrafficPercentage / 100;                //Assign the number of bytes to migrate

    GRBLinExpr trafficConstraint = 0.0;
    GRBLinExpr objective = 0.0;

    for (int source = 0; source < intersects_source_target_blocksn.size(); source++) {
        bool setSourceDel = false;
        for (int target = 0; target < intersects_source_target_blocksn[0].size(); target++) {
            std::vector<int> relevantIntersect = intersects_source_target_blocksn[source][target];
            if(source == target) {
                continue;
            }
            int previousIntersectIndex = -1;
            for(int j = 0; j < relevantIntersect.size(); j++) {
                int nextIntersectIndex = relevantIntersect[j];
                for(int i = previousIntersectIndex + 1; i < nextIntersectIndex; i++) {
                    trafficConstraint += C_i_s_t[i][source][target] * block_sizes[i];
                    objective += C_i_s_t[i][source][target] * block_sizes[i];
                    if(!setSourceDel) {
                        objective -= D_i_s[i][source] * block_sizes[i];
                    }
                }
                if(!setSourceDel) {
                    objective -= D_i_s[nextIntersectIndex][source] * block_sizes[nextIntersectIndex];
                }
                previousIntersectIndex = nextIntersectIndex;
            }
            for(int i = previousIntersectIndex + 1; i < C_i_s_t.size(); i++) {
                trafficConstraint += C_i_s_t[i][source][target] * block_sizes[i];
                objective += C_i_s_t[i][source][target] * block_sizes[i];
                if(!setSourceDel) {
                    objective -= D_i_s[i][source] * block_sizes[i];
                }
            }
            setSourceDel = true;
        }
    }
    model.addConstr(trafficConstraint <= maxTrafficInKbytes);
    model.setObjective(objective, GRB_MINIMIZE); 
}

// Contraint 11. For each volume v,
// (wv−μ)×Size(S′) ≤ size(v′) ≤ (wv+μ)×Size(S′),
// where size(v′) is the volume size after migration,
// i.e. The sum of its non-deleted blocks and blocks copied to it.
void addConstraint_loadBalancing(GRBModel & model, std::vector<std::vector<GRBVar*>>& C_i_s_t,  std::vector<GRBVar*> &D_i_s, std::vector<std::string> &source_volume_list, std::vector<double> &load_volume_list, std::vector<std::vector<std::vector<int>>> &intersects_source_target_blocksn, std::vector<double> &block_sizes, double Margin, int innerFilterSize) {

    GRBLinExpr TotalChange = 0.0;
    std::vector<std::vector<bool>> exist_s_i;
    std::vector<bool> includeInInnerFilter(C_i_s_t.size(), true);
    for (int source = 0; source < source_volume_list.size(); source++) {
        std::vector<bool> exist_s(C_i_s_t.size(), false);
        std::ifstream source_volume_stream(source_volume_list[source].c_str(), std::ifstream::in);
        if (!source_volume_stream.is_open())
        {
            std::cout << "error opening volume file - addConstraint_loadBalancing" << source_volume_list[source] << std::endl;
            exit(1);
        }

        std::string content;
        std::vector<std::string> splitted_content;
        while (std::getline(source_volume_stream, content))
        {
            splitted_content = split_string(content, ",");
            if (splitted_content[0] == "B")
            {
                exist_s[std::stoi(splitted_content[1])] = true;
                if(InnerFilterSize) {
                    std::string threeLastChars = splitted_content[2].substr( splitted_content[2].length() - 3 );
                    std::string threeLastCharsBits = std::string(hex_char_to_bin(threeLastChars[0])) + std::string(hex_char_to_bin(threeLastChars[1])) + std::string(hex_char_to_bin(threeLastChars[2]));;
                    std::string filterBits = threeLastCharsBits.substr( threeLastCharsBits.length() - InnerFilterSize );
                    if(filterBits.find("1") != std::string::npos) {
                        std::cout << splitted_content[2] << " " << splitted_content[2].substr( splitted_content[2].length() - 3 ) << " " << filterBits << std::endl;
                        includeInInnerFilter[std::stoi(splitted_content[1])] = false;
                    }
                }
            }
        }
        exist_s_i.push_back(exist_s);
        source_volume_stream.close();
	}

    std::cout << "total: " << C_i_s_t.size() << " in filter: " << count(includeInInnerFilter.begin(), includeInInnerFilter.end(), true) << std::endl;

    for (int source = 0; source < intersects_source_target_blocksn.size(); source++) {
        for (int target = 0; target < intersects_source_target_blocksn[0].size(); target++) {
            std::vector<int> relevantIntersect = intersects_source_target_blocksn[source][target];
            if(source == target) {
                continue;
            }
            for(int i = 0; i < C_i_s_t.size(); i++) {
                if(exist_s_i[target][i]) {
                    continue;
                } 
                TotalChange += C_i_s_t[i][source][target] * block_sizes[i];
            }
        }

        for(int i = 0; i < C_i_s_t.size(); i++) {
            if(includeInInnerFilter[i]) {
                TotalChange -= D_i_s[i][source] * block_sizes[i];
            }
        }
    }

    for (int target = 0; target < intersects_source_target_blocksn[0].size(); target++) {
        GRBLinExpr SumNewBlocks = 0.0;
        GRBLinExpr SumOldBlocks = 0.0;
        for(int i = 0; i < block_sizes.size(); i++) {
            if(includeInInnerFilter[i]) {
                if(exist_s_i[target][i]) {
                    SumOldBlocks += (1 - D_i_s[i][target]) * block_sizes[i];
                }
            }
        }
        for (int source = 0; source < intersects_source_target_blocksn.size(); source++) {
            std::vector<int> relevantIntersect = intersects_source_target_blocksn[source][target];
            if(source == target) {
                continue;
            }
            for(int i = 0; i < C_i_s_t.size(); i++) {
                if(exist_s_i[target][i]) {
                    continue;
                } 
                SumNewBlocks += C_i_s_t[i][source][target] * block_sizes[i];
            }            
        }
        model.addConstr(SumNewBlocks + SumOldBlocks <= (SOURCE_SIZE + TotalChange) * (load_volume_list[target] + Margin));
        model.addConstr(SumNewBlocks + SumOldBlocks >= (SOURCE_SIZE + TotalChange) * (load_volume_list[target] - Margin));
    }   
}

/**
 * ---------------------------------------------------------------------------------------------------------------------------------------------
 * Core logic
 **/

int main(int argc, char *argv[])
{
    const auto begin = std::chrono::high_resolution_clock::now(); //Start the stopwatch for the total time.

/**
 * Read configs
 **/
    if (argc != 10)                                               //Very specific argument format for the program.
    {
        std::cout
            << "arguments format is: {volume_list} {conclusion output file name} {MaxTraffic} {Margin} {innerFilterSize} {outdir} {model time limit in seconds} {seed} {threads}"
            << std::endl;
        return 0;
    }
    volume_list = std::string(argv[1]);
    conclusionFileName = std::string(argv[2]);
    T_percentage = std::stod(std::string(argv[3]));
    Margin = std::stod(std::string(argv[4]));
    InnerFilterSize = std::stoi(std::string(argv[5]));
    std::string outDir = std::string(argv[6]);
    model_time_limit = std::stod(std::string(argv[7]));
    seed = std::string(argv[8]);
    number_of_threads = std::string(argv[9]);

    auto splitvolumelist = split_string(volume_list, "/");
    std::string write_solution = std::string(outDir) +  splitvolumelist[splitvolumelist.size() - 1] + std::string("_T") + std::to_string(T_percentage) + std::string("_M") + std::to_string(Margin) + std::string("_seed") + seed.c_str() + std::string("_TimeLimit") + std::string(argv[7]) + "_innerFilter" + std::string(argv[5]) + std::string(".csv");

    std::ifstream volume_list_f(volume_list.c_str(), std::ifstream::in);
    if (!volume_list_f.is_open())
    {
        std::cout << "error opening volume list." << std::endl;
        exit(1);
    }
	
    std::string content;

    std::vector<std::string> source_volume_list;
    std::vector<std::string> target_volume_list;
    std::vector<double> load_volume_list;

    /**
     * Read input volumes and desired balance
     **/
    while (std::getline(volume_list_f, content)) 
    {
        auto splitContent = split_string(content, ",");
        auto file_path = splitContent[0];
        auto volume_type = splitContent[1];
        double load = std::stod(splitContent[2]);
        
        if(volume_type.compare("a")) {
            load_volume_list.push_back(load);
        }
        if(volume_type.compare("s") || volume_type.compare("a")) {
            source_volume_list.push_back(file_path);
        }
        if(volume_type.compare("t") || volume_type.compare("a")) {
            target_volume_list.push_back(file_path);
        }
    }

    /**
     * Create intersects object
     **/
    std::vector<std::vector<std::vector<int>>> intersects_source_target_blocksn;
    for (int i = 0; i < source_volume_list.size(); i++) {        
        std::vector<std::vector<int>> intersects_i;
        intersects_source_target_blocksn.push_back(intersects_i);
        for (int j = 0; j < target_volume_list.size(); j++) {
            std::vector<int> intersect_i_j = calcVolumeIntersect(source_volume_list[i], target_volume_list[j]);
            intersects_source_target_blocksn[i].push_back(intersect_i_j);
        }
    }

    /**
     * Count and save amounts for easier future parsing
     **/
    std::vector<std::pair<int, int>> num_of_blocks_and_files_sourceVolumes;
    std::vector<std::pair<int, int>> num_of_blocks_and_files_targetVolumes;
    for (int i = 0; i < source_volume_list.size(); i++) {        
        int num_of_metadata_lines = get_num_of_metadata_lines(source_volume_list[i]);
        std::pair<int, int> num_of_blocks_and_files_i = get_num_of_blocks_and_files(source_volume_list[i], num_of_metadata_lines);
        num_of_blocks_and_files_sourceVolumes.push_back(num_of_blocks_and_files_i);
    }
    for (int j = 0; j < target_volume_list.size(); j++) {
        int num_of_metadata_lines = get_num_of_metadata_lines(target_volume_list[j]);
        std::pair<int, int> num_of_blocks_and_files_j = get_num_of_blocks_and_files(target_volume_list[j], num_of_metadata_lines);  
        num_of_blocks_and_files_targetVolumes.push_back(num_of_blocks_and_files_j);
    }

    std::pair<int, int> lastSourceSn_block_file = getLastBlockAndFileSn(source_volume_list);
    std::vector<double> block_sizes = getBlockSizes(source_volume_list, lastSourceSn_block_file);
    
    /* Calculate the global variable of the starting source size */
    SOURCE_SIZE = 0.0;
    for(int source = 0 ; source < intersects_source_target_blocksn.size(); source ++) {
        for(int i = 0; i < intersects_source_target_blocksn[source][source].size(); i++) {
            SOURCE_SIZE += block_sizes[intersects_source_target_blocksn[source][source][i]];
        }
    }
    std::cout << SOURCE_SIZE;
    
    std::vector<std::set<int>> fileSnInVolume = getFileSnInVolumes(source_volume_list);

    /**
     * Initiate global variables
     **/
    GRBEnv *env = 0;
    GRBVar *blocks_copied = 0; //Cist
    GRBVar *blocks_deleted = 0; //Dis
    GRBVar *files = 0;  //Xist
    GRBConstr *constrains = 0;
    GRBConstr *constrains_hint = 0;
    bool need_to_free_hint_constrains = false;
    std::vector<GRBLinExpr> left_side;
    std::vector<GRBLinExpr> left_side_hint; //Files that do not have blocks should stay at source

    std::vector<std::vector<GRBVar*>> C_i_s_t;
    std::vector<GRBVar*> D_i_s;
    std::vector<std::vector<GRBVar*>> X_l_s_t;

    try
    {
        env = new GRBEnv(); //This may throw and error if there is no valid licence.
        GRBModel model = GRBModel(*env);
        model.set(GRB_StringAttr_ModelName, "GoSeed");
        std::cout << model_time_limit;
        model.set(GRB_DoubleParam_TimeLimit, model_time_limit); //Set time limit

        model.set("Seed", seed.c_str());
        model.set("Threads", number_of_threads.c_str());

        for (int i = 0; i < block_sizes.size(); i++) {
            std::vector<GRBVar*> C_s_t;

            std::string blocks_deleted_names[source_volume_list.size()] ;
            char block_types_delete[source_volume_list.size()] ;
            
            for (int source = 0; source < source_volume_list.size(); source++) {
                std::string delete_name = std::string("D_") + std::to_string(i) + std::string("_") + std::to_string(source);
                blocks_deleted_names[source] = delete_name;
                block_types_delete[source] = GRB_BINARY;
                
                std::string blocks_copied_names[target_volume_list.size()] ;
                char block_types[target_volume_list.size()] ;
                
                for(int target = 0 ; target < target_volume_list.size() ; target ++) {
                    std::string name = std::string("C_") + std::to_string(i) + std::string("_") + std::to_string(source) + std::string("_") + std::to_string(target);
                    blocks_copied_names[target] = name;
                    block_types[target] = GRB_BINARY;
                }

                blocks_copied = model.addVars(NULL, NULL, NULL, block_types, blocks_copied_names, target_volume_list.size());
                C_s_t.push_back(blocks_copied);
            }
            C_i_s_t.push_back(C_s_t);

            blocks_deleted = model.addVars(NULL, NULL, NULL, block_types_delete, blocks_deleted_names, source_volume_list.size());
            D_i_s.push_back(blocks_deleted);
        }
        model.update();

        for (int l = 0; l <= lastSourceSn_block_file.second; l++) {
            std::vector<GRBVar*> X_s_t;
            for (int source = 0; source < source_volume_list.size(); source++) {
                
                std::string file_names[target_volume_list.size()] ;
                char block_types[target_volume_list.size()] ;
                
                for(int target = 0 ; target < target_volume_list.size() ; target ++) {
                    file_names[target] = std::string("X_") + std::to_string(l) + std::string("_") + std::to_string(source) + std::string("_") + std::to_string(target);
                    block_types[target] = GRB_BINARY;
                }
                files = model.addVars(NULL, NULL, NULL, block_types, file_names, target_volume_list.size());
                
                X_s_t.push_back(files);
            }
            X_l_s_t.push_back(X_s_t);  
        }
        model.update();
        
        // Check if load balancing is necessary. 
        // If margin is 1 or more, no load balancing is necessary (error margin of above 100%)
        if(Margin < 1.0) {
            addConstraint_allIntersectsAreCopied(model, intersects_source_target_blocksn, C_i_s_t);
            model.update();
            addConstraint_blockNeedsToExistToBeCloneOrDeleated(model, C_i_s_t, D_i_s, source_volume_list, target_volume_list);
            model.update();
            addConstraint_RemapFilesToOnlyOneVolume(model, X_l_s_t, num_of_blocks_and_files_sourceVolumes, target_volume_list.size());
            model.update();
            addConstraint_DontRemapNonExistantFilesOfRemapToSelf(model, X_l_s_t, fileSnInVolume, target_volume_list.size());
            model.update();
            add_Three_constraints_BlockIsDeletedOnlyIfNoLocalFileUsingItRemains_BlockIsDeletedOnlyIfNoNewFileTransferredIsUsingIt_RemmapedFileHasAllItsBlocks(model, X_l_s_t, C_i_s_t, D_i_s, source_volume_list, target_volume_list.size()); //one parsing of everything instead of 3
            model.update();
            addConstraint_doNotCopyFromMoreThanOneSource(model, C_i_s_t, source_volume_list, target_volume_list.size());
            model.update();
            addConstraint_loadBalancing(model, C_i_s_t, D_i_s, source_volume_list, load_volume_list, intersects_source_target_blocksn, block_sizes, Margin, InnerFilterSize);
            model.update();
            addConstraint_traffic_and_setObjective(model, C_i_s_t, D_i_s, intersects_source_target_blocksn, block_sizes, T_percentage);
            model.update();
            addConstraint_no_orphans(model, X_l_s_t, C_i_s_t, D_i_s, source_volume_list, target_volume_list.size(), lastSourceSn_block_file.first, lastSourceSn_block_file.second); 
            model.update();
            addConstraint_no_extraBlocks(model, X_l_s_t, C_i_s_t, D_i_s, source_volume_list, target_volume_list.size(), lastSourceSn_block_file.first, lastSourceSn_block_file.second); 
            model.update();
        } else {
            //The case where no load balncing is needed
            addConstraint_allIntersectsAreCopied(model, intersects_source_target_blocksn, C_i_s_t);
            model.update();
            addConstraint_blockNeedsToExistToBeCloneOrDeleated(model, C_i_s_t, D_i_s, source_volume_list, target_volume_list);
            model.update();
            addConstraint_RemapFilesToOnlyOneVolume(model, X_l_s_t, num_of_blocks_and_files_sourceVolumes, target_volume_list.size());
            model.update();
            addConstraint_DontRemapNonExistantFilesOfRemapToSelf(model, X_l_s_t, fileSnInVolume, target_volume_list.size());
            model.update();
            add_Three_constraints_BlockIsDeletedOnlyIfNoLocalFileUsingItRemains_BlockIsDeletedOnlyIfNoNewFileTransferredIsUsingIt_RemmapedFileHasAllItsBlocks(model, X_l_s_t, C_i_s_t, D_i_s, source_volume_list, target_volume_list.size()); //one parsing of everything instead of 3
            model.update();
            addConstraint_traffic_and_setObjectiveFullMigration(model, C_i_s_t, D_i_s, intersects_source_target_blocksn, block_sizes, T_percentage);
            model.update();
        }
        model.write("debud.lp");

        double modelCreationTime = ingestionTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - begin).count();
        auto s1 = std::chrono::high_resolution_clock::now();

        //This is where we call the gurobi solver. Everything else was setup everything after is summary.
        model.optimize();
        
        double solver_time = optimizationTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - s1).count();

        variable_number = model.get(GRB_IntAttr_NumVars);
        constraint_number = model.get(GRB_IntAttr_NumConstrs);

        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL)
        {
            solution_status = "OPTIMAL";
        }
        if (status == GRB_INFEASIBLE)
        {
            solution_status = "INFEASIBLE";
        }
        if (status == GRB_TIME_LIMIT)
        {
            solution_status = "TIME_LIMIT";
        }
        std::cout << "done optimization" << std::endl
                  << std::flush;
        if (solution_status != "INFEASIBLE")
        {
            Kbytes_to_replicate = actual_R_Kbytes;

            //Print the results.
            try
            {
                saveSolution(X_l_s_t, C_i_s_t, target_volume_list.size(), write_solution);
                saveMigrationPlan(source_volume_list, fileSnInVolume, X_l_s_t, write_solution+"_migrationPlan", target_volume_list.size());
                saveFilterBalance(write_solution, C_i_s_t, D_i_s, source_volume_list, block_sizes, intersects_source_target_blocksn);
                saveInput(write_solution, volume_list, T_percentage, Margin, InnerFilterSize, model_time_limit, seed);
                saveRunSolvingData(write_solution, modelCreationTime, solver_time, solution_status);
            }
            catch (...)
            {
                std::cout << "Exception at print_results, probably can't read variables" << std::endl;
                solution_status = "TIME_LIMIT_AT_PRESOLVE";
                recordFail(write_solution, solution_status);
            }
        } else {
            saveInput(write_solution, volume_list, T_percentage, Margin, InnerFilterSize, model_time_limit, seed);
            saveRunSolvingData(write_solution, modelCreationTime, solver_time, solution_status);
        }
        double elapsed_secs = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - begin).count();

        computeConclusionFields(volume_list,  C_i_s_t, D_i_s, source_volume_list, block_sizes, intersects_source_target_blocksn);
        save_runConclusionLine();
    } catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        recordFail(write_solution, e.getMessage());

    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
        recordFail(write_solution, "Exception during optimization");
    }

    delete[] blocks_copied;
    delete[] blocks_deleted;
    delete[] files;
    delete env;
    return 0;
}
