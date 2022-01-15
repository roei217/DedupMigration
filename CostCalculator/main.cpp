#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <sqlite3.h>
#include <tuple>
#include <sys/file.h>
#include <unistd.h>
#include "CommandLineParser.h"
#include "Lock.h"

using namespace std;

/**
    Main data structure explanation:
    _____________________________________
    vector<double> requested_volumes_sizes;
        A vector containing the requested final volume sizes in percentages

    vector<tuple<string, set<int>, set<int>, bool >> workloads:
        A vector of [workload's path, moved files, received files, use for system size] for each workload

    vector<tuple<string, string, string, string>> to_hash:
        A vector of [volume name, flattened init files, flattened final finals, use for system size] needed to hash

    vector<double> volumes_sizes:
        A vector containing the sizes of each volume in percentage

    map<int, long long int> fingerprints_to_move_to_me;
        A map containing all the blocks that need to be in this volume after migration, the mapping is between those blocks
        to their sizes.

    set<int> to_delete:
        A set of blocks that no longer need after the migration in this volume.

    set<int> received:
        A set of blocks which are needed after the migration in this volume.

    set<int> counted_for_total_size:
        A set of blocks that already been calculated in the system's final size.

    set<int> volume_blocks:
        A set of blocks that already been calculated in the volume's final size.

    Main variables explanation:
    _____________________________________
    tuple<double, double> cached_traffic_deletion:
        A tuple containing {-1,-1} in case this plan is not in cache, or {traffic, deletion} in case it is.
    long long int initial_system_size, final_system_size:
        An integer containing 0 in case this plan is not in cache, or initial/final system's size in case it is.
    vector<long long int> initial_volumes_sizes, volumes_traffic, volumes_deletion:
        An empty vector in case this plan is not in cache, or initial volume size/volume traffic/volume deletion for
        each volume in case it is.

*/

/**
    Global variables and their default values
*/

tuple<double, double> cached_traffic_deletion;
long long int initial_system_size = 0, final_system_size = 0;
vector<long long int> initial_volumes_sizes, volumes_traffic, volumes_deletion;


/**
    Utility general functions
*/

/**
 * splits string into token by delimiter
 * @param line - string to split
 * @param delimiter - delimiter
 * @return a vector of those tokens
 */
vector<string> split_string(const string &line, char delimiter) {
    stringstream ss(line);
    string item;
    vector<string> elements;
    while (getline(ss, item, delimiter)) {
        elements.push_back(item);
    }
    return elements;
}

/**
 * flattens a vector of ints to a string seperated by spaces
 * @param vec - int vector
 * @return - the resulted string
 */
string flat_vector(vector<int> &vec) {
    string output;
    for (int i = 0; i < vec.size(); ++i) {
        output += to_string(vec[i]) + " ";
    }
    return output.substr(0, output.length() - 1);
}

/**
 * converts size_t to string
 * @param x - size_t value
 * @return corresponding string
 */
string size_t_to_string(size_t x) {
    stringstream ss;
    ss << x;
    return ss.str();
}


/**
 * the cache class
 */
class cache {
    // db instance and attempts counter
    sqlite3 *db;
    string cache_path;

    /**
     * just opens the db in the path
     */
    void open_db() {
        char *zErrMsg = 0;
        int rc = sqlite3_open(cache_path.c_str(), &db);

        if (rc) {
            cerr << "Can't open database: " << sqlite3_errmsg(db) << endl;
            sqlite3_free(zErrMsg);
            sqlite3_close(db);
            exit(1);
        }
        sqlite3_free(zErrMsg);
    }

    /**
     * creates our calculator cache table
     */
    void create_table() {
        char *zErrMsg = 0;

        string sql = "CREATE TABLE IF NOT EXISTS calc_cache(volumes_hash text, init_files_hash text, final_files_hash text,"
                     " included text, Traffic double, Deletion double, Init_sys_size INTEGER , final_sys_size INTEGER, ";
        for (int i = 1; i < 21; i++) {
            sql += "init_vol" + to_string(i) + "_size INTEGER, vol" + to_string(i) + "_traffic INTEGER, vol" +
                   to_string(i) + "_deletion INTEGER, ";
        }

        sql += "PRIMARY KEY(volumes_hash, init_files_hash, final_files_hash, included))";

        /* Execute SQL statement */
        int rc = sqlite3_exec(db, sql.c_str(), NULL, 0, &zErrMsg);

        if (rc != SQLITE_OK) {
            cerr << "SQL error: " << zErrMsg << endl;
            sqlite3_free(zErrMsg);
            exit(1);
        }
        sqlite3_free(zErrMsg);
    }


public:
    /**
     * init the calculator cache
     */
    cache(string cache_path) : cache_path(cache_path) {
        open_db();
        create_table();
        sqlite3_close(db);
    }

    ~cache() {
    }

    /**
     * inserts a record to our calculator cache table
     * @param volumes - the hash of the volumes
     * @param init - the hash of the initial files
     * @param final - the hash of the final files
     * @param included - include in the final calculation
     * @param traffic - traffic result
     * @param deletion - deletion result
     * @param init_sys_size - the initial size of the system in bytes
     * @param final_sys_size - the final size of the system in bytes
     * @param init_vol_sizes - the initial sizes of each volume in bytes
     * @param vol_traffics - the resulted traffic sizes of each volume in bytes
     * @param vol_deletions - the resulted deletion sizes of each volume in bytes
     */
    void insert_result(string volumes, string init, string final, string included, double traffic, double deletion,
                       long long int init_sys_size, long long int final_sys_size, vector<long long int> &init_vol_sizes,
                       vector<long long int> &vol_traffics, vector<long long int> &vol_deletions) {
        //open db again, we might use concurrency, so we are closing it each time
        open_db();

        char *zErrMsg = 0;

        string sql = "INSERT INTO calc_cache(volumes_hash, init_files_hash, final_files_hash,"
                     " included, Traffic, Deletion, Init_sys_size, final_sys_size, ";
        for (int i = 1; i < 21; i++) {
            sql += "init_vol" + to_string(i) + "_size, vol" + to_string(i) + "_traffic, vol" + to_string(i) +
                   "_deletion, ";
        }

        sql = sql.substr(0, sql.length() - 2);
        sql += ") VALUES('" + volumes + "','" + init + "','" + final + "','" + included + "'," + to_string(traffic)
               + "," + to_string(deletion) + "," + to_string(init_sys_size) + "," + to_string(final_sys_size);

        for (int i = 1; i < 21; i++) {
            if (i > init_vol_sizes.size())
                sql += ",NULL,NULL,NULL";
            else {
                sql += "," +
                       to_string(init_vol_sizes[i - 1]) + "," +
                       to_string(vol_traffics[i - 1]) + "," +
                       to_string(vol_deletions[i - 1]);
            }
        }
        sql += ");";

        /* Execute SQL statement */
        int rc = sqlite3_exec(db, sql.c_str(), NULL, 0, &zErrMsg);

        if (rc != SQLITE_OK) {
            cerr << "SQL error: " << zErrMsg << endl;
        }
        sqlite3_free(zErrMsg);
        sqlite3_close(db);
    }

    /**
     * calllback for the select query, this is applied on each row returned from the select query
     */
    static int callback(void *data, int argc, char **argv, char **azColName) {
        get<0>(cached_traffic_deletion) = stod(string(argv[4]));//traffic
        get<1>(cached_traffic_deletion) = stod(string(argv[5]));//deletion
        initial_system_size = stoll(string(argv[6]));
        final_system_size = stoll(string(argv[7]));
        for (int i = 8; i < 8 + 3 * 20; i += 3) {
            if (argv[i])
                initial_volumes_sizes.push_back(stoll(string(argv[i])));
        }
        for (int i = 9; i < 9 + 3 * 20; i += 3) {
            if (argv[i])
                volumes_traffic.push_back(stoll(string(argv[i])));
        }
        for (int i = 10; i < 10 + 3 * 20; i += 3) {
            if (argv[i])
                volumes_deletion.push_back(stoll(string(argv[i])));
        }
        return 0;
    }

    /**
     * retrieves the result from the database, based on the input
     * @param volumes - the hash of the volumes
     * @param init - the hash of the initial files
     * @param final - the hash of the final files
     * @param included - include in the final calculation
     */
    void get_result(string volumes, string init, string final, string included) {
        open_db();

        char *zErrMsg = 0;
        cached_traffic_deletion = {-1, -1};

        string sql =
                "SELECT * FROM calc_cache WHERE volumes_hash='" + volumes + "' AND init_files_hash='" + init +
                "' AND final_files_hash='" + final + "' AND included='" + included + "'";

        /* Execute SQL statement */
        int rc = sqlite3_exec(db, sql.c_str(), callback, NULL, &zErrMsg);

        if (rc != SQLITE_OK) {
            cerr << "SQL error: " << zErrMsg << endl;
        }
        sqlite3_free(zErrMsg);
        sqlite3_close(db);
    }
};


/**
 * To evaluate the nature of a specific migration plan, we implemented a calculator which outputs the cost and requirements for this plan.
 * The calculator yields the traffic consumption, score obtained, and deletion acquired from this plan while also showing the initial and final size of each volume with specific traffic and deletion.
 * The calculator iterates over each workload, calculates its incoming traffic and resulted deletion.
 * The final traffic is the sum of all the traffics calculated while the final deletion is the sum of all traffic minus the sum of all deletions.
 * In an attempt to save time, we added a results cache which stores our results, in case a migration plan is required to be calculated again, it does not, and instead the results are fetched from the cache.
 * The calculator can only run concurrently on multiply threads using our file-based lock.
 *
 * Run with:
 * -output -  output file name
 * -file - a csv file with a row for every workload:
 * [workload full path, files in volume at start, files in volume at end, use for system size] - files separated with "-"
 * For example, a file might look like this:
 * row1: [volume1.csv,0-1-2-3-4,0-1-2-3-6-7, 1]
 * row2: [volume2.csv,5-6-7-8-9-10,5-8-9-10-4, 1]
 * row3: [volume2.csv,11-12-13-14-15,11-12-13-14-15, 0]
 * For a five workloads run, there should be five rows with no header.
 * -no_cache does not use the cache
 * -cache_path needs the full path of the desired db, such as '/nfs_share/cache.db'
 * -lb_sizes calculates the score based on predefined sizes, for example: 10 10 20 30 30 where each double is the
 *      requested size of a volume, if no such parameter then perfect lb is required
 * Note: currently 'use for system size' is set to 1 and is not used
 */
int main(int argc, char **argv) {
    // init parser for argv
    CommandLineParser parser(argc, argv);
    parser.addConstraint("-file", CommandLineParser::ArgumentType::STRING, 1, false, "migration file");
    parser.addConstraint("-output", CommandLineParser::ArgumentType::STRING, 1, true, "output file name");
    parser.addConstraint("-no_cache", CommandLineParser::ArgumentType::BOOL, 0, true, "use cache or not");
    parser.addConstraint("-cache_path", CommandLineParser::ArgumentType::STRING, 1, true, "path to cache");
    parser.addConstraint("-lb_sizes", CommandLineParser::ArgumentType::DOUBLE, -1, true, "sizes for score calculation");
    try {
        parser.validateConstraintsHold();
    } catch (const exception &e) {
        parser.printUsageAndDescription();
        throw;
    }

    //parse our args
    bool use_cache = !parser.isTagExist("-no_cache");

    string filename = parser.getTag("-file")[0];
    ifstream file;
    file.open(filename);
    if (!file.is_open()) {
        cerr << "error opening file" << endl;
        exit(1);
    }

    vector<double> requested_volumes_sizes;
    if (parser.isTagExist("-lb_sizes")) {
        vector<string> current_arg = parser.getTag("-lb_sizes");
        std::transform(current_arg.begin(), current_arg.end(), back_inserter(requested_volumes_sizes),
                       [](const string &x) {
                           long long int d = stoll(x);
                           return d;
                       });
        sort(requested_volumes_sizes.begin(), requested_volumes_sizes.end());
    }


    //vector of [file descriptor, moved files, received files, use for system size] for each workload
    vector<tuple<string, set<int>, set<int>, bool >> workloads;
    vector<string> splitted_line_content;
    string line_content;
    vector<tuple<string, string, string, string>> to_hash;//[volumes, init_files, final_finals, use for system size] needed to hash
    hash<string> hash_function;
    string db_path = parser.getTag("-cache_path")[0];
    cache calc_cache(db_path);
    Lock lock(db_path.substr(0, db_path.length() - 3) + ".txt");

    //parse input file workload by workload
    while (getline(file, line_content)) {
        //split the input line by ',' and parse the values in it
        splitted_line_content = split_string(line_content, ',');
        string workload_path = splitted_line_content[0];
        bool calculate_initial_system_size = splitted_line_content[3].find('1') != string::npos;
        vector<string> initial_files_string = split_string(splitted_line_content[1], '-');
        vector<string> final_files_string = split_string(splitted_line_content[2], '-');
        vector<int> initial_files, final_files;

        //receive the initial volume files
        transform(initial_files_string.begin(), initial_files_string.end(), back_inserter(initial_files),
                  [](const string &f) { return stoi(f); });
        set<int> initial_files_set(initial_files.begin(), initial_files.end());

        //receive the final volume files
        transform(final_files_string.begin(), final_files_string.end(), back_inserter(final_files),
                  [](const string &f) { return stoi(f); });
        set<int> final_files_set(final_files.begin(), final_files.end());

        //calculate the files that are moved from this volume and removed/transferred off it
        set<int> moved, received;
        for (int f : initial_files_set) {
            if (final_files_set.find(f) == final_files_set.end())
                moved.insert(f);
        }
        for (int f : final_files_set) {
            if (initial_files_set.find(f) == initial_files_set.end())
                received.insert(f);
        }

        //add as a workload
        workloads.emplace_back(workload_path, moved, received, calculate_initial_system_size);

        //hashing stuff - sets the parameters to hash in order to save in the cache
        vector<int> init_files_to_hash, final_files_to_hash;
        copy(initial_files_set.begin(), initial_files_set.end(), back_inserter(init_files_to_hash));
        copy(final_files_set.begin(), final_files_set.end(), back_inserter(final_files_to_hash));
        to_hash.emplace_back(split_string(workload_path, '/').back(), flat_vector(init_files_to_hash),
                             flat_vector(final_files_to_hash), calculate_initial_system_size ? "T" : "F");
    }

    // parsed the entire file, we sort it now to maintain an order
    file.close();
    sort(to_hash.begin(), to_hash.end());

    //concat each value (volumes, sizes, ...) in order to hash later
    string volumes, init_files, final_files, sizes;
    for (int i = 0; i < to_hash.size(); ++i) {
        volumes += get<0>(to_hash[i]) + " ";
        init_files += get<1>(to_hash[i]) + " ";
        final_files += get<2>(to_hash[i]) + " ";
        sizes += get<3>(to_hash[i]) + " ";
    }
    //remove back space
    volumes = volumes.substr(0, volumes.length() - 1);
    init_files = init_files.substr(0, init_files.length() - 1);
    final_files = final_files.substr(0, final_files.length() - 1);
    sizes = sizes.substr(0, sizes.length() - 1);;

    //if we are not using the cache, or we can't acquire the lock, or there are more than the maximum of 20 volumes, continue without a cache
    if (!use_cache || workloads.size() > 20 || !lock.lock());
    else {//otherwise, try to hash our parameters and find them in the db, eventually release the lock
        calc_cache.get_result(size_t_to_string(hash_function(volumes)), size_t_to_string(hash_function(init_files)),
                              size_t_to_string(hash_function(final_files)), size_t_to_string(hash_function(sizes)));
        lock.unlock();
    }

    //prepare an output file name in case one is not provided
    string output_filename = volumes + "_" + final_files + "_" + sizes + ".csv";
    replace(output_filename.begin(), output_filename.end(), ' ', '_');

    //stats of our calculation
    if (get<0>(cached_traffic_deletion) >= 0) {//traffic is not negative, found in cache (traffic, deletion)
        ofstream unit_output(
                parser.isTagExist("-output") ? parser.getTag("-output")[0] : output_filename);//use this filename
        if (!unit_output.is_open()) {
            cerr << "Error opening output file" << endl;
        }
        // in case we do not wanna use the cache or can not acquire the lock for some reason, calculate without it
        if (!use_cache || !lock.lock())
            goto not_found;


        //add to balancing file and output stuff
        double traffic = get<0>(cached_traffic_deletion), deletion = get<1>(cached_traffic_deletion);

        //output to our csv, first header
        unit_output
                << "workload,initial size B,initial size %,traffic B,deletion B,final size B,final size %,total traffic B,total traffic %,total deletion B,total deletion %,lb score"
                << endl;

        long long int traffic_bytes = 0, deletion_bytes = 0;
        vector<double> volumes_sizes;
        //now output our data to the output file by the header above
        for (int workload = 0; workload < workloads.size(); ++workload) {
            unit_output << get<0>(to_hash[workload]) << ",";
            if (get<3>(workloads[workload]))//relevant only for source
                unit_output << initial_volumes_sizes[workload] << ","
                            << initial_volumes_sizes[workload] * 100.0 / initial_system_size << ",";
            unit_output << volumes_traffic[workload] << ",";
            unit_output << volumes_deletion[workload] << ",";

            traffic_bytes += volumes_traffic[workload];//traffic is addable
            deletion_bytes += volumes_deletion[workload] - volumes_traffic[workload];//so is the deletion
            if (get<3>(workloads[workload])) {//relevant only for source
                unit_output << (initial_volumes_sizes[workload] + volumes_traffic[workload] -
                                volumes_deletion[workload])
                            << "," <<
                            (initial_volumes_sizes[workload] + volumes_traffic[workload] - volumes_deletion[workload]) *
                            100.0 / final_system_size << endl;
                double perc_from_sys =
                        (initial_volumes_sizes[workload] + volumes_traffic[workload] - volumes_deletion[workload]) *
                        100.0 / final_system_size;
                volumes_sizes.push_back(perc_from_sys);
            }
        }

        lock.unlock();

        //calculate the lb score, min/max
        double score = -1;
        if (!parser.isTagExist("-lb_sizes")) {
            score = *min_element(volumes_sizes.begin(), volumes_sizes.end()) /
                    *max_element(volumes_sizes.begin(), volumes_sizes.end());
        } else {
            std::sort(volumes_sizes.begin(), volumes_sizes.end());
            vector<double> normalized_volumes_sizes;

            for (int i = 0; i < volumes_sizes.size(); ++i) {
                normalized_volumes_sizes.push_back(
                        (100.0 / volumes_sizes.size()) / requested_volumes_sizes[i] * volumes_sizes[i]);
            }
            score = *min_element(normalized_volumes_sizes.begin(), normalized_volumes_sizes.end()) /
                    *max_element(normalized_volumes_sizes.begin(), normalized_volumes_sizes.end());

        }
        unit_output << ",,,,,,," << traffic_bytes << "," << traffic << ",";
        unit_output << deletion_bytes << "," << deletion << "," << score << endl;
        unit_output.close();
        return 0;
    }
    // did not find a result in cache
    not_found:
    //we look on every volume by itself, iteration for each
    long long int traffic = 0;
    long long int deletion = 0;
    for (int workload = 0; workload < workloads.size(); ++workload) {
        map<int, long long int> fingerprints_to_move_to_me;
        long long int temporary_traffic = 0;
        long long int temporary_deletion = 0;
        set<int> to_delete;
        set<int> received = get<2>(workloads[workload]);
        //go over all e the other volumes
        for (int other_workload = 0; other_workload < workloads.size(); ++other_workload) {
            if (workload == other_workload)
                continue;

            //open workload's file
            ifstream f(get<0>(workloads[other_workload]));
            if (!f.is_open()) {
                cerr << "error opening workload " << get<0>(workloads[other_workload]) << endl;
                exit(1);
            }

            //find blocks that need to move to us
            while (getline(f, line_content)) {//for each line
                splitted_line_content = split_string(line_content, ',');
                if (splitted_line_content[0] == "F" &&
                    received.find(stoi(splitted_line_content[1])) != received.end()) {//check only files that we receive
                    int number_of_blocks_in_file_line = stoi(splitted_line_content[4]);
                    //for each fp for this file if the fp was selected, add to relevant data structures
                    for (int block = 0; block < 2 * number_of_blocks_in_file_line; block += 2) {
                        int block_sn = stoi(splitted_line_content[5 + block]);
                        long long int size_read = stoll(splitted_line_content[6 + block]);
                        fingerprints_to_move_to_me[block_sn] = size_read;
                    }
                }
                if (splitted_line_content[0] == "B")
                    break;
            }
            f.close();
        }

        //open this workload
        ifstream f(get<0>(workloads[workload]));
        if (!f.is_open()) {
            cerr << "error opening workload" << get<0>(workloads[workload]) << endl;
            exit(1);
        }

        // sets the size of a fp that moves to it and already is there to 0
        while (getline(f, line_content)) {//for each line
            splitted_line_content = split_string(line_content, ',');
            if (splitted_line_content[0] == "B") {//check only blocks
                //for traffic
                int block_sn = stoi(splitted_line_content[1]);
                if (fingerprints_to_move_to_me.find(block_sn) !=
                    fingerprints_to_move_to_me.end())//already exists, no need for traffic
                    fingerprints_to_move_to_me[block_sn] = 0;
            }
        }

        //each size we didn't zero, needs to be added to traffic
        for (auto const &fingerprint : fingerprints_to_move_to_me) {
            temporary_traffic += fingerprint.second;
        }

        f.clear();
        f.seekg(0);//set to start of workload

        //calculate deletion, fp that is no longer needed in this volume and is not moved to it, counts as deletion
        set<int> moved = get<1>(workloads[workload]);
        while (getline(f, line_content)) {//for each line
            splitted_line_content = split_string(line_content, ',');
            if (splitted_line_content[0] == "B") {//check only blocks
                int block_sn = stoi(splitted_line_content[1]);
                int files_having_block = stoi(splitted_line_content[3]);
                int i;
                for (i = 0; i < files_having_block; ++i) {//try to find if the block still has a link to this volume
                    int file_sn = stoi(splitted_line_content[4 + i]);
                    if (moved.find(file_sn) == moved.end()) {//stayed
                        break;//still has a link to a file in this volume, next block
                    }
                }
                //no link to the files here, check if any migrating file has a link to it
                if (i == files_having_block &&
                    (fingerprints_to_move_to_me.empty() || fingerprints_to_move_to_me.find(block_sn) ==
                                                           fingerprints_to_move_to_me.end()))//it does not migrate to this volume - hence deleted
                    to_delete.insert(block_sn);

            }
        }
        f.clear();
        f.seekg(0);//set to start of workload

        //calculate the volume's final size and the total deletion of it
        set<int> counted_for_total_size, volume_blocks;
        long long int volume_size = 0;
        while (getline(f, line_content)) {//for each line
            splitted_line_content = split_string(line_content, ',');
            if (splitted_line_content[0] == "F") {//check only files that we receive
                bool file_moved = moved.find(stoi(splitted_line_content[1])) != moved.end();
                int number_of_blocks_in_file_line = stoi(splitted_line_content[4]);
                //for each fp for this file if the fp was selected, add to relevant data structures
                for (int block = 0; block < 2 * number_of_blocks_in_file_line; block += 2) {
                    int block_sn = stoi(splitted_line_content[5 + block]);
                    long long int size_read = stoll(splitted_line_content[6 + block]);
                    if (file_moved && to_delete.find(block_sn) != to_delete.end()) {
                        temporary_deletion += size_read;
                        to_delete.erase(block_sn);
                    }
                    if (get<3>(workloads[workload]) && counted_for_total_size.find(block_sn) ==
                                                       counted_for_total_size.end()) {//only if suppose to calculate for system size
                        counted_for_total_size.insert(block_sn);
                        initial_system_size += size_read;
                    }
                    if (volume_blocks.find(block_sn) == volume_blocks.end()) {//for general volume size
                        volume_blocks.insert(block_sn);
                        volume_size += size_read;
                    }
                }
            } else if (splitted_line_content[0] == "B")
                break;
        }
        counted_for_total_size.clear();

        f.close();

        //traffic is just added
        traffic += temporary_traffic;
        //deletion is added, but traffic brings data in, so subtract it
        deletion += temporary_deletion - temporary_traffic;

        //save it in the results vector
        initial_volumes_sizes.push_back(volume_size);
        volumes_traffic.push_back(temporary_traffic);
        volumes_deletion.push_back(temporary_deletion);

    }

    // open the output file
    ofstream unit_output(parser.isTagExist("-output") ? parser.getTag("-output")[0] : output_filename);
    if (!unit_output.is_open()) {
        cerr << "Error opening output file" << endl;
    }

    //calculate the entire system's final size
    final_system_size = initial_system_size;
    for (int workload = 0; workload < workloads.size(); ++workload) {
        final_system_size += volumes_traffic[workload] - volumes_deletion[workload];
    }

    //output our result to the output file, same as if we found the solution in cache
    unit_output
            << "workload,initial size B,initial size %,traffic B,deletion B,final size B,final size %,total traffic B,"
               "total traffic %,total deletion B,total deletion %,lb score"
            << endl;
    vector<double> volumes_sizes;
    for (int workload = 0; workload < workloads.size(); ++workload) {

        unit_output << get<0>(to_hash[workload]) << ",";
        if (get<3>(workloads[workload]))//relevant only for source
            unit_output << initial_volumes_sizes[workload] << ","
                        << initial_volumes_sizes[workload] * 100.0 / initial_system_size << ",";
        unit_output << volumes_traffic[workload] << ",";
        unit_output << volumes_deletion[workload] << ",";
        if (get<3>(workloads[workload])) {//relevant only for source
            unit_output << (initial_volumes_sizes[workload] + volumes_traffic[workload] - volumes_deletion[workload])
                        << ","
                        << (initial_volumes_sizes[workload] + volumes_traffic[workload] - volumes_deletion[workload]) *
                           100.0 /
                           final_system_size << endl;

            double perc_from_sys =
                    (initial_volumes_sizes[workload] + volumes_traffic[workload] - volumes_deletion[workload]) *
                    100.0 /
                    final_system_size;
            volumes_sizes.push_back(perc_from_sys);
        }

    }

    if (get<0>(cached_traffic_deletion) < 0) {//we need to insert our result, since we didn't find it in cache
        if (use_cache && lock.lock()) {
            try {
                calc_cache.insert_result(size_t_to_string(hash_function(volumes)),
                                         size_t_to_string(hash_function(init_files)),
                                         size_t_to_string(hash_function(final_files)),
                                         size_t_to_string(hash_function(sizes)),
                                         traffic * 100.0 / initial_system_size,
                                         deletion * 100.0 / initial_system_size,
                                         initial_system_size, final_system_size, initial_volumes_sizes,
                                         volumes_traffic, volumes_deletion);
            }
            catch (...) { ;//couldn't add it for some reason, just calculate it next time
            }
            lock.unlock();
        }
    }

    //calculate score
    double score = -1;
    if (!parser.isTagExist("-lb_sizes")) {
        score = *min_element(volumes_sizes.begin(), volumes_sizes.end()) /
                *max_element(volumes_sizes.begin(), volumes_sizes.end());
    } else {
        std::sort(volumes_sizes.begin(), volumes_sizes.end());
        vector<double> normalized_volumes_sizes;

        for (int i = 0; i < volumes_sizes.size(); ++i) {
            normalized_volumes_sizes.push_back(
                    (100.0 / volumes_sizes.size()) / requested_volumes_sizes[i] * volumes_sizes[i]);
        }
        score = *min_element(normalized_volumes_sizes.begin(), normalized_volumes_sizes.end()) /
                *max_element(normalized_volumes_sizes.begin(), normalized_volumes_sizes.end());

    }

    unit_output << ",,,,,,," << traffic << "," << traffic * 100.0 / initial_system_size << ",";
    unit_output << deletion << "," << deletion * 100.0 / initial_system_size << "," << score << endl;
    unit_output.close();

    return 0;
}
