#include "HierarchicalClustering.hpp"
#include "CommandLineParser.hpp"

#include <iostream>
#include <cassert>

using namespace std;

static void setUpAndValidateParser(CommandLineParser& parser){
    parser.addConstraint("-workloads", CommandLineParser::ArgumentType::STRING, CommandLineParser::VARIABLE_NUM_OF_OCCURRENCES, false,
                         " workloads to load, list of strings");
    parser.addConstraint("-fps", CommandLineParser::ArgumentType::STRING, 1, false,
                         "number of min hash fingerprints (input 'all' for all fps) - int or 'all'");
    parser.addConstraint("-WT", CommandLineParser::ArgumentType::DOUBLE, CommandLineParser::VARIABLE_NUM_OF_OCCURRENCES, false,
                         "WT - list of int");
    parser.addConstraint("-lb", CommandLineParser::ArgumentType::BOOL,0, true,
                         "use load balancing (optional, default false)");
    parser.addConstraint("-seed", CommandLineParser::ArgumentType::INT, CommandLineParser::VARIABLE_NUM_OF_OCCURRENCES, false,
                         "seeds for the algorithm, list of int");
    parser.addConstraint("-gap", CommandLineParser::ArgumentType::DOUBLE, CommandLineParser::VARIABLE_NUM_OF_OCCURRENCES, false,
                         "pick random value from values in this gap - list of double");
    parser.addConstraint("-lb_sizes", CommandLineParser::ArgumentType::DOUBLE, CommandLineParser::VARIABLE_NUM_OF_OCCURRENCES, true,
                         "a list of clusters' requested sizes - list of int, must sum to 100 (optional, default is even distribution)");
    parser.addConstraint("-eps", CommandLineParser::ArgumentType::INT, 1, true,
                         "% to add to system's initial size at every iteration - int (mandatory if -lb is used, default is 5)");
    parser.addConstraint("-clusters", CommandLineParser::ArgumentType::INT, 1, true,
                         "number of clusters to output (optional, default is same number of input workloads)");
    parser.addConstraint("-output_path_prefix", CommandLineParser::ArgumentType::STRING, CommandLineParser::VARIABLE_NUM_OF_OCCURRENCES, true,
                         "output path prefix (optional, default is results/result in the working directory)");

    try {
        parser.validateConstraintsHold();
    }catch (const exception& e){
        parser.printUsageAndDescription();
        throw;
    }
}

static vector<string> validateAndGetSortedWorkloadsPaths(const CommandLineParser& parser){
    vector<string> workloads_paths = parser.getTag("-workloads");

    //sort in lexicographic order for convenience
    sort(workloads_paths.begin(), workloads_paths.end(), [](string &a, string &b) {
        return Utility::splitString(a, '/').back() < Utility::splitString(b, '/').back();
    });

    //check files exist
    for(const string& path : workloads_paths){
        if (!Utility::isFileExists(path))
            throw invalid_argument("The given workload's file path=" + path+ " doesn't exist");
    }

    return workloads_paths;
}

static int validateAndGetRequestedNumOfFingerprints(const CommandLineParser& parser){

    const vector<string> current_arg = parser.getTag("-fps");
    try {
        return stoi(current_arg.front());
    } catch (invalid_argument &e) {
        if (current_arg.front() != "all")
            throw invalid_argument("Please enter a valid fps_size");

        //default is set to 'all'
        static constexpr int ALL_FINGERPRINTS = -1;
        return ALL_FINGERPRINTS;
    }
}

static vector<double> validateAndGetTraffics(const CommandLineParser& parser){

    vector<double> traffics;

    for(const string& wt: parser.getTag("-WT")){
        try {
            const double converted_wt = stod(wt);
            if(converted_wt < 0 || converted_wt > 1)
                throw invalid_argument("WT value should be between 0 to 1");

            traffics.emplace_back(converted_wt * 100);
        } catch (invalid_argument &e) {
            throw invalid_argument("WT value should be between 0 to 1");
        }
    }

    return traffics;
}

static int validateAndGetNumOfClusters(const CommandLineParser& parser, const vector<string>& workloads_paths){
    if (parser.isTagExist("-clusters"))
        return stoi(parser.getTag("-clusters").front());

    //default
    return workloads_paths.size();
}

static int validateAndGetEps(const CommandLineParser& parser){
    if (parser.isTagExist("-eps"))
        return stoi(parser.getTag("-eps").front());

    static constexpr int EPS_DEFAULT = 5;
    return EPS_DEFAULT;
}

static vector<int> validateAndGetSeeds(const CommandLineParser& parser){
    vector<int> seeds;
    for(const string& seed: parser.getTag("-seed"))
        seeds.emplace_back(stoi(seed));

    return seeds;
}

static vector<double> validateAndGetGaps(const CommandLineParser& parser){
    vector<double> gaps;
    for(const string& gap: parser.getTag("-gap"))
        gaps.emplace_back(stod(gap));

    return gaps;
}

static vector<double> validateAndGetSortedLbSizes(const CommandLineParser& parser, const int number_of_clusters){
    if (!parser.isTagExist("-lb_sizes"))
        return vector<double>(number_of_clusters, 100.0/number_of_clusters);

    vector<double> lb_sizes;
    for(const string& lb_size : parser.getTag("-lb_sizes")){
        try {
            const double converted_lb_size = stod(lb_size);
            if(converted_lb_size < 0 || converted_lb_size > 100)
                throw invalid_argument("lb size values should be between 0 to 100");

            lb_sizes.emplace_back(converted_lb_size);
        } catch (invalid_argument &e) {
            throw invalid_argument("lb size values should be between 0 to 100");
        }
    }

    sort(lb_sizes.begin(), lb_sizes.end(), greater<double>());

    return lb_sizes;
}

static string validateAndGetOutputPath(const CommandLineParser& parser){
    static const string DEFAULT_OUTPUT_PATH_PREFIX = "results/result";
    if(!parser.isTagExist("-output_path_prefix"))
        return DEFAULT_OUTPUT_PATH_PREFIX;

    return parser.getTag("-output_path_prefix").front();
}

/**
 * 1. -workloads: workloads to load, list of strings
 * 2. -fps: number of min hash fingerprints (input 'all' for all fps) - int or 'all'
 * 3. -WT: WT - list of doubles
 * 4. -lb: use load balancing (optional, default false)
 * 5. -seed: seeds for the algorithm, list of int
 * 6. -gap: pick random value from values in this gap - list of double
 * 7. -lb_sizes: a list of clusters' requested sizes - list of int, must sum to 100 (optional, default is even distribution)
 * 8. -eps: % to add to system's initial size at every iteration - int (mandatory if -lb is used, default is 5)
 * 9. -clusters: number of clusters to output (optional, default is same number of input workloads)
 * 10. -output_path_prefix: output path prefix (optional, default is results/result in the working directory)
 */
int main(int argc, char **argv) {
    try {
        CommandLineParser parser(argc, argv);
        setUpAndValidateParser(parser);

        //parsing arguments
        const bool load_balance = parser.isTagExist("-lb");
        const vector<string> workloads_paths = validateAndGetSortedWorkloadsPaths(parser);
        const int requested_number_of_fingerprints = validateAndGetRequestedNumOfFingerprints(parser);
        const int number_of_clusters = validateAndGetNumOfClusters(parser, workloads_paths);
        const int eps = validateAndGetEps(parser);
        const string output_path = validateAndGetOutputPath(parser);
        const vector<double> traffics = validateAndGetTraffics(parser);
        const vector<int> seeds = validateAndGetSeeds(parser);
        const vector<double> gaps = validateAndGetGaps(parser);
        const vector<double> lb_sizes = validateAndGetSortedLbSizes(parser, number_of_clusters);

        //init matrices
        unique_ptr<AlgorithmDSManager> DSManager = make_unique<AlgorithmDSManager>(workloads_paths,
                                                                                   requested_number_of_fingerprints,
                                                                                   load_balance);

        //run HC
        HierarchicalClustering HC(DSManager, lb_sizes);

        HC.run(load_balance, eps, number_of_clusters, traffics, seeds, gaps, output_path);

        return EXIT_SUCCESS;
    }catch (const exception& e){
        cerr << "Got exception: "<< e.what()<< endl;
        exit(EXIT_FAILURE);
    }
}
