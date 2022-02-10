/**
Created by Roei Kisous
*/
#include "CommandLineParser.hpp"

#include <iostream>
#include <string>

CommandLineParser::CommandLineParser(int argc, char **argv){
    for (int i = 1; i < argc; ++i) {
        const string current_arg = argv[i];
        if (current_arg[0] == '-') {
            vector<string> current_tag_value;
            while (i < argc - 1 && argv[i+1][0] != '-')
            {
                current_tag_value.emplace_back(argv[i+1]);
                i++;
            }

            m_args.insert(make_pair(current_arg, current_tag_value));
        }
    }
}

void CommandLineParser::validateConstraintsHold() const{
    for (const ConstraintTuple& constraint : m_constraints)
        validateConstraintHolds(constraint);
}

void CommandLineParser::validateConstraintHolds(const ConstraintTuple& constraint) const{
    if (!constraint.is_optional & !isTagExist(constraint.tag))
        throw std::invalid_argument("ERROR: The param " + constraint.tag + " is missing");

    if(isTagExist(constraint.tag))
    {
        const vector<string>& tag_values = getTag(constraint.tag);
        if(tag_values.size() != constraint.occurrences && constraint.occurrences != -1)
            throw std::invalid_argument("ERROR: The param " + constraint.tag + " Need to receive " +
            to_string(constraint.occurrences) + " parameters");

        validateConstraintValueType(constraint, tag_values);
    }
}

vector<string> CommandLineParser::getTag(const string& option) const{
    const auto found = m_args.find(option);
    if (found != m_args.end())
        return found->second;

    throw std::invalid_argument("Tag does not exist");
}

void CommandLineParser::validateConstraintValueType(const ConstraintTuple& constraint, const vector<string>& values){
    switch (constraint.type) {
        case ArgumentType::INT:
            for (const auto& value : values){
                if (!isNumber(value))
                    throw std::invalid_argument("ERROR: Need to receive int for " + constraint.tag);
            }
            break;
        case ArgumentType::DOUBLE:
            for (const auto &value : values) {
                try {
                    stod(value);
                }
                catch (const std::exception &e){
                    throw std::invalid_argument("ERROR: Need to receive double for " + constraint.tag);
                }
            }
            break;

        default:
            break;
    }
}

bool CommandLineParser::isNumber(const string& s) {
    return !s.empty() && std::find_if(s.cbegin(),s.cend(),[](unsigned char c) { return !std::isdigit(c); }) == s.cend();
}

bool CommandLineParser::isTagExist(const string &option) const {
    return m_args.find(option) != m_args.cend();
}

void CommandLineParser::addConstraint(const string &tag, ArgumentType type, int occurrences, bool optional,
                                      const string & description) {
    m_constraints.emplace_back(tag, type, occurrences, optional, description);
}

void CommandLineParser::printUsageAndDescription() const{
    cout << "[USAGE]:"<<endl<<"./hc ";

    for (const auto& constraint : m_constraints) {
        string constraint_str = constraint.tag;
        if(constraint.type != ArgumentType::BOOL)
        {
            const string type_str = ArgumentTypeToString(constraint.type) +
                    (constraint.occurrences == CommandLineParser::VARIABLE_NUM_OF_OCCURRENCES ?
                        " - VARIABLE LENGTH LIST" : "*" + to_string(constraint.occurrences));

            constraint_str +=  "(TYPE=" + type_str +")";
        }

        constraint_str = constraint.is_optional? "[" + constraint_str + "] " : constraint_str +" ";
        cout << constraint_str;
    }

    cout <<endl << endl;

    printDescription();
}

string CommandLineParser::ArgumentTypeToString(const ArgumentType arg_type){
    switch (arg_type) {
        case ArgumentType::INT:
            return "INT";
        case ArgumentType::STRING:
            return "STRING";
        case ArgumentType::DOUBLE:
            return "DOUBLE";
        default:
            return "";
    }
}

void CommandLineParser::printDescription() const{
    cout <<"PARAMS DESCRIPTION:" <<endl;

    for (const auto& constraint : m_constraints)
        cout<< "\t" << constraint.tag << ": " << constraint.description << endl;

    cout <<endl;

}
