#pragma once

#include <algorithm>
#include <map>
#include <vector>
#include <string>

using namespace std;

class  CommandLineParser final
{
public:
    static constexpr int VARIABLE_NUM_OF_OCCURRENCES = -1;
    enum class ArgumentType{
        INT=0,
        STRING,
        DOUBLE,
        BOOL
    };

    struct ConstraintTuple{
        ConstraintTuple(const string &tag, ArgumentType type, int occurrences, bool optional, const string &description):
            tag(tag), type(type), occurrences(occurrences), is_optional(optional), description(description)
        {
        }

        string tag;
        string description;
        ArgumentType type;
        int occurrences;
        bool is_optional;
    };

    using TagValuesMap = map<string, vector<string>>;

public:
    /**
     * parses argv
     * @param argc - program's argc
     * @param argv - program's argv
     */
    explicit CommandLineParser(int argc, char **argv);

    /**
     * adds a constraint for command line
     * @param tag - in form of '-' and something, such as -t
     * @param type - INT/STRING/DOUBLE/BOOL
     * @param occurrences - number of occurrences after the tag (-1 for any size)
     * @param optional - is this tag optional or no
     * @param description - string describes the param, default is empty string
     */
    void addConstraint(const string& tag, const ArgumentType type, const int occurrences = 0,
                       const bool optional = true, const string & description="");

    /**
     * prints the USAGE and description
     */
    void printUsageAndDescription() const;

    /**
     * prints params description
     */
    void printDescription() const;

    /**
     * enforce those constrains
     * @throw std::invalid_argument if some constraint is not valid
     */
    void validateConstraintsHold() const;

    /**
     *
     * @param tag - command line tag
     * @return the corresponding arguments
     */
    vector<string> getTag(const string &tag) const;

    /**
     *
     * @param tag - command line tag
     * @return whether the tag exists
     */
    bool isTagExist(const string& tag) const;

private:
    /**
     * @param constraint - a given constraint
     * @throw std::invalid_argument if constraint is not valid
     */
    void validateConstraintHolds(const ConstraintTuple& constraint) const;

private:
    /**
     * @param constraint - a given constraint
     * @param values - a params list given for the constraint
     * @throw std::invalid_argument if the params list is not valid for the constraint
     */
    static void validateConstraintValueType(const ConstraintTuple& constraint, const vector<string>& values);

    /**
     * @param arg_type - a ArgumentType enum value
     * @return string represent the enum value
     */
    static string ArgumentTypeToString(const ArgumentType arg_type);

    /**
     * @param s - a string
     * @return whether the string represents a number
     */
    static bool isNumber(const string& s);

private:
    TagValuesMap m_args;
    vector<ConstraintTuple> m_constraints;
};
