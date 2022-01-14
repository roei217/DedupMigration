#pragma once

#include <cstddef>
#include <memory>
#include <set>

class Node final{
public:
    /**
     * creates an initial node
     * @param file - file
     * @param size - the given file's size
     */
    explicit Node(const int file,  const double size);

    /**
     * the next 8 functions are getters/setters
     */
    bool isActivated() const;

    void disableNode();

    const std::set<int>& getCurrentFiles() const;

    void addFiles(const std::set<int>& new_files);

    void setSize(const double new_size);

    double getSize();

    void clearFiles();

    void reset();

//m_size - cluster size in bytes
//m_orig_size - used for resetting node to its orig state
//m_orig_file - used for resetting node to its orig state
//m_is_activated - is this not is activated
//m_current_files - set of all file currently in cluster
private:
    double m_size;
    const double m_orig_size;
    const int m_orig_file;
    bool m_is_activated;
    std::unique_ptr<std::set<int>> m_current_files;
};
