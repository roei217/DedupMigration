#ifndef LOCK_H
#define LOCK_H

#include "string"

class Lock {
public:
    bool lock();

    void unlock();

    Lock(std::string lock_path);

private:
    int fd;
    std::string path;
};

#endif // LOCK_H
