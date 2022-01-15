#include "Lock.h"
#include <stdexcept>

#include <thread>
#include <chrono>

#include <sys/file.h>
#include <unistd.h>

Lock::Lock(std::string lock_path):path(lock_path), fd(-1){}

bool Lock::lock() {
    fd = open(this->path.c_str(), O_RDWR | O_CREAT, 0666);
    if (fd < 0) {
        throw std::runtime_error("Cannot open file lock");
    }
    int try_num = 0;
    while (try_num < 30 && flock(fd, LOCK_EX | LOCK_NB)) {
        try_num++;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
    if (try_num == 30) {
        close(fd);
        fd = -1;
        return false;
    }
    return true;
}

void Lock::unlock() {
    if (fd >= 0) {
        flock(fd, LOCK_UN);
        close(fd);
        fd = -1;
    }
}


