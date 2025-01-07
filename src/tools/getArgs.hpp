
#ifndef GETARGS_HPP
#define GETARGS_HPP

#include <string>

#include <cstring>
#include <tuple>

class argsController {
    int argc;
    char ** argv;
public:
    argsController(int a, char ** b) {
        argc = a;
        argv = b;
    }

    std::string get(const std::string & s) {
        for(int i = 0; i < argc; i++) {
            if(s == argv[i]) {
                return i + 1 == argc ? "" : argv[i+1];
            }
        }
        return "";
    }

    std::string operator [](const std::string & s) {
        return get(s);
    }

    bool exist(const std::string & s) {
        for(int i = 0; i < argc; i++) {
            if(s == argv[i]) return true;
        }
        return false;
    }
};

#endif