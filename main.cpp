#include <iostream>
#include "init_data.h"
#include "model.h"

int main(int argc, char *argv[]) {
    std::string path = "..\\ËãÀı\\Ëæ»úËãÀı\\";
    std::string file = "5-3-1(0).txt";
    TXTDatareader(path + file);

    modelSolver(path, file);

    return 0;
}


