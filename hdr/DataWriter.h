#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <filesystem>
#include "global.h"

namespace fs = std::filesystem;
using std::ofstream;
using std::string;
using std::vector;

struct DataWriter
{
public:
    DataWriter(string pathName_);
    fs::path createTimeDirectory(double time);
    void writeData(vector<macroParam> data, double time);
    void setDelta_h(double dh_);
private:
    double dh = 1;
    fs::path directory;
};
