//
// Created by Hugin on 2025.12.22.
//

#ifndef BISHE_MODEL_H
#define BISHE_MODEL_H

#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <tuple>
#include <numeric>
#include <map>
#include <gurobi_c++.h>
#include "init_data.h"

typedef std::vector<GRBVar> GRBVector;
typedef std::vector<GRBVector> GRBMatrix;
typedef std::vector<GRBMatrix> GRBCube;
typedef std::vector<GRBCube> GRBTensor4d;


void modelSolver(const std::string& path, const std::string& logname);


#endif //BISHE_MODEL_H
