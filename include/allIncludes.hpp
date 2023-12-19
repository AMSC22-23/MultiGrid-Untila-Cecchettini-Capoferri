#ifndef ALL_H
#define ALL_H

#include <iostream>
#include <cmath>
#include <tuple>
#include <vector>
#include <functional>
#include <numeric>
#include <memory>
#include <fstream>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif


#include "utilities.hpp"
#include "domain.hpp"
#include "linear_system.hpp"
#include "solvers.hpp"
#include "multigrid.hpp"
#include "Main.hpp"

#endif