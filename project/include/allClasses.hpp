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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utilities.hpp"
#include "domain.hpp"
#include "linearSystem.hpp"
#include "solvers.hpp"
#include "multigrid.hpp"

#endif