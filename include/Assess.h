#pragma once
#include <algorithm>
#include <vector>
#include <string>
#include <map>


#include "UniformityResult.h"
void Evaluate_Results (const std::string &generatorName, const std::map<std::string, std::vector<double>> &pvalues, double ALPHA);
void Summarize_pvalue_Series (const std::string &testName, const std::vector<double> &values, double ALPHA);
UniformityResult calculate_uniformity_p_value (const std::vector<double> &pvalues);
