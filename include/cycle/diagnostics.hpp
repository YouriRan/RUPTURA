#pragma once
#include "cycle/types.hpp"
namespace cycle
{
// Descriptive statistics only: cycle variation is not a confidence interval.
Json sample_statistics(const Vec& values);
Json assess_convergence(const Json& rows, const Json& config);
}  // namespace cycle
