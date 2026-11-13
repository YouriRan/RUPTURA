#pragma once
#include <string>
#include <vector>

#include "json.hpp"
namespace cycle
{
using Json = nlohmann::ordered_json;
using json = Json;
using Vec = std::vector<double>;
}  // namespace cycle
