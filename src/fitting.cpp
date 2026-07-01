#include "fitting.h"

#include <algorithm>
#include <bitset>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <print>
#include <sstream>
#include <unordered_set>

#include "json.h"
#include "random_numbers.h"
#include "special_functions.h"
#if __cplusplus >= 201703L && __has_include(<filesystem>)
#include <filesystem>
#elif __cplusplus >= 201703L && __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
#else
#include <sys/stat.h>
#endif

Fitting::Fitting(const InputReader& inputreader)
    : numberOfComponents(inputreader.components.size()),
      components(inputreader.components),
      displayName(inputreader.displayName),
      componentName(numberOfComponents),
      filename(numberOfComponents),
      isotherms(numberOfComponents),
      columnPressure(inputreader.columnPressure - 1),
      columnLoading(inputreader.columnLoading - 1),
      columnError(inputreader.columnError - 1),
      pressureScale(PressureScale(inputreader.pressureScale)),
      geneticAlgorithmSize(static_cast<size_t>(std::pow(2.0, 12.0))),
      geneticAlgorithmMutationRate(1.0 / 3.0),
      geneticAlgorithmEliteRate(0.15),
      geneticAlgorithmMotleyCrowdRate(0.25),
      geneticAlgorithmDisasterRate(0.001),
      geneticAlgorithmElitists(
          static_cast<size_t>(static_cast<double>(geneticAlgorithmSize) * geneticAlgorithmEliteRate)),
      geneticAlgorithmMotleists(
          static_cast<size_t>(static_cast<double>(geneticAlgorithmSize) * (1.0 - geneticAlgorithmMotleyCrowdRate))),
      alphaPopulation(static_cast<size_t>(std::pow(2.0, 12.0))),
      betaPopulation(static_cast<size_t>(std::pow(2.0, 12.0))),
      parents(alphaPopulation),
      children(betaPopulation)
{
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    componentName[i] = inputreader.components[i].name;
    filename[i] = inputreader.components[i].filename;
    isotherms[i] = inputreader.components[i].isotherm;
  }
}

bool startsWith(const std::string& str, const std::string& prefix)
{
  return str.size() >= prefix.size() && str.substr(0, prefix.size()) == prefix;
}

std::string trim(const std::string& s)
{
  auto start = s.begin();
  while (start != s.end() && std::isspace(*start))
  {
    start++;
  }

  auto end = s.end();
  do
  {
    end--;
  } while (std::distance(start, end) > 0 && std::isspace(*end));

  return std::string(start, end + 1);
}

namespace
{
std::string isothermTypeName(Isotherm::Type type)
{
  switch (type)
  {
    case Isotherm::Type::Langmuir:
      return "Langmuir";
    case Isotherm::Type::Anti_Langmuir:
      return "Anti-Langmuir";
    case Isotherm::Type::BET:
      return "BET";
    case Isotherm::Type::Henry:
      return "Henry";
    case Isotherm::Type::Freundlich:
      return "Freundlich";
    case Isotherm::Type::Sips:
      return "Sips";
    case Isotherm::Type::Langmuir_Freundlich:
      return "Langmuir-Freundlich";
    case Isotherm::Type::Redlich_Peterson:
      return "Redlich-Peterson";
    case Isotherm::Type::Toth:
      return "Toth";
    case Isotherm::Type::Unilan:
      return "Unilan";
    case Isotherm::Type::OBrien_Myers:
      return "O'Brien & Myers";
    case Isotherm::Type::Quadratic:
      return "Quadratic";
    case Isotherm::Type::Temkin:
      return "Temkin";
    case Isotherm::Type::BingelWalton:
      return "Bingel&Walton";
    default:
      throw std::runtime_error("Error: unknown isotherm type");
  }
}
}  // namespace

void Fitting::readData(size_t ID)
{
  std::ifstream fileInput{filename[ID]};
  std::string errorOpeningFile = "File '" + filename[ID] + "' exists, but error opening file";
  if (!fileInput) throw std::runtime_error(errorOpeningFile);

  std::print("Reading: {}\n", filename[ID]);

  std::string line{};

  maximumLoading = 0.0;
  rawData.clear();
  while (std::getline(fileInput, line))
  {
    std::string trimmedLine = trim(line);
    if (!startsWith(trimmedLine, "#"))
    {
      if (!line.empty())
      {
        std::istringstream iss(line);

        std::vector<std::string> results((std::istream_iterator<std::string>(iss)),
                                         std::istream_iterator<std::string>());
        if (columnPressure < results.size() && columnLoading < results.size())
        {
          double pressure;
          double loading;
          std::istringstream s(results[columnPressure]);
          s >> pressure;
          std::istringstream t(results[columnLoading]);
          t >> loading;
          if (loading > maximumLoading)
          {
            maximumLoading = loading;
          }
          rawData.push_back(std::make_pair(pressure, loading));
        }
      }
    }
  }

  if (rawData.empty())
  {
    throw std::runtime_error("Error: no pressure points found");
  }

  // sort the pressures
  std::sort(rawData.begin(), rawData.end());

  pressureRange = std::make_pair(rawData.front().first, rawData.back().first);
  logPressureRange = std::make_pair(std::log(pressureRange.first), std::log(pressureRange.second));

  std::print("Found {} data points\n", rawData.size());
  for (const std::pair<double, double>& data : rawData)
  {
    std::print("{} {}\n", data.first, data.second);
  }
  std::print("\n");
  std::print("Lowest pressure: {}\n", pressureRange.first);
  std::print("Highest pressure: {}\n", pressureRange.second);
  std::print("Log lowest pressure: {}\n", logPressureRange.first);
  std::print("Log highest pressure: {}\n", logPressureRange.second);

  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    std::print("Number of isotherm parameters: {}\n", isotherms[i].numberOfParameters);
    isotherms[i].print();
  }
}

void Fitting::writeComponentsJson(const std::string& path) const
{
  using nlohmann::ordered_json;

  ordered_json outputJson;
  outputJson["Components"] = ordered_json::array();

  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    const Component& component = components[i];
    const MultiSiteIsotherm& isotherm = isotherms[i];

    ordered_json componentJson;
    componentJson["Name"] = component.name;
    componentJson["GasPhaseMolFraction"] = component.initialGasMoleFraction;

    if (component.isCarrierGas)
    {
      componentJson["CarrierGas"] = true;
    }
    else
    {
      componentJson["MassTransferCoefficient"] = component.massTransferCoefficient;
      componentJson["AxialDispersionCoefficient"] = component.axialDispersionCoefficient;
      componentJson["NumberOfIsothermSites"] = isotherm.numberOfSites;
      componentJson["IsothermSites"] = ordered_json::array();

      for (size_t siteIndex = 0; siteIndex < isotherm.numberOfSites; ++siteIndex)
      {
        const Isotherm& site = isotherm.sites[siteIndex];

        ordered_json siteJson;
        siteJson["Type"] = isothermTypeName(site.type);
        siteJson["Parameters"] = ordered_json::array();

        for (size_t parameterIndex = 0; parameterIndex < site.numberOfParameters; ++parameterIndex)
        {
          siteJson["Parameters"].push_back(site.parameters[parameterIndex]);
        }

        componentJson["IsothermSites"].push_back(siteJson);
      }
    }

    outputJson["Components"].push_back(componentJson);
  }

  std::ofstream output(path);
  if (!output)
  {
    throw std::runtime_error("Error: could not open JSON output file '" + path + "'");
  }

  std::print(output, "{}\n", outputJson.dump(2));

  std::print("Wrote fitted components JSON: {}\n", path);
}

void Fitting::run()
{
  std::print("STARTING FITTING\n");

  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    if (components[i].isCarrierGas)
    {
      continue;
    }

    readData(i);

    const DNA bestCitizen = fit(i);
    const DNA optimizedBestCitizen = simplex(bestCitizen, 1.0);

    optimizedBestCitizen.phenotype.print();

    isotherms[i] = optimizedBestCitizen.phenotype;
    components[i].isotherm = optimizedBestCitizen.phenotype;

    // createPlotScripts(optimizedBestCitizen, i);
  }

  const std::string outputFilename = displayName.empty() ? "fitted.json" : displayName + "_fitted.json";

  writeComponentsJson(outputFilename);

  // createPlotScript();
}

// create a new citizen in the Ensemble
Fitting::DNA Fitting::newCitizen(size_t ID)
{
  DNA citizen;

  citizen.phenotype = isotherms[ID].randomized(maximumLoading);

  citizen.genotype.clear();
  citizen.genotype.reserve((sizeof(double) * CHAR_BIT) * citizen.phenotype.numberOfParameters);
  for (size_t i = 0; i < citizen.phenotype.numberOfParameters; ++i)
  {
    // convert from double to bitset
    uint64_t p;
    std::memcpy(&p, &citizen.phenotype.parameters(i), sizeof(double));
    std::bitset<sizeof(double) * CHAR_BIT> bitset(p);

    // add the bit-string to the genotype representation
    citizen.genotype += bitset.to_string();
  }

  citizen.hash = std::hash<MultiSiteIsotherm>{}(citizen.phenotype);
  citizen.fitness = fitness(citizen.phenotype);

  return citizen;
}

void Fitting::updateCitizen(DNA& citizen) { citizen.fitness = fitness(citizen.phenotype); }

inline bool my_isnan(double val)
{
  union
  {
    double f;
    uint64_t x;
  } u = {val};
  return (u.x << 1) > (0x7ff0000000000000u << 1);
}

double Fitting::fitness(const MultiSiteIsotherm& phenotype)
// For evaluating isotherm goodness-of-fit:
// Residual Root Mean Square Error (RMSE)
{
  double fitnessValue = phenotype.fitness();
  size_t m = rawData.size();                // number of observations
  size_t p = phenotype.numberOfParameters;  // number of adjustable parameters
  for (std::pair<double, double> dataPoint : rawData)
  {
    double pressure = dataPoint.first;
    double loading = dataPoint.second;
    double difference = loading - phenotype.value(pressure, 1.0);
    // double weight = 1.0/(1.0+loading);
    double weight = 1.0;
    fitnessValue += weight * difference * difference;
  }
  fitnessValue = sqrt(fitnessValue / static_cast<double>(m - p));

  if (my_isnan(fitnessValue)) fitnessValue = 99999999.999999;
  if (fitnessValue == 0.0000000000) fitnessValue = 99999999.999999;

  return fitnessValue;
}

double Fitting::RCorrelation(const MultiSiteIsotherm& phenotype)
{
  double RCorrelationValue = phenotype.fitness();
  size_t m = rawData.size();
  double loading_avg_o = 0.0;
  double loading_avg_e = 0.0;
  double tmp1 = 0.0;
  double tmp2 = 0.0;
  double tmp3 = 0.0;

  for (std::pair<double, double> dataPoint : rawData)
  {
    double pressure = dataPoint.first;
    double loading = dataPoint.second;
    loading_avg_o += loading / static_cast<double>(m);
    loading_avg_e += phenotype.value(pressure, 1.0) / static_cast<double>(m);
  }

  for (std::pair<double, double> dataPoint : rawData)
  {
    double pressure = dataPoint.first;
    double loading = dataPoint.second;
    tmp1 += (loading - loading_avg_o) * (phenotype.value(pressure, 1.0) - loading_avg_e);
    tmp2 += (loading - loading_avg_o) * (loading - loading_avg_o);
    tmp3 += (phenotype.value(pressure, 1.0) - loading_avg_e) * (phenotype.value(pressure, 1.0) - loading_avg_e);
  }
  RCorrelationValue = tmp1 / sqrt(tmp2 * tmp3);

  return RCorrelationValue;
}

size_t Fitting::biodiversity(const std::vector<DNA>& citizens)
{
  std::map<size_t, size_t> counts;
  for (const DNA& dna : citizens)
  {
    if (counts.find(dna.hash) != counts.end())
    {
      ++counts[dna.hash];
    }
    else
    {
      counts[dna.hash] = 1;
    }
  }
  size_t biodiversity = 0;
  for (const std::pair<size_t, size_t> value : counts)
  {
    if (value.second > 1)
    {
      biodiversity += value.second;
    }
  }

  return biodiversity;
}

void Fitting::nuclearDisaster(size_t ID)
{
  for (size_t i = 1; i < children.size(); ++i)
  {
    children[i] = newCitizen(ID);
  }
}

void Fitting::elitism()
{
  std::copy(parents.begin(), parents.begin() + static_cast<std::vector<DNA>::difference_type>(geneticAlgorithmElitists),
            children.begin());
}

void Fitting::mutate(DNA& mutant)
{
  mutant.genotype.clear();
  mutant.genotype.reserve((sizeof(double) * CHAR_BIT) * mutant.phenotype.numberOfParameters);
  for (size_t i = 0; i < mutant.phenotype.numberOfParameters; ++i)
  {
    // convert from double to bitset
    uint64_t p;
    std::memcpy(&p, &mutant.phenotype.parameters(i), sizeof(double));
    std::bitset<sizeof(double) * CHAR_BIT> bitset(p);

    // mutation: randomly flip bit
    bitset.flip(std::size_t((sizeof(double) * CHAR_BIT) * RandomNumber::Uniform()));

    // convert from bitset to double
    p = bitset.to_ullong();
    std::memcpy(&mutant.phenotype.parameters(i), &p, sizeof(double));

    // add the bit-string to the genotype representation
    mutant.genotype += bitset.to_string();
  }

  // calculate the hah-value from the entire bit-string
  mutant.hash = std::hash<MultiSiteIsotherm>{}(mutant.phenotype);
}

// [s1:s2] range of children
// [i1:i2] range of parent1
// [j1:j2] range of parent2
// One-point crossover
//----------------------------------------------
//  parent1    parent2                children
//    *          *                      *
//  00|000000  11|111111        ->    00|111111
//----------------------------------------------
// Two-point
//----------------------------------------------
//  parent1     parent2               children
//    *    *      *    *                *    *
//  00|0000|00  11|1111|11       ->   00|1111|00
//----------------------------------------------
void Fitting::crossover(size_t ID, size_t s1, size_t s2, size_t i1, size_t i2, size_t j1, size_t j2)
{
  size_t k1, k2;
  double tmp1;
  for (size_t i = s1; i < s2; ++i)
  {
    chooseRandomly(i1, i2, j1, j2, k1, k2);
    tmp1 = RandomNumber::Uniform();
    // choose between single cross-over using bit-strings or random parameter-swap
    if (tmp1 < 0.490)
    // One-point crossover:
    // --------------------
    {
      // remove the extreme values 0 and 32*Npar - 1 (they are not valid for crossover)
      size_t bitStringSize = (sizeof(double) * CHAR_BIT) * isotherms[ID].numberOfParameters;
      size_t spos = RandomNumber::Integer(1, bitStringSize - 2);
      children[i].genotype =
          parents[k1].genotype.substr(0, spos) + parents[k2].genotype.substr(spos, bitStringSize - spos);

      // convert the bit-strings to doubles
      for (size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
      {
        size_t pos = j * (sizeof(double) * CHAR_BIT);
        size_t size = sizeof(double) * CHAR_BIT;
        std::bitset<sizeof(double) * CHAR_BIT> bitset(children[i].genotype, pos, size);
        uint64_t p = bitset.to_ullong();
        std::memcpy(&children[i].phenotype.parameters(j), &p, sizeof(double));
      }
    }
    else if (tmp1 < 0.499)
    // Two-point crossover:
    // --------------------
    {
      size_t bitStringSize = (sizeof(double) * CHAR_BIT) * isotherms[ID].numberOfParameters;
      size_t spos1 = RandomNumber::Integer(1, bitStringSize - 3);
      size_t spos2 = RandomNumber::Integer(spos1, bitStringSize - 2);
      children[i].genotype = parents[k1].genotype.substr(0, spos1) + parents[k2].genotype.substr(spos1, spos2 - spos1) +
                             parents[k1].genotype.substr(spos2, bitStringSize - spos2);
      // convert the bit-strings to doubles
      for (size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
      {
        size_t pos = j * (sizeof(double) * CHAR_BIT);
        size_t size = sizeof(double) * CHAR_BIT;
        std::bitset<sizeof(double) * CHAR_BIT> bitset(children[i].genotype, pos, size);
        uint64_t p = bitset.to_ullong();
        std::memcpy(&children[i].phenotype.parameters(j), &p, sizeof(double));
      }
    }
    else if (tmp1 < 0.500)
    {
      // Uniform crossover:
      // ------------------
      size_t bitStringSize = (sizeof(double) * CHAR_BIT) * isotherms[ID].numberOfParameters;
      size_t rolling_k = k1;
      for (size_t j = 0; j < bitStringSize; j++)
      {
        if (RandomNumber::Uniform() < 0.25)
        {
          if (rolling_k == k1)
          {
            rolling_k = k2;
          }
          else
          {
            rolling_k = k1;
          }
        }
        children[i].genotype.substr(j, 1) = parents[rolling_k].genotype.substr(j, 1);
      }
      // convert the bit-strings to doubles
      for (size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
      {
        size_t pos = j * (sizeof(double) * CHAR_BIT);
        size_t size = sizeof(double) * CHAR_BIT;
        std::bitset<sizeof(double) * CHAR_BIT> bitset(children[i].genotype, pos, size);
        uint64_t p = bitset.to_ullong();
        std::memcpy(&children[i].phenotype.parameters(j), &p, sizeof(double));
      }
    }
    else
    {
      children[i].genotype.clear();
      for (size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
      {
        // randomly choose whether the parameter comes from parent k1 or k2
        if (RandomNumber::Uniform() < 0.5)
        {
          children[i].phenotype.parameters(j) = parents[k1].phenotype.parameters(j);
        }
        else
        {
          children[i].phenotype.parameters(j) = parents[k2].phenotype.parameters(j);
        }

        // convert from double to bitString
        uint64_t p;
        std::memcpy(&p, &children[i].phenotype.parameters(j), sizeof(double));
        std::bitset<sizeof(double) * CHAR_BIT> bitset(p);
        children[i].genotype += bitset.to_string();
      }
    }

    children[i].hash = std::hash<MultiSiteIsotherm>{}(children[i].phenotype);
  }
}

void Fitting::chooseRandomly(size_t kk1, size_t kk2, size_t jj1, size_t jj2, size_t& ii1, size_t& ii2)
{
  ii1 = RandomNumber::Integer(kk1, kk2);
  ii2 = RandomNumber::Integer(jj1, jj2);
  while (ii1 == ii2)
  {
    ii2 = RandomNumber::Integer(jj1, jj2);
  };
}

void Fitting::mate(size_t ID)
{
  // retain the first 25% of the children
  elitism();

  // mates from geneticAlgorithmElitists to (geneticAlgorithmSize - geneticAlgorithmElitists)
  crossover(ID, geneticAlgorithmElitists,
            geneticAlgorithmElitists + static_cast<size_t>(static_cast<double>(geneticAlgorithmMotleists) * 0.5), 0,
            geneticAlgorithmElitists, 0, geneticAlgorithmElitists);
  crossover(ID,
            geneticAlgorithmElitists + static_cast<size_t>(static_cast<double>(geneticAlgorithmMotleists) * 0.5) + 1,
            geneticAlgorithmSize - geneticAlgorithmElitists, 0, geneticAlgorithmElitists, geneticAlgorithmElitists,
            geneticAlgorithmSize - 1);

  // mutation from geneticAlgorithmElitists to (geneticAlgorithmSize - geneticAlgorithmElitists) with
  // "geneticAlgorithmMutationRate" probability
  for (size_t i = geneticAlgorithmElitists; i < geneticAlgorithmSize - geneticAlgorithmElitists; ++i)
  {
    if (RandomNumber::Uniform() < geneticAlgorithmMutationRate)
    {
      mutate(children[i]);
    }
    updateCitizen(children[i]);
  }

  // replace the last geneticAlgorithmElitists (the worst) of the children by new children
  for (size_t i = geneticAlgorithmSize - geneticAlgorithmElitists; i < geneticAlgorithmSize; ++i)
  {
    children[i] = newCitizen(ID);
  }

  // replace the last (geneticAlgorithmSize - 1) children by new children
  if (RandomNumber::Uniform() < geneticAlgorithmDisasterRate)
  {
    nuclearDisaster(ID);
  }
}

bool DNA_Fitness_Sorter(Fitting::DNA const& lhs, Fitting::DNA const& rhs) { return lhs.fitness < rhs.fitness; }

void Fitting::sortByFitness() { std::sort(parents.begin(), parents.end(), &DNA_Fitness_Sorter); }

void Fitting::writeCitizen(size_t citizen, size_t id, size_t step, size_t variety, size_t fullfilledCondition)
{
  char info[256];
  if (fullfilledCondition > 0)
  {
    snprintf(info, 256,
             "mol: %2ld  step: %5ld  Fitness: %10.6lf R^2: %10.6lf Similarity: %5ld/%-5ld Finishing: %3ld/%-3d\n", id,
             step, parents[citizen].fitness, pow(RCorrelation(parents[citizen].phenotype), 2), variety,
             geneticAlgorithmSize, fullfilledCondition, 100);
  }
  else
  {
    snprintf(info, 256, "mol: %2ld  step: %5ld  Fitness: %10.6lf R^2: %10.6lf Similarity: %5ld/%-5ld\n", id, step,
             parents[citizen].fitness, pow(RCorrelation(parents[citizen].phenotype), 2), variety, geneticAlgorithmSize);
  }
  std::print("{}", info);
  std::print("number of parameters: {}\n", parents[citizen].phenotype.numberOfParameters);
  for (size_t i = 0; i < parents[citizen].phenotype.numberOfParameters; ++i)
  {
    std::print("      genotype: {} parameter: {}\n", parents[citizen].genotype.substr(64 * i, 64),
               parents[citizen].phenotype.parameters(i));
  }
  std::print("\n");
}

Fitting::DNA Fitting::fit(size_t ID)
{
  size_t optimisationStep{0};
  const size_t maxOptimisationStep{1000};
  size_t fullFilledConditionStep{0};
  const size_t maxFullfilledConditionStep{100};
  const double minimumFitness{5.0e-1};
  double tempFitnessValue{999.0};
  size_t tempVarietyValue{0};
  const double toleranceEqualFitness{1e-3};
  const size_t minstep{10};

  fullFilledConditionStep = 0;
  optimisationStep = 0;

  for (size_t i = 0; i < alphaPopulation.size(); ++i)
  {
    alphaPopulation[i] = newCitizen(ID);
    betaPopulation[i] = newCitizen(ID);
  }

  parents = alphaPopulation;
  children = betaPopulation;

  sortByFitness();

  // print Initial (and unsorted) population
  writeCitizen(0, ID, 0, 0.0, 0);

  if (refittingFlag)
  {
    std::print("Refitting activated\n");
    for (size_t citizen = 0; citizen < 2; ++citizen)
    {
      parents[citizen].genotype.clear();
      parents[citizen].genotype.reserve((sizeof(double) * CHAR_BIT) * parents[citizen].phenotype.numberOfParameters);
      for (size_t i = 0; i < parents[citizen].phenotype.numberOfParameters; ++i)
      {
        // convert from double to bitset
        uint64_t p;
        std::memcpy(&p, &parents[citizen].phenotype.parameters(i), sizeof(double));
        std::bitset<sizeof(double) * CHAR_BIT> bitset(p);

        // add the bit-string to the genotype representation
        parents[citizen].genotype += bitset.to_string();
      }
      parents[citizen].hash = std::hash<MultiSiteIsotherm>{}(parents[citizen].phenotype);
      updateCitizen(parents[citizen]);
    }
    std::copy(parents.begin(), parents.end(), children.begin());

    mate(ID);

    std::swap(parents, children);

    sortByFitness();
  }

  isotherms[ID] = parents[0].phenotype;

  std::print("Starting Genetic Algorithm optimization\n");

  bool continueCondition = true;
  do
  {
    sortByFitness();

    tempVarietyValue = biodiversity(children);

    writeCitizen(0, ID, optimisationStep, tempVarietyValue, fullFilledConditionStep);

    if (optimisationStep >= minstep && parents[0].fitness <= minimumFitness &&
        std::abs(parents[0].fitness - tempFitnessValue) <= toleranceEqualFitness)
    {
      fullFilledConditionStep += 1;
    }
    else
    {
      fullFilledConditionStep = 0;
    }

    if (optimisationStep >= maxOptimisationStep || fullFilledConditionStep >= maxFullfilledConditionStep)
    {
      continueCondition = false;
    }

    // pairing
    mate(ID);

    std::swap(parents, children);

    // take the best fitness value:
    tempFitnessValue = parents[0].fitness;

    optimisationStep += 1;

  } while (continueCondition);

  writeCitizen(0, ID, optimisationStep, tempVarietyValue, fullFilledConditionStep);

  return parents[0];
}

// The Nelder-Mead method uses a simplex (a hyper-tetrahedron of n+1 vertices in n dimensions)
// Advantage: it does not use derivates, works well, can tolerate some noise
// Disadvantage: it is not garanteed to converge
// The steps of the method are:
// 1) Sort: according to the fitness
// 2) Reflect: get rid of the worst point, replace it by something better
// 3) Extend: if better then extend it even further
// 4) Contract: if not better then contract it
// 5) Shrink: if still not better we shrink towards the best performing point
// 6) Check convergence
const Fitting::DNA Fitting::simplex(DNA citizen, double scale)
{
  size_t n = citizen.phenotype.numberOfParameters;
  std::vector<std::vector<double>> v(n + 1, std::vector<double>(n));  // holds vertices of simplex
  std::vector<double> f(n + 1);                                       // value of function at each vertex
  std::vector<double> vr(n);                                          // reflection - coordinates
  std::vector<double> ve(n);                                          // expansion - coordinates
  std::vector<double> vc(n);                                          // contraction - coordinates
  std::vector<double> vm(n);                                          // centroid - coordinates
  std::vector<double> vtmp(n);                                        // temporary array passed to FUNC
  size_t vs;                                                          // vertex with the smallest value
  size_t vh;                                                          // vertex with next smallest value
  size_t vg;                                                          // vertex with largest value
  double fr;                                                          // value of function at reflection point
  double fe;                                                          // value of function at expansion point
  double fc;                                                          // value of function at contraction point
  size_t iprint{0};
  size_t MAX_IT = 1000000;
  double EPSILON = 1.0e-4;
  double ALPHA = 1.0;
  double BETA = 0.5;
  double GAMMA = 2.0;

  std::print("\nMinimising the cost function using the Nelder-Mead SIMPLEX method:\n\n");

  for (size_t i = 0; i < n; ++i)
  {
    v[0][i] = citizen.phenotype.parameters(i);
  }

  // values used to create initial simplex
  double pn = scale * (std::sqrt(static_cast<double>(n) + 1.0) - 1.0 + static_cast<double>(n)) /
              (static_cast<double>(n) * sqrt(2.0));
  double qn = scale * (std::sqrt(static_cast<double>(n) + 1.0) - 1.0) / (static_cast<double>(n) * sqrt(2.0));
  for (size_t i = 1; i <= n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      if (i - 1 == j)
      {
        v[i][j] = pn + citizen.phenotype.parameters(j);
      }
      else
      {
        v[i][j] = qn + citizen.phenotype.parameters(j);
      }
    }
  }

  for (size_t i = 0; i <= n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      citizen.phenotype.parameters(j) = v[i][j];
    }
    f[i] = fitness(citizen.phenotype);
  }

  // print out the initial simplex
  // print out the initial function values
  // find the index of the smallest value for printing
  if (iprint == 0)
  {
    vs = 0;
    for (size_t j = 0; j <= n; ++j)
    {
      if (f[j] < f[vs])
      {
        vs = j;
      }
    }

    std::print("Initial Values from genetic algorithm:\n");

    for (size_t j = 0; j < n; ++j)
    {
      std::print("{} ", v[vs][j]);
    }
    std::print("Fit: {}\n\n", f[vs]);
  }

  for (size_t itr = 1; itr <= MAX_IT; ++itr)
  {
    // Step 1: Sort
    // ====================================================================
    std::vector<size_t> sortIndexes = sort_indexes(f);
    vs = sortIndexes[0];      // index of smallest
    vg = sortIndexes[n];      // index of largest
    vh = sortIndexes[n - 1];  // index of second largest

    // calculate the center point of every point except for the worst one
    for (size_t j = 0; j < n; ++j)
    {
      double cent = 0.0;
      for (size_t i = 0; i <= n; ++i)
      {
        if (i != vg)
        {
          cent = cent + v[i][j];
        }
      }
      vm[j] = cent / static_cast<double>(n);
    }

    // Step 2: Reflect vg to new vertex vr
    // ====================================================================
    for (size_t j = 0; j < n; ++j)
    {
      vr[j] = (1.0 + ALPHA) * vm[j] - ALPHA * v[vg][j];
      citizen.phenotype.parameters(j) = vr[j];
    }
    fr = fitness(citizen.phenotype);

    if ((fr <= f[vh]) && (fr > f[vs]))
    {
      for (size_t j = 0; j < n; ++j)
      {
        v[vg][j] = vr[j];
      }
      f[vg] = fr;
    }

    // Step 3: Extend a step further in this direction
    // ====================================================================
    if (fr <= f[vs])
    {
      for (size_t j = 0; j < n; ++j)
      {
        ve[j] = GAMMA * vr[j] + (1.0 - GAMMA) * vm[j];
        citizen.phenotype.parameters(j) = ve[j];
      }
      fe = fitness(citizen.phenotype);

      // by making fe < fr as opposed to fe < f(vs), Rosenbrocks function
      // takes 62 iterations as opposed to 64.

      if (fe < fr)
      {
        for (size_t j = 0; j < n; ++j)
        {
          v[vg][j] = ve[j];
        }
        f[vg] = fe;
      }
      else
      {
        for (size_t j = 0; j < n; ++j)
        {
          v[vg][j] = vr[j];
        }
        f[vg] = fr;
      }
    }

    if (fr > f[vh])
    {
      // Step 4: Contraction
      // ====================================================================
      for (size_t j = 0; j < n; ++j)
      {
        vc[j] = BETA * v[vg][j] + (1.0 - BETA) * vm[j];
        citizen.phenotype.parameters(j) = vc[j];
      }
      fc = fitness(citizen.phenotype);
      if (fc < f[vg])
      {
        for (size_t j = 0; j < n; ++j)
        {
          v[vg][j] = vc[j];
        }
        f[vg] = fc;
      }
      else
      {
        // Step 4: Shrink
        // ====================================================================
        // at this point the contraction is not successful,
        // we must halve the distance from vs to all the
        // vertices of the simplex and then continue.
        for (size_t row = 0; row <= n; ++row)
        {
          if (row != vs)
          {
            for (size_t j = 0; j < n; ++j)
            {
              v[row][j] = v[vs][j] + 0.5 * (v[row][j] - v[vs][j]);
            }
          }
        }
        for (size_t m = 0; m < n; ++m)
        {
          vtmp[m] = v[vg][m];
          citizen.phenotype.parameters(m) = vtmp[m];
        }
        f[vg] = fitness(citizen.phenotype);

        for (size_t m = 0; m < n; ++m)
        {
          vtmp[m] = v[vh][m];
          citizen.phenotype.parameters(m) = vtmp[m];
        }
        f[vh] = fitness(citizen.phenotype);
      }
    }

    // Step 6: Test for convergence
    // ====================================================================
    double fsum = 0.0;
    for (size_t j = 0; j <= n; ++j)
    {
      fsum = fsum + f[j];
    }
    double favg = fsum / (static_cast<double>(n) + 1.0);

    if (favg < EPSILON || itr == MAX_IT)
    {
      // print out the value at each iteration

      if (itr != MAX_IT)
      {
        std::print("Nelder-Mead has converged: {} < {}\n\n", favg, EPSILON);
      }
      else
      {
        std::print("Reached maximum number of steps: {} = {}\n\n", itr, MAX_IT);
      }

      // find the index of the smallest value
      sortIndexes = sort_indexes(f);
      vs = sortIndexes[0];  // index of smallest
      for (size_t m = 0; m < n; ++m)
      {
        citizen.phenotype.parameters(m) = v[vs][m];
      }
      double min = fitness(citizen.phenotype);

      std::print("Final Values: \n");
      for (size_t j = 0; j < n; ++j)
      {
        std::print("{} ", v[vs][j]);
      }
      std::print("Fit: {} R2: {}\n\n", min, pow(RCorrelation(citizen.phenotype), 2));

      return citizen;
    }
  }

  return citizen;
}
