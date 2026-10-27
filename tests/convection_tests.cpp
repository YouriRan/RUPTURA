#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "column.h"
#include "component.h"
#include "inputreader.h"
#include "integrators/compute.h"
#include "json.h"
#include "temporary_directory.h"

namespace
{
std::string examplePath(const std::string& relativePath)
{
  return (std::filesystem::path{RUPTURA_SOURCE_DIR} / relativePath).string();
}

// values[0] at the inlet face z_0, values[i] the average of f over the cell (z_{i-1}, z_i], for f(z) = a + b z.
std::vector<double> linearCellValues(const std::vector<double>& distances, double a, double b)
{
  std::vector<double> values(distances.size());
  values[0] = a + b * distances[0];
  for (size_t i = 1; i < distances.size(); ++i) values[i] = a + b * 0.5 * (distances[i - 1] + distances[i]);
  return values;
}

std::vector<double> uniformDistances(size_t intervals, double length)
{
  std::vector<double> distances(intervals + 1);
  for (size_t i = 0; i <= intervals; ++i) distances[i] = length * static_cast<double>(i) / static_cast<double>(intervals);
  return distances;
}

// Largest face error of the uniform WENO3 reconstruction of exp(z) on [0, 1] over the faces z_0 ... z_{N-1}.
double exponentialFaceError(size_t intervals)
{
  const std::vector<double> distances = uniformDistances(intervals, 1.0);
  const double h = distances[1];
  std::vector<double> values(intervals + 1);
  values[0] = 1.0;
  for (size_t i = 1; i <= intervals; ++i) values[i] = (std::exp(distances[i]) - std::exp(distances[i - 1])) / h;
  std::vector<double> faces(intervals + 1);
  computeWENO(values, distances, faces);
  double error = 0.0;
  for (size_t i = 0; i < intervals; ++i) error = std::max(error, std::abs(faces[i] - std::exp(distances[i])));
  return error;
}

Component inertComponent(size_t id, double axialDispersion)
{
  return Component(id, "C" + std::to_string(id), std::vector<Isotherm>{}, 0.0, 0.0, axialDispersion);
}

// A single adsorbent over the whole column, carrying the components (only their axial dispersion is used).
std::vector<MixturePrediction> singleAdsorbent(const std::vector<Component>& components)
{
  return {MixturePrediction("inert", components, 0, 0, 300.0, 1.0e3, 1.0e5, 2, 0, 0, 0)};
}

nlohmann::json beaInput()
{
  std::ifstream input(examplePath("examples/BEA-alkanes-C7/breakthrough/simulation.json"));
  return nlohmann::json::parse(input);
}
}  // namespace

TEST(WENO3Reconstruction, IsExactForLinearProfilesOnUniformAndNonUniformGrids)
{
  const std::vector<double> uniform = uniformDistances(8, 0.4);
  const std::vector<double> stretched{0.0, 0.01, 0.025, 0.05, 0.06, 0.1, 0.16, 0.2, 0.21};
  for (const std::vector<double>& distances : {uniform, stretched})
  {
    const std::vector<double> values = linearCellValues(distances, 0.3, -1.7);
    std::vector<double> faces(values.size());
    computeWENO(values, distances, faces);
    // The outlet face is first order (last cell value); every other face is exact.
    for (size_t i = 0; i + 1 < values.size(); ++i) EXPECT_NEAR(faces[i], 0.3 - 1.7 * distances[i], 1.0e-14) << i;
  }
}

TEST(WENO3Reconstruction, IsThirdOrderForASmoothProfile)
{
  const double coarse = exponentialFaceError(20);
  const double medium = exponentialFaceError(40);
  const double fine = exponentialFaceError(80);
  EXPECT_GT(coarse / medium, 6.0);
  EXPECT_GT(medium / fine, 6.0);
}

TEST(WENO3Reconstruction, KeepsAStepWithinItsBounds)
{
  std::vector<double> values(21, 0.0);
  std::fill(values.begin(), values.begin() + 8, 1.0);
  std::vector<double> faces(values.size());
  computeWENO(values, uniformDistances(20, 1.0), faces);
  for (double face : faces)
  {
    EXPECT_GE(face, -1.0e-12);
    EXPECT_LE(face, 1.0 + 1.0e-12);
  }
  EXPECT_NEAR(faces[7], 1.0, 1.0e-12);  // last cell before the step
  EXPECT_NEAR(faces[8], 0.0, 1.0e-12);  // first cell after the step
}

TEST(MassBalanceWENO3, GasSpeciesFluxesSumToTheUpwindTotalFlux)
{
  // The velocity closures conserve the upwind total molar flux v_i c_tot,i, so the WENO3 species fluxes must sum to
  // it: sum_j dc_j/dt = -(v_i c_tot,i - v_{i-1} c_tot,i-1) / dz without dispersion and sinks.
  constexpr size_t gridPoints = 12;
  constexpr size_t componentCount = 3;
  constexpr double resolution = 0.01;
  const std::vector<Component> components{inertComponent(0, 0.0), inertComponent(1, 0.0), inertComponent(2, 0.0)};
  std::vector<double> velocity(gridPoints + 1);
  std::vector<double> concentration((gridPoints + 1) * componentCount);
  for (size_t grid = 0; grid <= gridPoints; ++grid)
  {
    const double z = static_cast<double>(grid) / static_cast<double>(gridPoints);
    velocity[grid] = 0.1 - 0.03 * z;
    const double front = 0.5 + 0.5 * std::tanh(15.0 * (0.4 - z));
    concentration[grid * componentCount + 0] = 30.0 + 2.0 * z;
    concentration[grid * componentCount + 1] = 5.0 * front;
    concentration[grid * componentCount + 2] = 3.0 * (1.0 - front) + 0.2;
  }
  std::vector<double> concentrationDot(concentration.size());
  const std::vector<double> sink(concentration.size(), 0.0);
  const std::vector<double> ones(gridPoints + 1, 1.0);
  computeMassDerivatives(singleAdsorbent(components), gridPoints, componentCount, 1,
                         uniformDistances(gridPoints, gridPoints * resolution), ones, velocity, concentration,
                         concentrationDot, sink, Column::FluidPhase::Gas, ones, Column::ConvectionScheme::WENO3);

  auto total = [&](size_t grid)
  {
    double sum = 0.0;
    for (size_t comp = 0; comp < componentCount; ++comp) sum += concentration[grid * componentCount + comp];
    return sum;
  };
  for (size_t grid = 1; grid <= gridPoints; ++grid)
  {
    double sumDot = 0.0;
    for (size_t comp = 0; comp < componentCount; ++comp) sumDot += concentrationDot[grid * componentCount + comp];
    const double upwindTotal = -(velocity[grid] * total(grid) - velocity[grid - 1] * total(grid - 1)) / resolution;
    EXPECT_NEAR(sumDot, upwindTotal, 1.0e-10 * std::abs(velocity[0] * total(0) / resolution)) << grid;
  }
}

TEST(MassBalanceWENO3, IsExactForLinearProfilesInGasAndLiquid)
{
  // Linear mole fractions at constant c_tot and v (gas), and linear concentrations at constant v (liquid): the
  // convective term is -v dc/dz at every node except the first-order outlet node.
  constexpr size_t gridPoints = 10;
  constexpr double resolution = 0.02;
  constexpr double v = 0.05;
  const std::vector<Component> components{inertComponent(0, 0.0), inertComponent(1, 0.0)};
  const std::vector<double> distances = uniformDistances(gridPoints, gridPoints * resolution);
  const std::vector<double> fraction = linearCellValues(distances, 0.8, -2.0);
  std::vector<double> concentration(2 * (gridPoints + 1));
  for (size_t grid = 0; grid <= gridPoints; ++grid)
  {
    concentration[2 * grid] = 40.0 * fraction[grid];
    concentration[2 * grid + 1] = 40.0 * (1.0 - fraction[grid]);
  }
  const std::vector<double> velocity(gridPoints + 1, v);
  const std::vector<double> sink(concentration.size(), 0.0);
  const std::vector<double> ones(gridPoints + 1, 1.0);
  for (Column::FluidPhase phase : {Column::FluidPhase::Gas, Column::FluidPhase::Liquid})
  {
    std::vector<double> concentrationDot(concentration.size());
    computeMassDerivatives(singleAdsorbent(components), gridPoints, 2, 1, distances, ones, velocity, concentration,
                           concentrationDot, sink, phase, ones, Column::ConvectionScheme::WENO3);
    for (size_t grid = 1; grid < gridPoints; ++grid)
    {
      EXPECT_NEAR(concentrationDot[2 * grid], -v * 40.0 * -2.0, 1.0e-10) << grid;
      EXPECT_NEAR(concentrationDot[2 * grid + 1], -v * 40.0 * 2.0, 1.0e-10) << grid;
    }
  }
}

TEST(ConvectionSchemeInput, DefaultsToWENO3AndUpwindForSIRK3)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_convection_scheme.json";
  nlohmann::json data = beaInput();
  std::ofstream(path) << data;
  {
    const InputReader reader(path.string());
    EXPECT_EQ(reader.convectionScheme, 1U);
    EXPECT_EQ(Column(reader).convectionScheme, Column::ConvectionScheme::WENO3);
  }

  data["ConvectionScheme"] = "Upwind";
  std::ofstream(path) << data;
  EXPECT_EQ(Column(InputReader(path.string())).convectionScheme, Column::ConvectionScheme::Upwind);

  data["ConvectionScheme"] = "Central";
  std::ofstream(path) << data;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);

  data.erase("ConvectionScheme");
  data["BreakthroughIntegrator"] = "SIRK3";
  std::ofstream(path) << data;
  EXPECT_EQ(InputReader(path.string()).convectionScheme, 0U);

  data["ConvectionScheme"] = "WENO3";
  std::ofstream(path) << data;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);
}
