#include "jacobian.h"

#include <algorithm>
#include <array>
#include <utility>

#include "compute.h"
#include "sorption.h"
#include "utils.h"

SorptionJacobian::SorptionJacobian(const Column& column)
    : numberOfComponents(column.numberOfComponents),
      numberOfNodes(column.numberOfGridPoints + 1),
      size(2 * column.numberOfComponents * (column.numberOfGridPoints + 1)),
      // A WENO3 row reaches the concentrations two nodes upstream, up to 5N - 1 positions before the diagonal.
      bandwidth((column.convectionScheme == Column::ConvectionScheme::WENO3 ? 5 : 2) * column.numberOfComponents),
      leadingDimension(3 * bandwidth + 1),
      band(leadingDimension * size, 0.0),
      equilibriumSensitivity((column.numberOfGridPoints + 1) * column.numberOfComponents * column.numberOfComponents),
      mixtureSensitivity(column.numberOfComponents * column.numberOfComponents)
{
}

void SorptionJacobian::gather(std::span<const double> concentration, std::span<const double> loading,
                              std::span<double> nodeMajor) const noexcept
{
  for (size_t node = 0; node < numberOfNodes; ++node)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      nodeMajor[concentrationIndex(node, comp)] = concentration[node * numberOfComponents + comp];
      nodeMajor[loadingIndex(node, comp)] = loading[node * numberOfComponents + comp];
    }
  }
}

void SorptionJacobian::scatter(std::span<const double> nodeMajor, std::span<double> concentration,
                               std::span<double> loading) const noexcept
{
  for (size_t node = 0; node < numberOfNodes; ++node)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      concentration[node * numberOfComponents + comp] = nodeMajor[concentrationIndex(node, comp)];
      loading[node * numberOfComponents + comp] = nodeMajor[loadingIndex(node, comp)];
    }
  }
}

namespace
{
/**
 * \brief Adds the WENO3 convection of computeMassDerivatives to the concentration rows: the derivatives of the face
 *        fluxes with respect to the concentrations of the nodes i - 1, i, and i + 1 of each face i with v_i >= 0, and
 *        of the nodes i, i + 1, and i + 2 of each face with v_i < 0.
 *
 * In the gas phase F_ik = eps_i v_i C_u f_ik / S_i, with C_u the total concentration of the node upstream of the face,
 * f_ik the face value of the mole fraction y_k = c_k / C, and S_i = sum_k f_ik, so that every component couples to
 * every other; in the liquid phase F_ik = eps_i v_i f_ik with f the face value of c_k. Node 0 holds the inlet-face
 * value, and node N + 1 stands for the value at the face z_N for inflow there (outletFace), which changes with c_N by
 * outletSlope.
 */
void addWenoConvection(const Column& column, std::span<const double> concentration, std::span<const double> outletFace,
                       std::span<const double> outletSlope, SorptionJacobian& jacobian)
{
  const size_t n = column.numberOfComponents;
  const size_t nodes = column.numberOfGridPoints + 1;
  const size_t outlet = column.numberOfGridPoints;
  const bool gas = column.fluidPhase == Column::FluidPhase::Gas;

  // The node concentrations with node N + 1 appended, the reconstructed node values u (mole fractions or
  // concentrations), their face values, and the slopes d f_ik / d u_jk with respect to the three source nodes of the
  // face, in the order of their index.
  std::vector<double> extended(concentration.begin(), concentration.begin() + static_cast<std::ptrdiff_t>(nodes * n));
  extended.insert(extended.end(), outletFace.begin(), outletFace.end());
  std::vector<double> total(nodes + 1, 0.0);
  for (size_t node = 0; node <= nodes; ++node)
  {
    for (size_t comp = 0; comp < n; ++comp) total[node] += std::max(0.0, extended[node * n + comp]);
  }
  std::vector<double> nodeValues((nodes + 1) * n);
  for (size_t node = 0; node <= nodes; ++node)
  {
    for (size_t comp = 0; comp < n; ++comp)
    {
      const double value = extended[node * n + comp];
      nodeValues[node * n + comp] = gas ? (total[node] > 0.0 ? std::max(0.0, value) / total[node] : 0.0) : value;
    }
  }
  const bool anyReverseFace = std::any_of(column.interstitialGasVelocity.begin(), column.interstitialGasVelocity.end(),
                                          [](double velocity) { return velocity < 0.0; });
  std::vector<double> faceValues(nodes * n);
  std::vector<double> faceSlopes(3 * nodes * n);
  std::vector<double> values(nodes), faces(nodes), lower(nodes), centre(nodes), upper(nodes);
  std::vector<double> reverseFaces(nodes), downstream(nodes), reverseCentre(nodes), upstream(nodes);
  for (size_t comp = 0; comp < n; ++comp)
  {
    for (size_t node = 0; node < nodes; ++node) values[node] = nodeValues[node * n + comp];
    computeWENO(values, column.columnDistances, faces);
    computeWENOLinearization(values, column.columnDistances, lower, centre, upper);
    if (anyReverseFace)
    {
      const double inflow = nodeValues[nodes * n + comp];
      computeWENOReverse(values, inflow, column.columnDistances, reverseFaces);
      computeWENOReverseLinearization(values, inflow, column.columnDistances, downstream, reverseCentre, upstream);
    }
    for (size_t node = 0; node < nodes; ++node)
    {
      const bool reverseFace = column.interstitialGasVelocity[node] < 0.0;
      faceValues[node * n + comp] = reverseFace ? reverseFaces[node] : faces[node];
      faceSlopes[3 * (node * n + comp)] = reverseFace ? downstream[node] : lower[node];
      faceSlopes[3 * (node * n + comp) + 1] = reverseFace ? reverseCentre[node] : centre[node];
      faceSlopes[3 * (node * n + comp) + 2] = reverseFace ? upstream[node] : upper[node];
    }
  }

  std::vector<double> block(n * n);
  for (size_t face = 0; face < nodes; ++face)
  {
    double faceSum = 0.0;
    for (size_t comp = 0; comp < n; ++comp) faceSum += faceValues[face * n + comp];
    const double velocity = column.interstitialGasVelocity[face];
    const double superficialVelocity = column.totalVoidFraction[face] * velocity;
    // The source nodes are face - 1 ... face + 1 (slots 0 ... 2) for v >= 0, face ... face + 2 for v < 0; the node
    // upstream of the face is the centre one.
    const size_t first = velocity < 0.0 ? face : face - 1;  // wraps for face 0, whose slot 0 is skipped
    const size_t upwind = velocity < 0.0 ? face + 1 : face;
    for (size_t slot = 0; slot < 3; ++slot)
    {
      const size_t source = first + slot;
      if ((velocity >= 0.0 && face == 0 && slot == 0) || source > nodes) continue;
      if (velocity >= 0.0 && source >= nodes) continue;
      // block = dF_face,k / dc_source,m
      std::fill(block.begin(), block.end(), 0.0);
      bool nonzero = false;
      if (gas)
      {
        if (faceSum <= 0.0 || total[source] <= 0.0) continue;
        const double scale = superficialVelocity * total[upwind] / (faceSum * total[source]);
        double weightedSlope = 0.0;
        for (size_t comp = 0; comp < n; ++comp)
        {
          weightedSlope += faceSlopes[3 * (face * n + comp) + slot] * nodeValues[source * n + comp];
        }
        for (size_t k = 0; k < n; ++k)
        {
          const double normalized = faceValues[face * n + k] / faceSum;
          const double slope = faceSlopes[3 * (face * n + k) + slot];
          for (size_t m = 0; m < n; ++m)
          {
            if (extended[source * n + m] <= 0.0) continue;  // the clamp max(0, c)
            const double fractionSlope = slope * ((k == m ? 1.0 : 0.0) - nodeValues[source * n + k]);
            const double sumSlope = faceSlopes[3 * (face * n + m) + slot] - weightedSlope;
            block[k + m * n] =
                (source == upwind ? superficialVelocity * normalized : 0.0) +
                scale * (fractionSlope - normalized * sumSlope);
            nonzero = true;
          }
        }
      }
      else
      {
        for (size_t k = 0; k < n; ++k)
        {
          block[k + k * n] = superficialVelocity * faceSlopes[3 * (face * n + k) + slot];
          nonzero = true;
        }
      }
      if (!nonzero) continue;

      // The flux through face i leaves cell i and enters cell i + 1; node N + 1 is c_N through outletSlope.
      const size_t target = std::min(source, outlet);
      for (const auto& [node, sign] : {std::pair{face, -1.0}, std::pair{face + 1, 1.0}})
      {
        if (node == 0 || node >= nodes) continue;
        const double factor =
            sign / std::max(1e-30, gridSpacing(column.columnDistances, node) * column.totalVoidFraction[node]);
        for (size_t k = 0; k < n; ++k)
        {
          for (size_t m = 0; m < n; ++m)
          {
            if (block[k + m * n] == 0.0) continue;
            const double chain = source == nodes ? outletSlope[m] : 1.0;
            jacobian(jacobian.concentrationIndex(node, k), jacobian.concentrationIndex(target, m)) +=
                factor * block[k + m * n] * chain;
          }
        }
      }
    }
  }
}
}  // namespace

void assembleSorptionJacobian(Column& column, std::span<const double> concentration,
                              std::span<const double> gasTemperature, SorptionJacobian& jacobian)
{
  const size_t n = column.numberOfComponents;
  const size_t outlet = column.numberOfGridPoints;
  const bool gas = column.fluidPhase == Column::FluidPhase::Gas;
  const bool weno = column.convectionScheme == Column::ConvectionScheme::WENO3;
  const bool reverse = column.flowDirection == Column::FlowDirection::Reverse;
  const bool hasSensitivities =
      std::none_of(column.physisorptionMixtures.begin(), column.physisorptionMixtures.end(),
                   [](const MixturePrediction& mixture)
                   { return mixture.predictionMethod == MixturePrediction::PredictionMethod::MPD; });
  if (hasSensitivities)
  {
    computePhysisorptionEquilibriumSensitivities(
        column.physisorptionMixtures, column.numberOfGridPoints, n, column.numberOfAdsorbents,
        column.fractionOfAdsorbent, column.hasAdsorbentOfType, column.maxIsothermTerms, column.idealGasMolFractions,
        jacobian.mixtureSensitivity, column.totalPressure, column.cachedPressure, column.moleFraction, gasTemperature,
        gas ? MixturePrediction::DrivingForceInput::MoleFraction : MixturePrediction::DrivingForceInput::Concentration,
        concentration, column.pH, jacobian.equilibriumSensitivity);
  }
  else
  {
    std::fill(jacobian.equilibriumSensitivity.begin(), jacobian.equilibriumSensitivity.end(), 0.0);
  }
  std::fill(jacobian.band.begin(), jacobian.band.end(), 0.0);

  // Danckwerts inlet: computeMassDerivatives transports the inlet-face value c_b = (u c_in + g c_1) / (u + g) in place
  // of the feed c_in at node 0, with g = D over half the first cell (WENO3) or the first dispersion distance (upwind).
  // In a reverse step node 0 is replaced by the first cell. inletWeights holds u, g, and u + g, with dc_b/dc_in =
  // u / (u + g) and dc_b/dc_1 = g / (u + g).
  std::vector<double> transported(concentration.begin(), concentration.end());
  std::vector<std::array<double, 3>> inletWeights(n, {1.0, 0.0, 1.0});
  const bool danckwerts = column.inletCondition == Column::InletCondition::Danckwerts;
  const double inflowVelocity =
      std::max(0.0, column.totalVoidFraction[0] * column.interstitialGasVelocity[0] / column.totalVoidFraction[1]);
  if (reverse)
  {
    for (size_t comp = 0; comp < n; ++comp)
    {
      transported[comp] = concentration[n + comp];
      inletWeights[comp] = {0.0, 1.0, 1.0};
    }
  }
  else if (danckwerts)
  {
    const double firstSpacing = gridSpacing(column.columnDistances, 1);
    const double inletDistance = weno         ? 0.5 * firstSpacing
                                 : outlet > 1 ? 0.5 * (firstSpacing + gridSpacing(column.columnDistances, 2))
                                              : firstSpacing;
    for (size_t comp = 0; comp < n; ++comp)
    {
      const double inletConductance = column.axialDispersionCoefficient[n + comp] / std::max(1e-30, inletDistance);
      const double sum = inflowVelocity + inletConductance;
      if (sum > 0.0)
      {
        transported[comp] = (inflowVelocity * concentration[comp] + inletConductance * concentration[n + comp]) / sum;
        inletWeights[comp] = {inflowVelocity, inletConductance, sum};
      }
    }
  }

  // Gas that flows in at z = L, as in computeMassDerivatives: the face value c_b of the product-end gas (or of the last
  // cell), with outletSlope = dc_b/dc_N, and the dispersion across that face.
  const double lastSpacing = gridSpacing(column.columnDistances, outlet);
  const double outletDistance = weno         ? 0.5 * lastSpacing
                                : outlet > 1 ? 0.5 * (lastSpacing + gridSpacing(column.columnDistances, outlet - 1))
                                             : lastSpacing;
  const double outletInflowVelocity = std::max(0.0, -column.interstitialGasVelocity[outlet]);
  const bool outletDispersion = danckwerts ? outletInflowVelocity > 0.0 : reverse;
  const bool followsLastCell = column.productEndConcentration.empty();
  std::vector<double> outletFace(transported.begin() + static_cast<std::ptrdiff_t>(outlet * n), transported.end());
  std::vector<double> outletSlope(n, 1.0);
  if (outletInflowVelocity > 0.0 || reverse)
  {
    for (size_t comp = 0; comp < n; ++comp)
    {
      const double last = transported[outlet * n + comp];
      const double reservoir = followsLastCell ? last : column.productEndConcentration[comp];
      if (!danckwerts)
      {
        outletFace[comp] = reservoir;
        outletSlope[comp] = followsLastCell ? 1.0 : 0.0;
        continue;
      }
      if (outletInflowVelocity <= 0.0) continue;
      const double conductance = column.axialDispersionCoefficient[outlet * n + comp] / std::max(1e-30, outletDistance);
      const double sum = outletInflowVelocity + conductance;
      outletFace[comp] = (outletInflowVelocity * reservoir + conductance * last) / sum;
      outletSlope[comp] = followsLastCell ? 1.0 : conductance / sum;
    }
  }

  for (size_t node = 0; node <= outlet; ++node)
  {
    const std::span<const double> sensitivity(jacobian.equilibriumSensitivity.data() + node * n * n, n * n);
    for (size_t comp = 0; comp < n; ++comp)
    {
      // LDF rate r = k (q* - q); the concentration rows lose loadingPrefactor r.
      const size_t c = jacobian.concentrationIndex(node, comp);
      const size_t q = jacobian.loadingIndex(node, comp);
      const double rate = column.massTransferCoefficient[node * n + comp];
      jacobian(q, q) = -rate;
      for (size_t other = 0; other < n; ++other)
      {
        jacobian(q, jacobian.concentrationIndex(node, other)) = rate * sensitivity[comp + other * n];
      }
      if (node == 0)
      {
        // The inlet concentrations are prescribed; those of a forward step follow a pressure ramp.
        if (!reverse) jacobian(c, c) = column.relativePressureRate;
        continue;
      }

      const double prefactor = column.loadingPrefactor[node];
      jacobian(c, q) = prefactor * rate;
      for (size_t other = 0; other < n; ++other)
      {
        jacobian(c, jacobian.concentrationIndex(node, other)) = -prefactor * rate * sensitivity[comp + other * n];
      }

      // Upwind convection with the velocity of the last precompute; fluxes are superficial. Each face takes the node
      // upstream of it, the value at z = L for inflow there. WENO3 convection follows below.
      const double hm = gridSpacing(column.columnDistances, node);
      const size_t upstream = jacobian.concentrationIndex(node - 1, comp);
      if (!weno)
      {
        const double denominator = std::max(1e-30, hm * column.totalVoidFraction[node]);
        const double lowerFlux =
            column.totalVoidFraction[node - 1] * column.interstitialGasVelocity[node - 1] / denominator;
        const double upperFlux = column.totalVoidFraction[node] * column.interstitialGasVelocity[node] / denominator;
        jacobian(c, column.interstitialGasVelocity[node - 1] >= 0.0 ? upstream : c) += lowerFlux;
        if (column.interstitialGasVelocity[node] >= 0.0)
        {
          jacobian(c, c) -= upperFlux;
        }
        else if (node < outlet)
        {
          jacobian(c, jacobian.concentrationIndex(node + 1, comp)) -= upperFlux;
        }
        else
        {
          jacobian(c, c) -= upperFlux * outletSlope[comp];
        }
      }

      // Axial dispersion: three-point second difference, zero gradient at the outlet; with WENO3 the inlet value lies
      // on the inlet face, half a cell from the centre of the first cell.
      const double dispersion = column.axialDispersionCoefficient[node * n + comp];
      if (node == outlet)
      {
        const double inverseSpacing = 1.0 / std::max(1e-30, hm);
        jacobian(c, upstream) += dispersion * inverseSpacing * inverseSpacing;
        jacobian(c, c) -= dispersion * inverseSpacing * inverseSpacing;
        if (outletDispersion)
        {
          jacobian(c, c) += dispersion * (outletSlope[comp] - 1.0) / std::max(1e-30, lastSpacing * outletDistance);
        }
      }
      else if (weno && node == 1)
      {
        const double hp = gridSpacing(column.columnDistances, node + 1);
        const double downstreamCoefficient = dispersion / std::max(1e-30, 0.5 * (hm + hp) * hm);
        const double upstreamCoefficient = dispersion / std::max(1e-30, 0.5 * hm * hm);
        jacobian(c, jacobian.concentrationIndex(node + 1, comp)) += downstreamCoefficient;
        jacobian(c, c) -= downstreamCoefficient + upstreamCoefficient;
        jacobian(c, upstream) += upstreamCoefficient;
      }
      else
      {
        const double hp = gridSpacing(column.columnDistances, node + 1);
        const double denominator = std::max(1e-30, hm * hp * (hm + hp));
        jacobian(c, upstream) += 2.0 * dispersion * hp / denominator;
        jacobian(c, c) -= 2.0 * dispersion * (hm + hp) / denominator;
        jacobian(c, jacobian.concentrationIndex(node + 1, comp)) += 2.0 * dispersion * hm / denominator;
      }
    }
  }
  if (weno) addWenoConvection(column, transported, outletFace, outletSlope, jacobian);

  // The entries for node 0 are those of its stencil value: c_b of a Danckwerts inlet, with dc_b/dc_in = u / (u + g) and
  // dc_b/dc_1 = g / (u + g), or the first cell in a reverse step.
  if (danckwerts || reverse)
  {
    for (size_t node = 1; node <= std::min<size_t>(2, outlet); ++node)
    {
      for (size_t row = 0; row < n; ++row)
      {
        const size_t c = jacobian.concentrationIndex(node, row);
        for (size_t comp = 0; comp < n; ++comp)
        {
          const size_t inlet = jacobian.concentrationIndex(0, comp);
          const auto& [first, second, sum] = inletWeights[comp];
          if (jacobian.entry(c, inlet) == 0.0) continue;
          const double inletCoefficient = jacobian(c, inlet);
          jacobian(c, inlet) = inletCoefficient * first / sum;
          jacobian(c, jacobian.concentrationIndex(1, comp)) += inletCoefficient * second / sum;
        }
      }
    }
  }
}

void assembleThermalJacobian(const Column& column, BlockBandedMatrix& jacobian)
{
  // Node blocks (T_g, T_s, T_w); see computeEnergyDerivatives for the coefficients.
  constexpr size_t gasT = 0;
  constexpr size_t solidT = 1;
  constexpr size_t wallT = 2;
  jacobian.setZero();
  const size_t outlet = column.numberOfGridPoints;
  const auto& areas = column.geometries.front().contactAreas;
  const double gasWallPrefactor = areas.fluidWallPerFluidVolume * column.heatTransferGasWall / column.heatCapacityGas;
  const double inverseWallHeatDensity = 1.0 / (column.heatCapacityWall * column.wallDensity);
  const double wallGas = column.heatTransferGasWall * areas.wallInnerPerWallVolume * inverseWallHeatDensity;
  const double wallExternal = column.heatTransferWallExternal * areas.wallOuterPerWallVolume * inverseWallHeatDensity;
  const double wallDiffusivity = column.wallThermalConductivity * inverseWallHeatDensity;
  const bool danckwerts = column.inletCondition == Column::InletCondition::Danckwerts;
  const bool reverse = column.flowDirection == Column::FlowDirection::Reverse;

  // Gas that flows in at z = L, as in GasTemperatureBalance::rate: the face temperature beyond the last node changes
  // with T_N by outletSlope, and the conduction across that face (conductsAtOutlet) over outletDistance.
  const bool followsLastCell = !(column.productEndTemperature > 0.0);
  const double lastSpacing = gridSpacing(column.columnDistances, outlet);
  const double outletDistance =
      outlet > 1 ? 0.5 * (lastSpacing + gridSpacing(column.columnDistances, outlet - 1)) : lastSpacing;
  const double outletInflowVelocity =
      outlet > 0 ? std::max(0.0, -column.interstitialGasVelocity[outlet - 1] * column.totalVoidFraction[outlet - 1] /
                                     column.totalVoidFraction[outlet])
                 : 0.0;
  double outletSlope = 1.0;
  bool conductsAtOutlet = false;
  if (danckwerts && outletInflowVelocity > 0.0)
  {
    const double conductance = column.coeffDiffusion[outlet] / std::max(1e-30, outletDistance);
    outletSlope = followsLastCell ? 1.0 : conductance / (outletInflowVelocity + conductance);
    conductsAtOutlet = true;
  }
  else if (!danckwerts)
  {
    outletSlope = followsLastCell ? 1.0 : 0.0;
    conductsAtOutlet = reverse;
  }

  for (size_t node = 0; node <= outlet; ++node)
  {
    const double solidGas = column.solidFluidContactArea[node] * column.heatTransferGasSolid /
                            (column.heatCapacitySolid * std::max(1e-10, column.particleDensity[node]));
    jacobian(node, node, solidT, solidT) = -solidGas;
    jacobian(node, node, solidT, gasT) = solidGas;
    jacobian(node, node, wallT, wallT) = -wallGas - wallExternal;
    jacobian(node, node, wallT, gasT) = wallGas;
    if (node == 0) continue;  // the inlet gas temperature is prescribed, and the inlet wall node does not conduct

    const double inverseGasDensity = 1.0 / std::max(1e-10, column.gasDensity[node]);
    const double gasSolid = column.fluidSolidContactArea[node] * column.heatTransferGasSolid * inverseGasDensity /
                            column.heatCapacityGas;
    const double gasWall = gasWallPrefactor * inverseGasDensity;
    jacobian(node, node, gasT, gasT) = -gasSolid - gasWall;
    jacobian(node, node, gasT, solidT) = gasSolid;
    jacobian(node, node, gasT, wallT) = gasWall;

    // Convection -max(v_i, 0) (T_i - T_{i-1}) / h_i + min(w_{i-1}, 0) (T_i - T_{i+1}) / h_i, with
    // w_{i-1} = eps_{i-1} v_{i-1} / eps_i the velocity of node i at its lower face.
    const double hm = gridSpacing(column.columnDistances, node);
    const double convection = std::max(0.0, column.interstitialGasVelocity[node]) / std::max(1e-30, hm);
    jacobian(node, node, gasT, gasT) -= convection;
    jacobian(node, node - 1, gasT, gasT) += convection;
    if (column.interstitialGasVelocity[node - 1] < 0.0)
    {
      const double backward = column.interstitialGasVelocity[node - 1] * column.totalVoidFraction[node - 1] /
                              (column.totalVoidFraction[node] * std::max(1e-30, hm));
      jacobian(node, node, gasT, gasT) += backward;
      if (node < outlet)
      {
        jacobian(node, node + 1, gasT, gasT) -= backward;
      }
      else
      {
        jacobian(node, node, gasT, gasT) -= backward * outletSlope;
      }
    }

    // Conduction of the gas and of the wall: three-point second difference, zero gradient at the outlet.
    for (const auto& [temperature, diffusivity] :
         {std::pair{gasT, column.coeffDiffusion[node]}, std::pair{wallT, wallDiffusivity}})
    {
      if (node == outlet)
      {
        const double coefficient = diffusivity / std::max(1e-30, hm * hm);
        jacobian(node, node - 1, temperature, temperature) += coefficient;
        jacobian(node, node, temperature, temperature) -= coefficient;
        if (temperature == gasT && conductsAtOutlet)
        {
          jacobian(node, node, gasT, gasT) +=
              diffusivity * (outletSlope - 1.0) / std::max(1e-30, hm * outletDistance);
        }
        continue;
      }
      const double hp = gridSpacing(column.columnDistances, node + 1);
      const double denominator = std::max(1e-30, hm * hp * (hm + hp));
      jacobian(node, node - 1, temperature, temperature) += 2.0 * diffusivity * hp / denominator;
      jacobian(node, node, temperature, temperature) -= 2.0 * diffusivity * (hm + hp) / denominator;
      jacobian(node, node + 1, temperature, temperature) += 2.0 * diffusivity * hm / denominator;
    }
  }

  // In a reverse step node 1 takes its own temperature in place of T_0 (zero gradient at z = 0).
  if (reverse && outlet > 0)
  {
    jacobian(1, 1, gasT, gasT) += jacobian(1, 0, gasT, gasT);
    jacobian(1, 0, gasT, gasT) = 0.0;
  }
  // Danckwerts inlet: node 1 transports T_b = (v_1 T_in + g T_1) / (v_1 + g) in place of T_in.
  else if (danckwerts && outlet > 0)
  {
    const double firstSpacing = gridSpacing(column.columnDistances, 1);
    const double inletDistance =
        outlet > 1 ? 0.5 * (firstSpacing + gridSpacing(column.columnDistances, 2)) : firstSpacing;
    const double velocity = std::max(0.0, column.interstitialGasVelocity[1]);
    const double conductance = column.coeffDiffusion[1] / std::max(1e-30, inletDistance);
    if (velocity + conductance > 0.0)
    {
      const double inletCoefficient = jacobian(1, 0, gasT, gasT);
      jacobian(1, 0, gasT, gasT) = inletCoefficient * velocity / (velocity + conductance);
      jacobian(1, 1, gasT, gasT) += inletCoefficient * conductance / (velocity + conductance);
    }
  }
}
