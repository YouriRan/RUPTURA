#include "component.h"

#include <format>
#include <iostream>
#include <print>
#include <string>
#include <utility>
#include <vector>

#include "isotherm.h"

Component::Component(size_t id, std::string name) : id(id), name(std::move(name)) {}

Component::Component(size_t id, std::string name, std::vector<Isotherm> isotherms, double initialGasMoleFraction,
                     double massTransferCoefficient, double axialDispersionCoefficient, bool isCarrierGas,
                     double molecularWeight, double heatOfAdsorption)
    : id(id),
      name(std::move(name)),
      initialGasMoleFraction(initialGasMoleFraction),
      massTransferCoefficient(massTransferCoefficient),
      axialDispersionCoefficient(axialDispersionCoefficient),
      heatOfAdsorption(heatOfAdsorption),
      isCarrierGas(isCarrierGas),
      molecularWeight(molecularWeight)
{
  for (const Isotherm& site : isotherms)
  {
    isotherm.add(site);
  }
}

Component::Component(size_t id, std::string name, MultiSiteIsotherm isotherm, double initialGasMoleFraction,
                     double massTransferCoefficient, double axialDispersionCoefficient, bool isCarrierGas,
                     double molecularWeight, double heatOfAdsorption)
    : id(id),
      name(std::move(name)),
      isotherm(std::move(isotherm)),
      initialGasMoleFraction(initialGasMoleFraction),
      massTransferCoefficient(massTransferCoefficient),
      axialDispersionCoefficient(axialDispersionCoefficient),
      heatOfAdsorption(heatOfAdsorption),
      isCarrierGas(isCarrierGas),
      molecularWeight(molecularWeight)
{
}

void Component::print() const { std::print("{}", repr()); }

std::string Component::repr() const
{
  std::string s;
  auto out = std::back_inserter(s);

  // The temperature laws of the isotherm parameters are printed with the isotherm sites.
  auto appendHeatOfAdsorption = [&]()
  {
    if (heatOfAdsorption == 0.0) return;

    std::format_to(out, "    heat of adsorption:       {:8.4e} [J/mol]\n", heatOfAdsorption);
  };

  s += "Component id: " + std::to_string(id) + " [" + name + "]:\n";
  if (isCarrierGas)
  {
    s += "    carrier-gas\n";
    appendHeatOfAdsorption();
    s += isotherm.repr();
  }
  s += "    mol-fraction in the gas:   " + std::to_string(initialGasMoleFraction) + " [-]\n";
  if (!isCarrierGas)
  {
    if (isWater) s += "    water (sets the relative humidity)\n";
    s += "    mas-transfer coefficient: " + std::to_string(massTransferCoefficient) + " [1/s]\n";
    s += "    diffusion coefficient:     " + std::to_string(axialDispersionCoefficient) + " [m^2/s]\n";
    appendHeatOfAdsorption();
    s += isotherm.repr();
  }
  if (chemisorption.enabled())
  {
    s += chemisorption.repr();
  }
  return s;
}
