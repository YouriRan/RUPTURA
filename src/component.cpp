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
  isotherm.numberOfSites = isotherms.size();
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

  auto appendThermalScaling = [&]()
  {
    if (!nonIsothermal) return;

    std::format_to(out, "    temperature scaling:      {}\n",
                   referenceTemperature.has_value() ? "exp((beta - beta_ref) H)" : "exp(beta H)");
    std::format_to(out, "        H:       {:8.4e} [J/mol]\n", heatOfAdsorption);
    if (referenceTemperature.has_value())
    {
      std::format_to(out, "        T_ref:   {:8.4e} [K]\n", referenceTemperature.value());
    }
  };

  s += "Component id: " + std::to_string(id) + " [" + name + "]:\n";
  if (isCarrierGas)
  {
    s += "    carrier-gas\n";
    appendThermalScaling();
    s += isotherm.repr();
  }
  s += "    mol-fraction in the gas:   " + std::to_string(initialGasMoleFraction) + " [-]\n";
  if (!isCarrierGas)
  {
    s += "    mas-transfer coefficient: " + std::to_string(massTransferCoefficient) + " [1/s]\n";
    s += "    diffusion coefficient:     " + std::to_string(axialDispersionCoefficient) + " [m^2/s]\n";
    appendThermalScaling();
    s += isotherm.repr();
  }
  return s;
}
