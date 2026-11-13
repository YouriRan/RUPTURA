#pragma once
#include "cycle/types.hpp"
#include "component.h"
#include "geometry.h"
#include "reaction.h"
#include <memory>
#include <optional>

namespace cycle {
// Adapter to the original main material parser and physical rate laws. The
// original JSON names and parameter order are retained at this boundary.
class MainMaterialModels {
 public:
  MainMaterialModels(const Json& input,const std::string& material="");
  ~MainMaterialModels();
  MainMaterialModels(const MainMaterialModels&);
  MainMaterialModels& operator=(const MainMaterialModels&);
  std::vector<Component> components;
  std::vector<Reaction> reactions;
  ::Geometry geometry{makeGeometry(PackedBedTubeSpec{})};
  std::size_t chemical_sites=0;
  bool surface_pore=false;
  double particle_density=1000.;
  Vec equilibrium(const Vec& driving,const Vec& temperature,std::size_t nz) const;
  std::optional<Vec> mpd_maximum_loadings() const;
  Vec chemical_equilibrium(const Vec& concentrations,double gas_temperature,double solid_temperature) const;
  double pH(const Vec& concentrations) const;
 private:
  struct Impl;
  std::unique_ptr<Impl> p;
};
}
