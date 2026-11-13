#pragma once
#include <string>
#include <utility>
#include <vector>

#include "cycle/types.hpp"
namespace cycle
{
std::pair<std::string, std::string> spatial_settings(const Json& cfg);
std::size_t stencil_radius(const Json& cfg);
Vec interface_locations(const Json& cfg);
Vec control_volume_faces(const Vec& z);
// Field-major data[field*nz+node] and output[field*(nz+1)+face].
Vec reconstructed_face_values(const Vec& data, const Vec& z, const Vec& flow, const Vec& left, const Vec& right,
                              const std::string& scheme = "UPWIND", const Vec& interfaces = {},
                              const std::string& method = "FVM");
Vec reconstruct(const Vec& data, std::size_t nfields, std::size_t nz, const Vec& z, const Vec& flow, const Vec& left,
                const Vec& right, const std::string& scheme = "UPWIND", const Vec& interfaces = {});
// Return dc/dt and component face flux, both component-major.
std::pair<Vec, Vec> conservative_species_face_flux_rhs(
    const Vec& concentration, const Vec& face_volume_flux, const Vec& sink, const Vec& fluid_capacity,
    const Vec& dispersion, const Vec& left_reservoir_concentration, const Vec& right_reservoir_concentration,
    const Vec& axial_coordinates, const Vec& face_fluid_capacity = {}, const std::string& scheme = "UPWIND",
    const Vec& interfaces = {}, const Vec& cell_widths = {}, const std::string& dispersion_closure = "concentration",
    const Vec& left_component_flux = {}, const Vec& right_component_flux = {},
    const std::string& method = "FVM");
}  // namespace cycle
