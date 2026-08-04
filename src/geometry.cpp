#include "geometry.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>

namespace
{
constexpr double tinyLength = 1.0e-30;

std::string normalized(std::string_view value)
{
  std::string result{value};
  result.erase(std::remove_if(result.begin(), result.end(),
                              [](unsigned char c) { return std::isspace(c) != 0 || c == '_' || c == '-'; }),
               result.end());
  std::transform(result.begin(), result.end(), result.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return result;
}

void requirePositive(double value, const char* name)
{
  if (value <= 0.0)
  {
    throw std::runtime_error(std::string("Error: Geometry ") + name + " must be positive");
  }
}

void requireNonNegative(double value, const char* name)
{
  if (value < 0.0)
  {
    throw std::runtime_error(std::string("Error: Geometry ") + name + " must be non-negative");
  }
}

struct ChannelShapeTerms
{
  double area;
  double perimeter;
  double hydraulicDiameter;
};

ChannelShapeTerms channelShapeTerms(ChannelShape channelShape, double internalChannelDimension)
{
  requirePositive(internalChannelDimension, "InternalChannelDimension");

  switch (channelShape)
  {
    case ChannelShape::Triangular:
      return {std::sqrt(3.0) * internalChannelDimension * internalChannelDimension / 4.0,
              3.0 * internalChannelDimension, internalChannelDimension / std::sqrt(3.0)};
    case ChannelShape::Square:
      return {internalChannelDimension * internalChannelDimension, 4.0 * internalChannelDimension,
              internalChannelDimension};
    case ChannelShape::Hexagonal:
      return {3.0 * std::sqrt(3.0) * internalChannelDimension * internalChannelDimension / 2.0,
              6.0 * internalChannelDimension, std::sqrt(3.0) * internalChannelDimension};
    case ChannelShape::Circular:
      return {std::numbers::pi * internalChannelDimension * internalChannelDimension / 4.0,
              std::numbers::pi * internalChannelDimension, internalChannelDimension};
  }

  throw std::runtime_error("Error: invalid monolith channel shape");
}
}  // namespace

ChannelShape parseChannelShape(std::string_view value)
{
  const std::string shape = normalized(value);
  if (shape == "tri" || shape == "triangle" || shape == "triangular")
  {
    return ChannelShape::Triangular;
  }
  if (shape == "sq" || shape == "square")
  {
    return ChannelShape::Square;
  }
  if (shape == "hex" || shape == "hexagon" || shape == "hexagonal")
  {
    return ChannelShape::Hexagonal;
  }
  if (shape == "circ" || shape == "circle" || shape == "circular")
  {
    return ChannelShape::Circular;
  }

  throw std::runtime_error("Error: Geometry ChannelShape must be one of triangular, square, hexagonal, circular");
}

std::string_view channelShapeName(ChannelShape channelShape)
{
  switch (channelShape)
  {
    case ChannelShape::Triangular:
      return "triangular";
    case ChannelShape::Square:
      return "square";
    case ChannelShape::Hexagonal:
      return "hexagonal";
    case ChannelShape::Circular:
      return "circular";
  }

  throw std::runtime_error("Error: invalid monolith channel shape");
}

Geometry makeGeometry(const PackedBedTubeSpec& specification)
{
  if (specification.voidFraction <= 0.0 || specification.voidFraction >= 1.0)
  {
    throw std::runtime_error("Error: Geometry ColumnVoidFraction must be between 0 and 1");
  }
  requirePositive(specification.particleDiameter, "ParticleDiameter");
  requireNonNegative(specification.internalDiameter, "InternalDiameter");
  requireNonNegative(specification.outerDiameter, "OuterDiameter");
  if (specification.internalDiameter > 0.0 && specification.outerDiameter > 0.0 &&
      specification.outerDiameter <= specification.internalDiameter)
  {
    throw std::runtime_error("Error: Geometry OuterDiameter must be larger than InternalDiameter");
  }

  const double solidToFluidVolumeRatio = (1.0 - specification.voidFraction) / specification.voidFraction;
  const double solidAreaPerSolidVolume = 6.0 / specification.particleDiameter;
  const double fluidArea = specification.internalDiameter > 0.0 ? std::numbers::pi * specification.internalDiameter *
                                                                      specification.internalDiameter / 4.0
                                                                : 0.0;
  const double outerArea = specification.outerDiameter > 0.0
                               ? std::numbers::pi * specification.outerDiameter * specification.outerDiameter / 4.0
                               : 0.0;
  const double wallArea = outerArea > fluidArea ? outerArea - fluidArea : 0.0;

  Geometry geometry{GeometryKind::HollowTube, specification.voidFraction, solidToFluidVolumeRatio, {}, {}, {}};
  geometry.contactAreas.fluidSolidPerFluidVolume = solidToFluidVolumeRatio * solidAreaPerSolidVolume;
  geometry.contactAreas.solidFluidPerSolidVolume = solidAreaPerSolidVolume;
  geometry.contactAreas.fluidWallPerFluidVolume =
      specification.internalDiameter > 0.0 ? 4.0 / specification.internalDiameter : 0.0;
  geometry.contactAreas.wallInnerPerWallVolume =
      wallArea > 0.0 ? std::numbers::pi * specification.internalDiameter / wallArea : 0.0;
  geometry.contactAreas.wallOuterPerWallVolume =
      wallArea > 0.0 ? std::numbers::pi * specification.outerDiameter / wallArea : 0.0;

  geometry.dimensions.channelShape = ChannelShape::Circular;
  geometry.dimensions.numberOfChannels = 1;
  geometry.dimensions.hydraulicDiameter = specification.internalDiameter;
  geometry.dimensions.solidCharacteristicLength = specification.particleDiameter;
  geometry.dimensions.poreDiffusionLength = 0.5 * specification.particleDiameter;
  geometry.dimensions.internalDiameter = specification.internalDiameter;
  geometry.dimensions.outerDiameter = specification.outerDiameter;
  geometry.dimensions.channelArea = fluidArea;
  geometry.dimensions.channelPerimeter =
      specification.internalDiameter > 0.0 ? std::numbers::pi * specification.internalDiameter : 0.0;
  geometry.dimensions.outerArea = outerArea;
  geometry.dimensions.openArea = fluidArea;
  geometry.dimensions.solidArea = wallArea;

  // Packed-bed pressure-drop coefficients.
  geometry.pressureDrop.viscousCoefficient = 150.0 * solidToFluidVolumeRatio * solidToFluidVolumeRatio;
  geometry.pressureDrop.inertialCoefficient = 1.75 * solidToFluidVolumeRatio;
  return geometry;
}

Geometry makeGeometry(const MonolithSpec& specification)
{
  requirePositive(specification.internalChannelDimension, "InternalChannelDimension");
  requirePositive(specification.outerDiameter, "OuterDiameter");
  if (specification.numberOfChannels == 0)
  {
    throw std::runtime_error("Error: Geometry NumberOfChannels must be positive");
  }
  requireNonNegative(specification.washcoatThickness, "WashcoatThickness");
  if (specification.washcoatVolumePerChannelVolume.has_value())
  {
    requireNonNegative(*specification.washcoatVolumePerChannelVolume, "WashcoatVolumePerChannelVolume");
  }

  const ChannelShapeTerms channel =
      channelShapeTerms(specification.channelShape, specification.internalChannelDimension);
  const double outerArea = std::numbers::pi * specification.outerDiameter * specification.outerDiameter / 4.0;
  const double openArea = static_cast<double>(specification.numberOfChannels) * channel.area;
  if (openArea >= outerArea)
  {
    throw std::runtime_error(
        "Error: invalid Monolith geometry; total channel open area must be smaller than outer area");
  }

  const double solidArea = outerArea - openArea;
  const double voidFraction = openArea / outerArea;
  const double fullSolidToFluidRatio = solidArea / openArea;
  const double fluidSolidAreaPerFluidVolume = channel.perimeter / channel.area;
  const double defaultWashcoatVolumeRatio = fluidSolidAreaPerFluidVolume * specification.washcoatThickness;
  const double activeWashcoatVolumeRatio =
      specification.washcoatVolumePerChannelVolume.value_or(defaultWashcoatVolumeRatio);
  const double solidToFluidVolumeRatio =
      activeWashcoatVolumeRatio > 0.0 ? activeWashcoatVolumeRatio : fullSolidToFluidRatio;
  const double solidFluidAreaPerSolidVolume =
      fluidSolidAreaPerFluidVolume / std::max(solidToFluidVolumeRatio, tinyLength);
  const double equivalentSolidLength = 6.0 / std::max(solidFluidAreaPerSolidVolume, tinyLength);

  Geometry geometry{GeometryKind::Monolith, voidFraction, solidToFluidVolumeRatio, {}, {}, {}};
  geometry.contactAreas.fluidSolidPerFluidVolume = fluidSolidAreaPerFluidVolume;
  geometry.contactAreas.solidFluidPerSolidVolume = solidFluidAreaPerSolidVolume;
  geometry.contactAreas.fluidWallPerFluidVolume = fluidSolidAreaPerFluidVolume;
  geometry.contactAreas.wallInnerPerWallVolume =
      static_cast<double>(specification.numberOfChannels) * channel.perimeter / solidArea;
  geometry.contactAreas.wallOuterPerWallVolume = std::numbers::pi * specification.outerDiameter / solidArea;

  geometry.dimensions.channelShape = specification.channelShape;
  geometry.dimensions.numberOfChannels = specification.numberOfChannels;
  geometry.dimensions.hydraulicDiameter = channel.hydraulicDiameter;
  geometry.dimensions.solidCharacteristicLength = equivalentSolidLength;
  geometry.dimensions.poreDiffusionLength =
      specification.washcoatThickness > 0.0 ? specification.washcoatThickness : 0.5 * equivalentSolidLength;
  geometry.dimensions.internalDiameter = specification.internalChannelDimension;
  geometry.dimensions.outerDiameter = specification.outerDiameter;
  geometry.dimensions.channelArea = channel.area;
  geometry.dimensions.channelPerimeter = channel.perimeter;
  geometry.dimensions.outerArea = outerArea;
  geometry.dimensions.openArea = openArea;
  geometry.dimensions.solidArea = solidArea;
  geometry.dimensions.washcoatThickness = specification.washcoatThickness;
  geometry.dimensions.washcoatVolumePerChannelVolume = activeWashcoatVolumeRatio;

  // Hydraulic-diameter laminar-channel approximation.
  geometry.pressureDrop.viscousCoefficient =
      32.0 / std::max(channel.hydraulicDiameter * channel.hydraulicDiameter, tinyLength);
  geometry.pressureDrop.inertialCoefficient = 0.0;
  return geometry;
}
