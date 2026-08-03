#include "geometry.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <numbers>
#include <stdexcept>

namespace
{
constexpr double tinyLength = 1.0e-30;

std::string normalized(std::string value)
{
  value.erase(std::remove_if(value.begin(), value.end(),
                             [](unsigned char c) { return std::isspace(c) != 0 || c == '_' || c == '-'; }),
              value.end());
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
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

ChannelShapeTerms channelShapeTerms(Monolith::ChannelShape channelShape, double internalChannelDimension)
{
  requirePositive(internalChannelDimension, "InternalChannelDimension");

  switch (channelShape)
  {
    case Monolith::ChannelShape::Triangular:
      return {std::sqrt(3.0) * internalChannelDimension * internalChannelDimension / 4.0,
              3.0 * internalChannelDimension,
              internalChannelDimension / std::sqrt(3.0)};
    case Monolith::ChannelShape::Square:
      return {internalChannelDimension * internalChannelDimension,
              4.0 * internalChannelDimension,
              internalChannelDimension};
    case Monolith::ChannelShape::Hexagonal:
      return {3.0 * std::sqrt(3.0) * internalChannelDimension * internalChannelDimension / 2.0,
              6.0 * internalChannelDimension,
              std::sqrt(3.0) * internalChannelDimension};
    case Monolith::ChannelShape::Circular:
      return {std::numbers::pi * internalChannelDimension * internalChannelDimension / 4.0,
              std::numbers::pi * internalChannelDimension,
              internalChannelDimension};
  }

  throw std::runtime_error("Error: invalid monolith channel shape");
}
}  // namespace

HollowTube::HollowTube(double voidFraction, double particleDiameter,
                       double internalDiameter, double outerDiameter)
{
  if (voidFraction <= 0.0 || voidFraction >= 1.0)
  {
    throw std::runtime_error("Error: Geometry ColumnVoidFraction must be between 0 and 1");
  }
  requirePositive(particleDiameter, "ParticleDiameter");
  requireNonNegative(internalDiameter, "InternalDiameter");
  requireNonNegative(outerDiameter, "OuterDiameter");
  if (internalDiameter > 0.0 && outerDiameter > 0.0 && outerDiameter <= internalDiameter)
  {
    throw std::runtime_error("Error: Geometry OuterDiameter must be larger than InternalDiameter");
  }

  const double solidToFluidVolumeRatio = (1.0 - voidFraction) / voidFraction;
  const double solidAreaPerSolidVolume = 6.0 / particleDiameter;
  const double fluidArea = internalDiameter > 0.0
                               ? std::numbers::pi * internalDiameter * internalDiameter / 4.0
                               : 0.0;
  const double outerArea = outerDiameter > 0.0
                               ? std::numbers::pi * outerDiameter * outerDiameter / 4.0
                               : 0.0;
  const double wallArea = outerArea > fluidArea ? outerArea - fluidArea : 0.0;

  parameters.geometryKind = ShapeParameters::GeometryKind::HollowTube;
  parameters.geometryName = "HollowTube";
  parameters.channelShape = "circular";
  parameters.numberOfChannels = 1;
  parameters.voidFraction = voidFraction;
  parameters.solidToFluidVolumeRatio = solidToFluidVolumeRatio;
  parameters.fluidSolidContactAreaPerFluidVolume = solidToFluidVolumeRatio * solidAreaPerSolidVolume;
  parameters.solidFluidContactAreaPerSolidVolume = solidAreaPerSolidVolume;
  parameters.fluidWallContactAreaPerFluidVolume = internalDiameter > 0.0 ? 4.0 / internalDiameter : 0.0;
  parameters.wallInnerContactAreaPerWallVolume =
      wallArea > 0.0 ? std::numbers::pi * internalDiameter / wallArea : 0.0;
  parameters.wallOuterContactAreaPerWallVolume =
      wallArea > 0.0 ? std::numbers::pi * outerDiameter / wallArea : 0.0;
  parameters.hydraulicDiameter = internalDiameter;
  parameters.solidCharacteristicLength = particleDiameter;
  parameters.poreDiffusionLength = 0.5 * particleDiameter;
  parameters.internalDiameter = internalDiameter;
  parameters.outerDiameter = outerDiameter;
  parameters.channelArea = fluidArea;
  parameters.channelPerimeter = internalDiameter > 0.0 ? std::numbers::pi * internalDiameter : 0.0;
  parameters.outerArea = outerArea;
  parameters.openArea = fluidArea;
  parameters.solidArea = wallArea;

  // Preserve the legacy packed-bed pressure-drop coefficients.
  parameters.viscousPressureDropCoefficient = 150.0 * solidToFluidVolumeRatio * solidToFluidVolumeRatio;
  parameters.inertialPressureDropCoefficient = 1.75 * solidToFluidVolumeRatio;
}

Monolith::Monolith(ChannelShape channelShape, double internalChannelDimension,
                   double outerDiameter, std::size_t numberOfChannels,
                   double washcoatThickness, double washcoatVolumePerChannelVolume)
{
  requirePositive(internalChannelDimension, "InternalChannelDimension");
  requirePositive(outerDiameter, "OuterDiameter");
  if (numberOfChannels == 0)
  {
    throw std::runtime_error("Error: Geometry NumberOfChannels must be positive");
  }
  requireNonNegative(washcoatThickness, "WashcoatThickness");
  if (washcoatVolumePerChannelVolume < 0.0 && washcoatVolumePerChannelVolume != -1.0)
  {
    throw std::runtime_error("Error: Geometry WashcoatVolumePerChannelVolume must be non-negative");
  }

  const ChannelShapeTerms channel = channelShapeTerms(channelShape, internalChannelDimension);
  const double outerArea = std::numbers::pi * outerDiameter * outerDiameter / 4.0;
  const double openArea = static_cast<double>(numberOfChannels) * channel.area;
  if (openArea >= outerArea)
  {
    throw std::runtime_error(
        "Error: invalid Monolith geometry; total channel open area must be smaller than outer area");
  }

  const double solidArea = outerArea - openArea;
  const double voidFraction = openArea / outerArea;
  const double fullSolidToFluidRatio = solidArea / openArea;
  const double fluidSolidAreaPerFluidVolume = channel.perimeter / channel.area;
  const double defaultWashcoatVolumeRatio = fluidSolidAreaPerFluidVolume * washcoatThickness;
  const bool hasWashcoatRatioOverride = washcoatVolumePerChannelVolume >= 0.0;
  const double activeWashcoatVolumeRatio = hasWashcoatRatioOverride
                                               ? washcoatVolumePerChannelVolume
                                               : defaultWashcoatVolumeRatio;
  const double solidToFluidVolumeRatio =
      activeWashcoatVolumeRatio > 0.0 ? activeWashcoatVolumeRatio : fullSolidToFluidRatio;
  const double solidFluidAreaPerSolidVolume =
      fluidSolidAreaPerFluidVolume / std::max(solidToFluidVolumeRatio, tinyLength);
  const double equivalentSolidLength = 6.0 / std::max(solidFluidAreaPerSolidVolume, tinyLength);

  parameters.geometryKind = ShapeParameters::GeometryKind::Monolith;
  parameters.geometryName = "Monolith";
  parameters.channelShape = channelShapeName(channelShape);
  parameters.numberOfChannels = numberOfChannels;
  parameters.voidFraction = voidFraction;
  parameters.solidToFluidVolumeRatio = solidToFluidVolumeRatio;
  parameters.fluidSolidContactAreaPerFluidVolume = fluidSolidAreaPerFluidVolume;
  parameters.solidFluidContactAreaPerSolidVolume = solidFluidAreaPerSolidVolume;
  parameters.fluidWallContactAreaPerFluidVolume = fluidSolidAreaPerFluidVolume;
  parameters.wallInnerContactAreaPerWallVolume =
      static_cast<double>(numberOfChannels) * channel.perimeter / solidArea;
  parameters.wallOuterContactAreaPerWallVolume = std::numbers::pi * outerDiameter / solidArea;
  parameters.hydraulicDiameter = channel.hydraulicDiameter;
  parameters.solidCharacteristicLength = equivalentSolidLength;
  parameters.poreDiffusionLength = washcoatThickness > 0.0 ? washcoatThickness : 0.5 * equivalentSolidLength;
  parameters.internalDiameter = internalChannelDimension;
  parameters.outerDiameter = outerDiameter;
  parameters.channelArea = channel.area;
  parameters.channelPerimeter = channel.perimeter;
  parameters.outerArea = outerArea;
  parameters.openArea = openArea;
  parameters.solidArea = solidArea;
  parameters.washcoatThickness = washcoatThickness;
  parameters.washcoatVolumePerChannelVolume = activeWashcoatVolumeRatio;

  // Hydraulic-diameter laminar-channel approximation.
  parameters.viscousPressureDropCoefficient =
      32.0 / std::max(channel.hydraulicDiameter * channel.hydraulicDiameter, tinyLength);
  parameters.inertialPressureDropCoefficient = 0.0;
}

Monolith::ChannelShape Monolith::parseChannelShape(const std::string& value)
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

  throw std::runtime_error(
      "Error: Geometry ChannelShape must be one of triangular, square, hexagonal, circular");
}

std::string Monolith::channelShapeName(ChannelShape channelShape)
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

const ShapeParameters& Geometry::shapeParameters() const noexcept
{
  return std::visit([](const auto& value) -> const ShapeParameters& { return value.shapeParameters(); },
                    value);
}

double Geometry::pressureGradient(double dynamicViscosity, double density, double velocity) const noexcept
{
  const ShapeParameters& shape = shapeParameters();
  return shape.viscousPressureDropCoefficient * dynamicViscosity * velocity +
         shape.inertialPressureDropCoefficient * density * velocity * velocity;
}
