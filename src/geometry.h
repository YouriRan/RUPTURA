#pragma once

#include <cstddef>
#include <string>
#include <variant>

/**
 * \brief Precomputed geometry terms used by transport and energy balances.
 *
 * Area densities are normalized by the phase volume named in the field. For
 * example, fluidSolidContactAreaPerFluidVolume is the gas/liquid-solid contact
 * area per flowing-channel volume.
 */
struct ShapeParameters
{
  enum struct GeometryKind
  {
    HollowTube = 0,
    Monolith = 1
  };

  GeometryKind geometryKind{GeometryKind::HollowTube};
  std::string geometryName{"HollowTube"};
  std::string channelShape{"circular"};
  std::size_t numberOfChannels{1};

  double voidFraction{0.4};                         ///< Flowing-volume fraction.
  double solidToFluidVolumeRatio{1.5};              ///< Solid volume divided by flowing volume.
  double fluidSolidContactAreaPerFluidVolume{0.0};  ///< Gas/liquid-solid area density, 1/m.
  double solidFluidContactAreaPerSolidVolume{0.0};  ///< Gas/liquid-solid area per solid volume, 1/m.
  double fluidWallContactAreaPerFluidVolume{0.0};   ///< Gas/liquid-wall area density, 1/m.
  double wallInnerContactAreaPerWallVolume{0.0};    ///< Inner wall area per wall volume, 1/m.
  double wallOuterContactAreaPerWallVolume{0.0};    ///< Outer wall area per wall volume, 1/m.

  double hydraulicDiameter{0.0};       ///< Flow-channel hydraulic diameter, m.
  double solidCharacteristicLength{0.0};  ///< Bead diameter or equivalent solid length, m.
  double poreDiffusionLength{0.0};     ///< Characteristic pore/washcoat diffusion length, m.
  double internalDiameter{0.0};        ///< Hollow tube inner diameter or channel characteristic size, m.
  double outerDiameter{0.0};           ///< Column/monolith outer diameter, m.

  double channelArea{0.0};       ///< Single-channel open area, m^2.
  double channelPerimeter{0.0};  ///< Single-channel wetted perimeter, m.
  double outerArea{0.0};         ///< Whole bundle/tube outer cross-section area, m^2.
  double openArea{0.0};          ///< Total flowing cross-section area, m^2.
  double solidArea{0.0};         ///< Solid/wall cross-section area, m^2.
  double washcoatThickness{0.0};  ///< Optional washcoat thickness, m.
  double washcoatVolumePerChannelVolume{0.0};  ///< Optional washcoat volume ratio.

  double viscousPressureDropCoefficient{0.0};  ///< Multiplies mu * u in dP/dz.
  double inertialPressureDropCoefficient{0.0};  ///< Multiplies rho * u^2 in dP/dz.

  [[nodiscard]] double loadingPrefactor(double particleDensity) const noexcept
  {
    return particleDensity * solidToFluidVolumeRatio;
  }
};

/**
 * \brief Packed beads in one hollow tube.
 */
struct HollowTube
{
  HollowTube(double voidFraction = 0.4, double particleDiameter = 1.0e-3,
             double internalDiameter = 0.0, double outerDiameter = 0.0);

  [[nodiscard]] const ShapeParameters& shapeParameters() const noexcept { return parameters; }

  ShapeParameters parameters;
};

/**
 * \brief Monolith bundle with many identical coated flow channels.
 */
struct Monolith
{
  enum struct ChannelShape
  {
    Triangular = 0,
    Square = 1,
    Hexagonal = 2,
    Circular = 3
  };

  Monolith(ChannelShape channelShape, double internalChannelDimension, double outerDiameter,
           std::size_t numberOfChannels, double washcoatThickness = 0.0,
           double washcoatVolumePerChannelVolume = -1.0);

  [[nodiscard]] const ShapeParameters& shapeParameters() const noexcept { return parameters; }

  [[nodiscard]] static ChannelShape parseChannelShape(const std::string& value);
  [[nodiscard]] static std::string channelShapeName(ChannelShape channelShape);

  ShapeParameters parameters;
};

/**
 * \brief Value-semantic geometry object stored by a Column.
 */
struct Geometry
{
  Geometry() = default;
  Geometry(const HollowTube& hollowTube) : value(hollowTube) {}
  Geometry(const Monolith& monolith) : value(monolith) {}

  [[nodiscard]] const ShapeParameters& shapeParameters() const noexcept;
  [[nodiscard]] double pressureGradient(double dynamicViscosity, double density, double velocity) const noexcept;

  std::variant<HollowTube, Monolith> value{HollowTube{}};
};
