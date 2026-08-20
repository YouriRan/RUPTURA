#pragma once

#include <cstddef>
#include <cmath>
#include <optional>
#include <string_view>

enum struct GeometryKind
{
  PackedBed = 0,
  Monolith = 1
};

enum struct ChannelShape
{
  Triangular = 0,
  Square = 1,
  Hexagonal = 2,
  Circular = 3
};

/**
 * \brief Parses a user-facing monolith channel-shape name.
 */
[[nodiscard]] ChannelShape parseChannelShape(std::string_view value);

/**
 * \brief Returns the canonical user-facing name for a channel shape.
 */
[[nodiscard]] std::string_view channelShapeName(ChannelShape channelShape);

/**
 * \brief Area densities used by transport and energy balances.
 *
 * Each value is normalized by the phase volume named in the field.
 */
struct ContactAreaDensities
{
  double fluidSolidPerFluidVolume{0.0};  ///< Gas/liquid-solid area per flowing volume, 1/m.
  double solidFluidPerSolidVolume{0.0};  ///< Gas/liquid-solid area per active solid volume, 1/m.
  double fluidWallPerFluidVolume{0.0};   ///< Gas/liquid-wall area per flowing volume, 1/m.
  double wallInnerPerWallVolume{0.0};    ///< Inner wall area per wall volume, 1/m.
  double wallOuterPerWallVolume{0.0};    ///< Outer wall area per wall volume, 1/m.
};

/**
 * \brief Physical dimensions and diagnostic cross-section terms.
 */
struct GeometryDimensions
{
  ChannelShape channelShape{ChannelShape::Circular};
  std::size_t numberOfChannels{1};

  double hydraulicDiameter{0.0};          ///< Flow-channel hydraulic diameter, m.
  double solidCharacteristicLength{0.0};  ///< Bead diameter or equivalent solid length, m.
  double poreDiffusionLength{0.0};        ///< Characteristic pore/washcoat diffusion length, m.
  double internalDiameter{0.0};           ///< Tube diameter or monolith channel dimension, m.
  double outerDiameter{0.0};              ///< Column/monolith outer diameter, m.

  double channelArea{0.0};                     ///< Single-channel open area, m^2.
  double channelPerimeter{0.0};                ///< Single-channel wetted perimeter, m.
  double outerArea{0.0};                       ///< Whole bundle/tube outer cross-section area, m^2.
  double openArea{0.0};                        ///< Total flowing cross-section area, m^2.
  double solidArea{0.0};                       ///< Solid/wall cross-section area, m^2.
  double washcoatThickness{0.0};               ///< Optional washcoat thickness, m.
  double washcoatVolumePerChannelVolume{0.0};  ///< Active washcoat volume ratio.
};

/**
 * \brief Coefficients for a value-semantic pressure-drop law.
 */
struct PressureDropLaw
{
  double viscousCoefficient{0.0};   ///< Multiplies mu * u in dP/dz.
  double inertialCoefficient{0.0};  ///< Multiplies rho * u^2 in dP/dz.

  [[nodiscard]] double gradient(double dynamicViscosity, double density, double velocity) const noexcept
  {
    return viscousCoefficient * dynamicViscosity * velocity +
           inertialCoefficient * density * velocity * std::abs(velocity);
  }

  [[nodiscard]] double velocity(double dynamicViscosity, double density, double pressureGradient) const noexcept
  {
    const double drivingGradient = -pressureGradient;
    const double viscousResistance = viscousCoefficient * dynamicViscosity;
    const double inertialResistance = inertialCoefficient * density;
    const double magnitude = std::abs(drivingGradient);

    if (magnitude == 0.0) return 0.0;
    if (inertialResistance <= 0.0)
    {
      return viscousResistance > 0.0 ? drivingGradient / viscousResistance : 0.0;
    }
    if (viscousResistance <= 0.0)
    {
      return std::copysign(std::sqrt(magnitude / inertialResistance), drivingGradient);
    }

    // Stable positive root of B |v|^2 + A |v| - |dP/dz| = 0.
    const double speed = 2.0 * magnitude /
                         (viscousResistance +
                          std::sqrt(viscousResistance * viscousResistance + 4.0 * inertialResistance * magnitude));
    return std::copysign(speed, drivingGradient);
  }
};

/**
 * \brief Fully derived geometry stored and consumed by a simulation.
 *
 * Geometry is an ordinary copyable value. Geometry-specific validation and
 * formulas live in the makeGeometry overloads; downstream calculations do not
 * branch on the source geometry type.
 */
struct Geometry
{
  Geometry(GeometryKind kind, double voidFraction, double solidToFluidVolumeRatio, ContactAreaDensities contactAreas,
           GeometryDimensions dimensions, PressureDropLaw pressureDrop) noexcept
      : kind(kind),
        voidFraction(voidFraction),
        solidToFluidVolumeRatio(solidToFluidVolumeRatio),
        contactAreas(contactAreas),
        dimensions(dimensions),
        pressureDrop(pressureDrop)
  {
  }

  GeometryKind kind;
  double voidFraction;             ///< Flowing-volume fraction.
  double solidToFluidVolumeRatio;  ///< Active solid volume divided by flowing volume.

  ContactAreaDensities contactAreas;
  GeometryDimensions dimensions;
  PressureDropLaw pressureDrop;

  [[nodiscard]] double loadingPrefactor(double particleDensity) const noexcept
  {
    return particleDensity * solidToFluidVolumeRatio;
  }
};

/**
 * \brief User inputs for a packed bed.
 */
struct PackedBedTubeSpec
{
  double voidFraction{0.4};
  double particleDiameter{1.0e-3};
  double internalDiameter{0.0};
  double outerDiameter{0.0};
};

/**
 * \brief User inputs for a monolith bundle with identical coated channels.
 */
struct MonolithSpec
{
  ChannelShape channelShape{ChannelShape::Circular};
  double internalChannelDimension{0.0};
  double outerDiameter{0.0};
  std::size_t numberOfChannels{0};
  double washcoatThickness{0.0};
  std::optional<double> washcoatVolumePerChannelVolume;
  double forchheimerCoefficient{0.0};  ///< Channel Forchheimer coefficient b_M, 1/m.
};

[[nodiscard]] Geometry makeGeometry(const PackedBedTubeSpec& specification);
[[nodiscard]] Geometry makeGeometry(const MonolithSpec& specification);
