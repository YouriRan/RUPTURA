#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>

/**
 * \brief A parameter value with an optional temperature dependence.
 *
 * Input files write a temperature-scalable number either as a plain number (Constant) or as a law object,
 * {"Value": V, "Law": "<name>", "<coefficient key>": c, "ReferenceTemperature": T_ref}. The value V is the
 * parameter at T_ref. The laws are
 *
 * - VantHoff:    V exp[Q/R (1/T - 1/T_ref)],  coefficient "Energy" Q in J/mol (Q > 0 grows as T drops);
 * - Arrhenius:   V exp[-E/R (1/T - 1/T_ref)], coefficient "ActivationEnergy" E in J/mol (E > 0 grows with T);
 * - Exponential: V exp[chi (1 - T/T_ref)],    coefficient "Chi" [-];
 * - Linear:      V + a (T - T_ref),           coefficient "Slope" a in units of V per K;
 * - Inverse:     V + alpha (1 - T_ref/T),     coefficient "Alpha" alpha in units of V.
 *
 * Without a reference temperature, VantHoff and Arrhenius use the uncentred forms V exp[Q/(RT)] and
 * V exp[-E/(RT)], with V the pre-exponential factor; the other laws require one. The exponents are bounded to
 * [-700, 700], so a law never returns inf or 0 from overflow.
 */
struct TemperatureLaw
{
  /**
   * \brief The temperature relation of the parameter.
   */
  enum class Type
  {
    Constant = 0,     ///< V
    VantHoff = 1,     ///< V exp[Q/R (1/T - 1/T_ref)]
    Arrhenius = 2,    ///< V exp[-E/R (1/T - 1/T_ref)]
    Exponential = 3,  ///< V exp[chi (1 - T/T_ref)]
    Linear = 4,       ///< V + a (T - T_ref)
    Inverse = 5       ///< V + alpha (1 - T_ref/T)
  };

  TemperatureLaw() noexcept = default;

  /**
   * \brief Constructs a constant law, the meaning of a plain number in the input (implicit, like the input).
   */
  TemperatureLaw(double constant) noexcept : value(constant) {}  // NOLINT(google-explicit-constructor)

  TemperatureLaw(Type type, double value, double coefficient,
                 std::optional<double> referenceTemperature = std::nullopt) noexcept
      : type(type), value(value), coefficient(coefficient), referenceTemperature(referenceTemperature)
  {
  }

  Type type{Type::Constant};  ///< The temperature relation.
  double value{0.0};          ///< Value at the reference temperature (pre-exponential factor without one).
  double coefficient{0.0};    ///< Q or E in J/mol, chi [-], slope [unit/K], or alpha [unit of the value].
  std::optional<double> referenceTemperature{std::nullopt};  ///< T_ref in K.
  /// The VantHoff energy was taken from the heat of adsorption (or reaction); it follows that heat when an
  /// adsorbent override changes it.
  bool energyInherited{false};

  /**
   * \brief Returns whether the value depends on temperature.
   */
  [[nodiscard]] bool isConstant() const noexcept { return type == Type::Constant; }

  /**
   * \brief Evaluates the parameter at a temperature in K.
   */
  [[nodiscard]] double operator()(double temperature) const;

  /**
   * \brief Returns a one-line description, e.g. "VantHoff, Q = 2.4860e+04 J/mol, T_ref = 3.1800e+02 K".
   */
  [[nodiscard]] std::string repr() const;

  /**
   * \brief Returns the bit of a law type in a mask of allowed laws.
   */
  static constexpr unsigned mask(Type type) noexcept { return 1u << static_cast<unsigned>(type); }

  /**
   * \brief Returns the input name of a law type ("VantHoff", ...).
   */
  static std::string_view name(Type type) noexcept;

  /**
   * \brief Returns the input key of the law's coefficient ("Energy", ...), or an empty view for Constant.
   */
  static std::string_view coefficientKey(Type type) noexcept;

  /**
   * \brief Returns the law type with the given input name (case-insensitive), if any.
   */
  static std::optional<Type> fromName(std::string_view name) noexcept;

  /**
   * \brief Returns whether the law needs a reference temperature (all but Constant, VantHoff, and Arrhenius).
   */
  static constexpr bool requiresReferenceTemperature(Type type) noexcept
  {
    return type == Type::Exponential || type == Type::Linear || type == Type::Inverse;
  }

  /**
   * \brief Returns the names of the laws in a mask, e.g. "VantHoff" or "Linear or Inverse".
   */
  static std::string describeMask(unsigned allowedLaws);
};

/**
 * \brief Describes whether, and how, one parameter position may depend on temperature.
 *
 * The isotherm tables (Isotherm::temperatureParameters) give one entry per entry of the Parameters array, and
 * the kinetic and reaction keys use the same description. A mask of zero means the parameter is constant.
 */
struct TemperatureParameterSpec
{
  std::string_view name;   ///< Name used in messages and the graphical interfaces, e.g. "affinity".
  unsigned allowedLaws{0};  ///< Mask (TemperatureLaw::mask) of the laws accepted besides a plain number.
  /// A VantHoff law without "Energy" inherits the heat of adsorption (the affinity of a single-affinity isotherm,
  /// or the equilibrium constant of a reaction).
  bool heatLinked{false};
};
