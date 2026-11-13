#include "cycle/spatial.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>
namespace cycle {
namespace {
std::string normalize(std::string v, bool upper = true) {
  auto first = v.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  v = v.substr(first, v.find_last_not_of(" \t\r\n") - first + 1);
  std::transform(v.begin(), v.end(), v.begin(), [upper](unsigned char c) {
    return upper ? std::toupper(c) : (c == '-' ? '_' : std::tolower(c));
  });
  return v;
}
std::string canonical_scheme(std::string s) {
  s = normalize(s);
  if (s == "FIRST_ORDER") s = "UPWIND";
  if (s == "MUSCL") s = "TVD";
  if (s == "WENO5") s = "WENO";
  if (s != "UPWIND" && s != "TVD" && s != "WENO" && s != "WENO3")
    throw std::invalid_argument("ConvectionScheme must be Upwind, TVD, WENO3 or WENO5");
  return s;
}
void finite(const Vec& values, const char* message) {
  for (double x : values)
    if (!std::isfinite(x)) throw std::invalid_argument(message);
}
bool crosses(const Vec& edges, std::size_t first, std::size_t last, const Vec& interfaces) {
  for (double x : interfaces)
    if (edges[first] < x && x < edges[last + 1]) return true;
  return false;
}
double minmod(double a, double b) {
  return a * b > 0 ? std::copysign(std::min(std::abs(a), std::abs(b)), a) : 0;
}
using Matrix = std::vector<Vec>;
Matrix inverse(Matrix a) {
  auto n = a.size();
  Matrix out(n, Vec(n));
  for (std::size_t i = 0; i < n; ++i) out[i][i] = 1;
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t pivot = k;
    for (std::size_t i = k + 1; i < n; ++i)
      if (std::abs(a[i][k]) > std::abs(a[pivot][k])) pivot = i;
    if (!std::isfinite(a[pivot][k]) || std::abs(a[pivot][k]) < 1e-300)
      throw std::invalid_argument("Singular WENO geometry");
    std::swap(a[pivot], a[k]);
    std::swap(out[pivot], out[k]);
    double d = a[k][k];
    for (std::size_t j = 0; j < n; ++j) {
      a[k][j] /= d;
      out[k][j] /= d;
    }
    for (std::size_t i = 0; i < n; ++i)
      if (i != k) {
        double f = a[i][k];
        for (std::size_t j = 0; j < n; ++j) {
          a[i][j] -= f * a[k][j];
          out[i][j] -= f * out[k][j];
        }
      }
  }
  return out;
}
double norm2(const Matrix& a) {
  auto n = a.size();
  Matrix s(n, Vec(n));
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t k = 0; k < n; ++k) s[i][j] += a[k][i] * a[k][j];
  // Jacobi rotations evaluate the spectral norm for the small WENO matrices.
  for (std::size_t iteration = 0; iteration < 100 * n * n; ++iteration) {
    std::size_t p = 0, q = 1;
    double largest = 0;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = i + 1; j < n; ++j)
        if (std::abs(s[i][j]) > largest) {
          largest = std::abs(s[i][j]);
          p = i;
          q = j;
        }
    double scale = 0;
    for (std::size_t i = 0; i < n; ++i) scale = std::max(scale, std::abs(s[i][i]));
    if (largest <= 1e-15 * std::max(scale, 1e-300)) break;
    double angle = .5 * std::atan2(2 * s[p][q], s[q][q] - s[p][p]), c = std::cos(angle),
           t = std::sin(angle), pp = s[p][p], qq = s[q][q], pq = s[p][q];
    for (std::size_t k = 0; k < n; ++k)
      if (k != p && k != q) {
        double kp = s[k][p], kq = s[k][q];
        s[k][p] = s[p][k] = c * kp - t * kq;
        s[k][q] = s[q][k] = t * kp + c * kq;
      }
    s[p][p] = c * c * pp - 2 * c * t * pq + t * t * qq;
    s[q][q] = t * t * pp + 2 * c * t * pq + c * c * qq;
    s[p][q] = s[q][p] = 0;
  }
  double maximum = 0;
  for (std::size_t i = 0; i < n; ++i) maximum = std::max(maximum, s[i][i]);
  return std::sqrt(maximum);
}
Matrix average_matrix(const Vec& edges, int first, int count, double origin, double scale) {
  Matrix a(count, Vec(count));
  for (int i = 0; i < count; ++i) {
    double lo = (edges[first + i] - origin) / scale, hi = (edges[first + i + 1] - origin) / scale;
    for (int k = 0; k < count; ++k)
      a[i][k] = (std::pow(hi, k + 1) - std::pow(lo, k + 1)) / ((k + 1) * (hi - lo));
  }
  return a;
}
struct WenoRecord {
  bool active = false;
  int first = 0, full = 5, candidates = 3, local = 3;
  std::array<std::array<double, 3>, 3> maps{};
  std::array<std::array<std::array<double, 3>, 3>, 3> beta{};
  std::array<double, 3> ideal{};
};
std::vector<WenoRecord> build_weno(const Vec& z, int order, bool nodal) {
  auto edges = control_volume_faces(z);
  int n = static_cast<int>(z.size());
  // FD-WENO reconstructs a numerical flux from point values. On a uniform
  // nodal grid its deconvolution moments equal FV average moments on cells
  // centred at those nodes. The physical end control volumes remain half-width.
  // Nonuniform FD stencils use the limited nodal reconstruction instead.
  if (nodal) {
    edges.front() = z.front() - .5 * (z[1] - z[0]);
    edges.back() = z.back() + .5 * (z[n-1] - z[n-2]);
  }
  const int local_count = order == 3 ? 2 : 3, full_count = 2 * local_count - 1;
  std::vector<WenoRecord> result(2 * (n + 1));
  for (int face = 1; face < n; ++face)
    for (int direction : {1, -1}) {
      int up = direction > 0 ? face - 1 : face, first = up - local_count + 1;
      if (first < 0 || first + full_count > n) continue;
      if (nodal) {
        const double h = z[first+1] - z[first];
        bool uniform = true;
        for (int j=1;j+1<full_count;++j)
          uniform &= std::abs((z[first+j+1]-z[first+j])-h) <= 1e-9*h;
        if (!uniform) continue;
      }
      double origin = edges[face], scale = edges[up + 1] - edges[up];
      auto a = average_matrix(edges, first, full_count, origin, scale), inv = inverse(a);
      if (norm2(a) * norm2(inv) > 1e13)
        throw std::invalid_argument(
            "WENO stencil is ill-conditioned; reduce adjacent grid-size ratios");
      // A^T optimal=e0 gives the first row of A^{-1}.
      std::array<double, 5> optimal{};
      for (int i = 0; i < full_count; ++i) optimal[i] = inv[0][i];
      std::array<std::array<double, 3>, 5> coefficients{};
      WenoRecord rec;
      rec.active = true;
      rec.first = first;
      rec.full = full_count; rec.candidates = local_count; rec.local = local_count;
      for (int start = 0; start < local_count; ++start) {
        auto local = inverse(average_matrix(edges, first + start, local_count, origin, scale));
        for (int j = 0; j < local_count; ++j) {
          rec.maps[start][j] = local[0][j];
          coefficients[start + j][start] = local[0][j];
        }
        double lo = (edges[up] - origin) / scale, hi = (edges[up + 1] - origin) / scale;
        std::array<std::array<double, 3>, 3> smooth{};
        for (int derivative = 1; derivative < local_count; ++derivative)
          for (int p = derivative; p < local_count; ++p)
            for (int q = derivative; q < local_count; ++q) {
              int exponent = p + q - 2 * derivative + 1;
              double factor = 1;
              for (int d = 0; d < derivative; ++d) factor *= (p - d) * (q - d);
              smooth[p][q] += factor * (std::pow(hi, exponent) - std::pow(lo, exponent)) / exponent;
            }
        for (int i = 0; i < local_count; ++i)
          for (int j = 0; j < local_count; ++j)
            for (int p = 0; p < local_count; ++p)
              for (int q = 0; q < local_count; ++q)
                rec.beta[start][i][j] += local[p][i] * smooth[p][q] * local[q][j];
      }
      Matrix normal(local_count, Vec(local_count));
      std::array<double, 3> rhs{};
      for (int i = 0; i < local_count; ++i)
        for (int row = 0; row < full_count; ++row) {
          rhs[i] += coefficients[row][i] * optimal[row];
          for (int j = 0; j < local_count; ++j) normal[i][j] += coefficients[row][i] * coefficients[row][j];
        }
      auto ni = inverse(normal);
      for (int i = 0; i < local_count; ++i)
        for (int j = 0; j < local_count; ++j) rec.ideal[i] += ni[i][j] * rhs[j];
      for (int row = 0; row < full_count; ++row) {
        double actual = 0;
        for (int k = 0; k < local_count; ++k) actual += coefficients[row][k] * rec.ideal[k];
        if (std::abs(actual - optimal[row]) > 1e-8)
          throw std::invalid_argument("Could not construct consistent nonuniform WENO weights");
      }
      result[2 * face + (direction < 0)] = rec;
    }
  return result;
}
const std::vector<WenoRecord>& weno_records(const Vec& z, int order, bool nodal) {
  using Key=std::tuple<Vec,int,bool>;
  static thread_local std::map<Key, std::vector<WenoRecord>> cache;
  Key key{z,order,nodal};
  auto it = cache.find(key);
  if (it != cache.end()) return it->second;
  if (cache.size() >= 32) cache.erase(cache.begin());
  return cache.emplace(std::move(key), build_weno(z,order,nodal)).first->second;
}
std::array<double, 3> nonlinear_weights(const std::array<double, 3>& ideal,
                                        const std::array<double, 3>& beta) {
  std::array<double, 3> positive{}, negative{}, result{};
  for (int k = 0; k < 3; ++k) {
    positive[k] = .5 * (ideal[k] + 3 * std::abs(ideal[k]));
    negative[k] = positive[k] - ideal[k];
  }
  for (int sign : {1, -1}) {
    auto& part = sign > 0 ? positive : negative;
    double total = std::accumulate(part.begin(), part.end(), 0.0);
    if (total == 0) continue;
    std::array<double, 3> alpha{};
    double alpha_total = 0;
    for (int k = 0; k < 3; ++k) {
      alpha[k] = (part[k] / total) / std::pow(1e-12 + beta[k], 2);
      alpha_total += alpha[k];
    }
    for (int k = 0; k < 3; ++k) result[k] += sign * total * alpha[k] / alpha_total;
  }
  return result;
}
}  // namespace
std::pair<std::string, std::string> spatial_settings(const Json& cfg) {
  std::string method = normalize(cfg.value("SPATIAL_METHOD", std::string("FVM")));
  if (method == "FV") method = "FVM";
  if (method == "FD") method = "FDM";
  if (method != "FVM" && method != "FDM") throw std::invalid_argument("SpatialDiscretization must be FDM or FVM");
  return {method, canonical_scheme(cfg.value("ADVECTION_SCHEME", std::string("UPWIND")))};
}
std::size_t stencil_radius(const Json& cfg) {
  auto scheme = spatial_settings(cfg).second;
  return scheme == "WENO" ? 3 : (scheme == "TVD" || scheme == "WENO3") ? 2 : 1;
}
Vec interface_locations(const Json& cfg) {
  Vec points, lengths, interfaces;
  if (cfg.contains("SPATIAL_INTERFACE_LOCATIONS") &&
      !cfg.at("SPATIAL_INTERFACE_LOCATIONS").is_null())
    points = cfg.at("SPATIAL_INTERFACE_LOCATIONS").get<Vec>();
  if (cfg.contains("L_ads") && !cfg.at("L_ads").is_null()) lengths = cfg.at("L_ads").get<Vec>();
  if (cfg.contains("L_int") && !cfg.at("L_int").is_null()) interfaces = cfg.at("L_int").get<Vec>();
  double location = 0;
  for (std::size_t i = 0; i + 1 < lengths.size(); ++i) {
    location += lengths[i];
    points.push_back(location);
    double thickness = i < interfaces.size() ? interfaces[i] : 0;
    if (thickness > 0) {
      location += thickness;
      points.push_back(location);
    }
  }
  finite(points, "SPATIAL_INTERFACE_LOCATIONS must be finite");
  std::sort(points.begin(), points.end());
  points.erase(std::unique(points.begin(), points.end()), points.end());
  return points;
}
Vec control_volume_faces(const Vec& z) {
  if (z.size() < 2)
    throw std::invalid_argument("Axial coordinates must be finite and strictly increasing");
  finite(z, "Axial coordinates must be finite and strictly increasing");
  Vec edges(z.size() + 1);
  edges.front() = z.front();
  edges.back() = z.back();
  for (std::size_t i = 1; i < z.size(); ++i) {
    if (z[i] <= z[i - 1])
      throw std::invalid_argument("Axial coordinates must be finite and strictly increasing");
    edges[i] = .5 * (z[i - 1] + z[i]);
  }
  return edges;
}
Vec reconstructed_face_values(const Vec& data, const Vec& z, const Vec& flow, const Vec& left,
                              const Vec& right, const std::string& requested_scheme,
                              const Vec& interfaces, const std::string& method) {
  auto n = z.size();
  auto edges = control_volume_faces(z);
  if (data.size() % n) throw std::invalid_argument("values must have shape (n_fields, Nz)");
  auto fields = data.size() / n;
  if (flow.size() != n + 1 || left.size() != fields || right.size() != fields)
    throw std::invalid_argument("Invalid face flow or reservoir value shape");
  for (const auto* a : {&data, &flow, &left, &right})
    finite(*a, "face reconstruction requires finite values");
  std::string scheme = canonical_scheme(requested_scheme);
  Vec centers(n), slope(data.size(), 0);
  const bool nodal = normalize(method) == "FDM";
  for (std::size_t j = 0; j < n; ++j) centers[j] = nodal ? z[j] : .5 * (edges[j] + edges[j + 1]);
  if (scheme != "UPWIND" && n > 2)
    for (std::size_t j = 1; j + 1 < n; ++j)
      if (!crosses(edges, j - 1, j + 1, interfaces))
        for (std::size_t field = 0; field < fields; ++field) {
          auto o = field * n;
          slope[o + j] = minmod((data[o + j] - data[o + j - 1]) / (centers[j] - centers[j - 1]),
                                (data[o + j + 1] - data[o + j]) / (centers[j + 1] - centers[j]));
        }
  const std::vector<WenoRecord>* records = (scheme == "WENO" || scheme == "WENO3")
      ? &weno_records(z,scheme == "WENO3" ? 3 : 5,nodal) : nullptr;
  Vec output(fields * (n + 1));
  for (std::size_t face = 1; face < n; ++face) {
    int direction = flow[face] >= 0 ? 1 : -1;
    std::size_t up = direction > 0 ? face - 1 : face;
    const WenoRecord* rec = records ? &records->at(2 * face + (direction < 0)) : nullptr;
    bool use = rec && rec->active && !crosses(edges, rec->first, rec->first + rec->full - 1, interfaces);
    for (std::size_t field = 0; field < fields; ++field) {
      auto o = field * n;
      double candidate = data[o + up] + slope[o + up] * (edges[face] - centers[up]);
      if (use) {
        std::array<double, 5> stencil{};
        double magnitude = 1e-100;
        for (int j = 0; j < rec->full; ++j) {
          stencil[j] = data[o + rec->first + j];
          magnitude = std::max(magnitude, std::abs(stencil[j]));
        }
        for (double& v : stencil) v /= magnitude;
        std::array<double, 3> estimates{}, beta{};
        for (int k = 0; k < rec->candidates; ++k) {
          for (int i = 0; i < rec->local; ++i) {
            estimates[k] += rec->maps[k][i] * stencil[k + i];
            for (int j = 0; j < rec->local; ++j)
              beta[k] += stencil[k + i] * rec->beta[k][i][j] * stencil[k + j];
          }
          beta[k] = std::max(beta[k], 0.0);
        }
        auto weights = nonlinear_weights(rec->ideal, beta);
        candidate = 0;
        for (int k = 0; k < 3; ++k) candidate += weights[k] * estimates[k];
        candidate *= magnitude;
      }
      output[field * (n + 1) + face] =
          std::min(std::max(candidate, 0.0), 2 * std::max(data[o + up], 0.0));
    }
  }
  for (std::size_t field = 0; field < fields; ++field) {
    output[field * (n + 1)] = flow.front() >= 0 ? left[field] : data[field * n];
    output[field * (n + 1) + n] = flow.back() < 0 ? right[field] : data[field * n + n - 1];
  }
  return output;
}
Vec reconstruct(const Vec& data, std::size_t fields, std::size_t nz, const Vec& z, const Vec& flow,
                const Vec& left, const Vec& right, const std::string& scheme,
                const Vec& interfaces) {
  if (z.size() != nz || data.size() != fields * nz)
    throw std::invalid_argument("Invalid reconstruction data shape");
  return reconstructed_face_values(data, z, flow, left, right, scheme, interfaces);
}
std::pair<Vec, Vec> conservative_species_face_flux_rhs(
    const Vec& c, const Vec& flow, const Vec& sink, const Vec& eps, const Vec& dispersion,
    const Vec& left, const Vec& right, const Vec& z, const Vec& eps_face, const std::string& scheme,
    const Vec& interfaces, const Vec& supplied_widths, const std::string& requested_closure,
    const Vec& left_flux, const Vec& right_flux, const std::string& method) {
  auto edges = control_volume_faces(z);
  auto n = z.size();
  if (c.size() % n) throw std::invalid_argument("concentration must have shape (n_comp, Nz)");
  auto nc = c.size() / n;
  if (eps.size() != n || sink.size() != c.size() || flow.size() != n + 1)
    throw std::invalid_argument("Invalid transport input shape");
  if (dispersion.size() != 1 && dispersion.size() != nc && dispersion.size() != c.size())
    throw std::invalid_argument("dispersion must be scalar, (n_comp,), or (n_comp, Nz)");
  if (!supplied_widths.empty() && supplied_widths.size() != n)
    throw std::invalid_argument("cell_widths must match axial control-volume geometry");
  Vec widths(n), centers(n);
  for (std::size_t j = 0; j < n; ++j) {
    widths[j] = edges[j + 1] - edges[j];
    centers[j] = normalize(method) == "FDM" ? z[j] : .5 * (edges[j] + edges[j + 1]);
    if (!supplied_widths.empty() &&
        std::abs(widths[j] - supplied_widths[j]) > 1e-14 + 1e-10 * std::abs(widths[j]))
      throw std::invalid_argument("cell_widths must match axial control-volume geometry");
    if (!(eps[j] > 0))
      throw std::invalid_argument("fluid_capacity must contain Nz positive values");
  }
  if (!eps_face.empty()) {
    if (eps_face.size() != flow.size())
      throw std::invalid_argument("face_fluid_capacity must contain Nz+1 positive values");
    for (double v : eps_face)
      if (!(v > 0) || !std::isfinite(v))
        throw std::invalid_argument("face_fluid_capacity must contain Nz+1 positive values");
  }
  for (const auto* a : {&c, &eps, &sink, &flow, &dispersion})
    finite(*a, "transport inputs must be finite");
  for (double d : dispersion)
    if (d < 0) throw std::invalid_argument("dispersion must be nonnegative");
  auto flux = reconstructed_face_values(c, z, flow, left, right, scheme, interfaces, method);
  for (std::size_t i = 0; i < nc; ++i)
    for (std::size_t f = 0; f < n + 1; ++f) flux[i * (n + 1) + f] *= flow[f];
  Vec resistance(nc * (n - 1), std::numeric_limits<double>::infinity());
  auto diffusion = [&](std::size_t i, std::size_t j) {
    return dispersion.size() == 1    ? dispersion[0]
           : dispersion.size() == nc ? dispersion[i]
                                     : dispersion[i * n + j];
  };
  for (std::size_t i = 0; i < nc; ++i)
    for (std::size_t j = 0; j + 1 < n; ++j) {
      double kl = eps[j] * diffusion(i, j), kr = eps[j + 1] * diffusion(i, j + 1);
      if (kl > 0 && kr > 0)
        resistance[i * (n - 1) + j] =
            (edges[j + 1] - centers[j]) / kl + (centers[j + 1] - edges[j + 1]) / kr;
    }
  auto closure = normalize(requested_closure, false);
  Vec dispersive(nc * (n - 1), 0);
  if (closure == "concentration" || closure == "concentration_gradient")
    for (std::size_t i = 0; i < nc; ++i)
      for (std::size_t j = 0; j + 1 < n; ++j)
        dispersive[i * (n - 1) + j] =
            -(c[i * n + j + 1] - c[i * n + j]) / resistance[i * (n - 1) + j];
  else if (closure == "mole_fraction" || closure == "mixture" || closure == "mixture_averaged") {
    Vec total(n, 0), fractions(c.size(), 0);
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t i = 0; i < nc; ++i) total[j] += std::max(c[i * n + j], 0.0);
      if (total[j] > 0)
        for (std::size_t i = 0; i < nc; ++i)
          fractions[i * n + j] = std::max(c[i * n + j], 0.0) / total[j];
    }
    for (std::size_t j = 0; j + 1 < n; ++j) {
      double face_total = .5 * (total[j] + total[j + 1]), fraction_sum = 0, flux_sum = 0;
      for (std::size_t i = 0; i < nc; ++i) {
        double d = -(fractions[i * n + j + 1] - fractions[i * n + j]) * face_total /
                   resistance[i * (n - 1) + j];
        dispersive[i * (n - 1) + j] = d;
        flux_sum += d;
        fraction_sum += .5 * (fractions[i * n + j] + fractions[i * n + j + 1]);
      }
      if (fraction_sum > 0)
        for (std::size_t i = 0; i < nc; ++i)
          dispersive[i * (n - 1) + j] -=
              .5 * (fractions[i * n + j] + fractions[i * n + j + 1]) / fraction_sum * flux_sum;
    }
  } else
    throw std::invalid_argument("dispersion_closure must be concentration or mole_fraction");
  for (std::size_t i = 0; i < nc; ++i)
    for (std::size_t j = 0; j + 1 < n; ++j)
      flux[i * (n + 1) + j + 1] += dispersive[i * (n - 1) + j];
  for (int side = 0; side < 2; ++side) {
    const auto& values = side ? right_flux : left_flux;
    if (values.empty()) continue;
    if (values.size() != nc)
      throw std::invalid_argument("Boundary component flux must have n_comp finite entries");
    finite(values, "Boundary component flux must have n_comp finite entries");
    for (std::size_t i = 0; i < nc; ++i) flux[i * (n + 1) + (side ? n : 0)] = values[i];
  }
  Vec rhs(c.size());
  for (std::size_t i = 0; i < nc; ++i)
    for (std::size_t j = 0; j < n; ++j)
      rhs[i * n + j] = -(flux[i * (n + 1) + j + 1] - flux[i * (n + 1) + j]) / (eps[j] * widths[j]) -
                       sink[i * n + j] / eps[j];
  return {std::move(rhs), std::move(flux)};
}
}  // namespace cycle
