#include "cycle/npz.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <regex>
#include <sstream>
#include <zlib.h>
namespace cycle {
namespace {
using Bytes = std::vector<unsigned char>;
constexpr std::uint64_t zip_limit = std::numeric_limits<std::uint32_t>::max();
void put16(Bytes& b, std::uint16_t x) {
  b.push_back(x & 255);
  b.push_back((x >> 8) & 255);
}
void put32(Bytes& b, std::uint32_t x) {
  for (int i = 0; i < 4; ++i) b.push_back((x >> (8 * i)) & 255);
}
void put64(Bytes& b, std::uint64_t x) {
  for (int i = 0; i < 8; ++i) b.push_back((x >> (8 * i)) & 255);
}
void putdouble(Bytes& b, double x) {
  static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
  std::uint64_t bits;
  std::memcpy(&bits, &x, 8);
  put64(b, bits);
}
std::size_t count(const std::vector<std::size_t>& shape) {
  std::size_t n = 1;
  for (auto s : shape) {
    if (s && n > std::numeric_limits<std::size_t>::max() / s)
      throw std::invalid_argument("NPY shape product overflow");
    n *= s;
  }
  return n;
}
std::vector<std::uint32_t> unicode(const std::string& text) {
  std::vector<std::uint32_t> out;
  for (std::size_t i = 0; i < text.size();) {
    auto lead = static_cast<unsigned char>(text[i++]);
    std::uint32_t c = 0, minimum = 0;
    int more = 0;
    if (lead < 128)
      c = lead;
    else if (lead >= 0xc2 && lead <= 0xdf) {
      c = lead & 31;
      more = 1;
      minimum = 0x80;
    } else if (lead >= 0xe0 && lead <= 0xef) {
      c = lead & 15;
      more = 2;
      minimum = 0x800;
    } else if (lead >= 0xf0 && lead <= 0xf4) {
      c = lead & 7;
      more = 3;
      minimum = 0x10000;
    } else
      throw std::invalid_argument("Invalid UTF-8 in NPZ text");
    if (i + more > text.size()) throw std::invalid_argument("Truncated UTF-8 in NPZ text");
    for (int j = 0; j < more; ++j) {
      auto next = static_cast<unsigned char>(text[i++]);
      if ((next & 0xc0) != 0x80) throw std::invalid_argument("Invalid UTF-8 continuation");
      c = (c << 6) | (next & 63);
    }
    if (c < minimum || c > 0x10ffff || (c >= 0xd800 && c <= 0xdfff))
      throw std::invalid_argument("Invalid Unicode code point");
    out.push_back(c);
  }
  return out;
}
Bytes header(const std::string& dtype, const std::vector<std::size_t>& shape, std::size_t size) {
  std::string dims = "(";
  for (std::size_t i = 0; i < shape.size(); ++i) {
    if (i) dims += ", ";
    dims += std::to_string(shape[i]);
  }
  if (shape.size() == 1) dims += ",";
  dims += ")";
  std::string dict = "{'descr': '" + dtype + "', 'fortran_order': False, 'shape': " + dims + ", }";
  auto padded = [&](std::size_t prefix) {
    return ((prefix + dict.size() + 1 + 63) / 64) * 64 - prefix;
  };
  int major = 1;
  std::size_t prefix = 10, length = padded(prefix);
  if (length > 65535) {
    major = 2;
    prefix = 12;
    length = padded(prefix);
  }
  if (length > zip_limit || static_cast<std::uint64_t>(size) + prefix + length > zip_limit)
    throw std::invalid_argument("NPY entry exceeds ZIP32 limit; split archive");
  Bytes b;
  b.reserve(prefix + length + size);
  b.insert(b.end(), {0x93, 'N', 'U', 'M', 'P', 'Y', static_cast<unsigned char>(major), 0});
  if (major == 1)
    put16(b, static_cast<std::uint16_t>(length));
  else
    put32(b, static_cast<std::uint32_t>(length));
  b.insert(b.end(), dict.begin(), dict.end());
  b.insert(b.end(), length - dict.size() - 1, ' ');
  b.push_back('\n');
  return b;
}
Bytes encode(const NpyArray& a) {
  auto n = count(a.shape);
  if (n != a.data.size()) throw std::invalid_argument("NPY data length does not match shape");
  if (n > zip_limit / 8) throw std::invalid_argument("NPY double data exceeds ZIP32 limit");
  auto b = header("<f8", a.shape, n * 8);
  for (double x : a.data) putdouble(b, x);
  return b;
}
struct Leaves {
  std::vector<const Json*> values;
  bool text = false, numeric = false, floating = false, integral = false, boolean = false,
       negative = false, large_unsigned = false;
  void add(const Json& v) {
    values.push_back(&v);
    if (v.is_string())
      text = true;
    else if (v.is_boolean()) {
      numeric = true;
      boolean = true;
    } else if (v.is_number_float()) {
      numeric = true;
      floating = true;
    } else if (v.is_number_unsigned()) {
      numeric = true;
      integral = true;
      large_unsigned |= v.get<std::uint64_t>() >
                        static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
    } else if (v.is_number_integer()) {
      numeric = true;
      integral = true;
      negative |= v.get<std::int64_t>() < 0;
    } else
      throw std::invalid_argument("NPZ leaves must be numbers, bool or text");
  }
};
std::vector<std::size_t> infer(const Json& a, Leaves& leaves) {
  if (!a.is_array()) {
    leaves.add(a);
    return {};
  }
  std::vector<std::size_t> shape{a.size()}, sub;
  bool first = true;
  for (const auto& v : a) {
    auto s = infer(v, leaves);
    if (first) {
      sub = std::move(s);
      first = false;
    } else if (s != sub)
      throw std::invalid_argument("NPZ JSON array is ragged");
  }
  shape.insert(shape.end(), sub.begin(), sub.end());
  return shape;
}
Bytes encode(const Json& a) {
  Leaves l;
  auto shape = infer(a, l);
  auto n = count(shape);
  if (n != l.values.size()) throw std::invalid_argument("NPZ shape mismatch");
  if (l.text && l.numeric)
    throw std::invalid_argument("NPZ array cannot mix text and numeric values");
  if (l.text) {
    std::vector<std::vector<std::uint32_t>> strings;
    std::size_t width = 1;
    for (auto v : l.values) {
      strings.push_back(unicode(v->get<std::string>()));
      width = std::max(width, strings.back().size());
    }
    if (width > zip_limit / 4 || (n && width > zip_limit / 4 / n))
      throw std::invalid_argument("NPZ Unicode data exceeds ZIP32 limit");
    auto b = header("<U" + std::to_string(width), shape, n * width * 4);
    for (const auto& s : strings) {
      for (auto c : s) put32(b, c);
      for (std::size_t i = s.size(); i < width; ++i) put32(b, 0);
    }
    return b;
  }
  if (l.boolean && !l.integral && !l.floating) {
    auto b = header("|b1", shape, n);
    for (auto v : l.values) b.push_back(v->get<bool>() ? 1 : 0);
    return b;
  }
  if (n > zip_limit / 8) throw std::invalid_argument("NPZ numeric data exceeds ZIP32 limit");
  if (l.floating || l.values.empty()) {
    auto b = header("<f8", shape, n * 8);
    for (auto v : l.values)
      putdouble(b, v->is_boolean() ? (v->get<bool>() ? 1.0 : 0.0) : v->get<double>());
    return b;
  }
  if (l.large_unsigned && l.negative)
    throw std::invalid_argument(
        "NPZ array mixes negative integers and uint64 values above INT64_MAX");
  auto b = header(l.large_unsigned ? "<u8" : "<i8", shape, n * 8);
  for (auto v : l.values) {
    if (v->is_boolean())
      put64(b, v->get<bool>() ? 1 : 0);
    else if (v->is_number_unsigned())
      put64(b, v->get<std::uint64_t>());
    else
      put64(b, static_cast<std::uint64_t>(v->get<std::int64_t>()));
  }
  return b;
}
std::uint32_t crc32(const Bytes& b) {
  static const std::array<std::uint32_t, 256> table = []() {
    std::array<std::uint32_t, 256> t{};
    for (std::uint32_t n = 0; n < 256; ++n) {
      auto c = n;
      for (int i = 0; i < 8; ++i) c = (c & 1) ? 0xedb88320U ^ (c >> 1) : c >> 1;
      t[n] = c;
    }
    return t;
  }();
  std::uint32_t crc = 0xffffffffU;
  for (auto x : b) crc = table[(crc ^ x) & 255] ^ (crc >> 8);
  return crc ^ 0xffffffffU;
}
class Zip {
  struct Entry {
    std::string name;
    std::uint32_t crc, size, offset, packed;
    std::uint16_t method;
  };
  std::ofstream stream;
  std::vector<Entry> entries;
  std::uint64_t position = 0;
  void write(const Bytes& b) {
    if (position + b.size() > zip_limit)
      throw std::invalid_argument("NPZ archive exceeds ZIP32 4 GiB limit; split output");
    stream.write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size()));
    position += b.size();
  }

 public:
  explicit Zip(const std::filesystem::path& path)
      : stream(path, std::ios::binary | std::ios::trunc) {
    stream.exceptions(std::ios::badbit | std::ios::failbit);
  }
  void add(const std::string& key, const Bytes& data, bool compress = false) {
    if (key.empty() || key.front() == '/' || key.front() == '\\' ||
        key.find('\0') != std::string::npos || key.find('\\') != std::string::npos)
      throw std::invalid_argument("Invalid NPZ array key");
    std::size_t start = 0;
    for (std::size_t i = 0; i <= key.size(); ++i)
      if (i == key.size() || key[i] == '/') {
        if (key.substr(start, i - start) == "..")
          throw std::invalid_argument("NPZ key contains directory traversal");
        start = i + 1;
      }
    unicode(key);
    auto name = key + ".npy";
    if (name.size() > 65535 || entries.size() >= 65535 || data.size() > zip_limit)
      throw std::invalid_argument("NPZ ZIP32 limits exceeded");
    Bytes packed;
    if (compress) {
      z_stream zs{};
      if (deflateInit2(&zs, 3, Z_DEFLATED, -MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK)
        throw std::runtime_error("Cannot initialize NPZ compression");
      packed.resize(deflateBound(&zs, data.size()));
      zs.next_in = const_cast<Bytef*>(data.data()); zs.avail_in = data.size();
      zs.next_out = packed.data(); zs.avail_out = packed.size();
      const int status = deflate(&zs, Z_FINISH);
      const auto size = zs.total_out;
      deflateEnd(&zs);
      if (status != Z_STREAM_END) throw std::runtime_error("NPZ compression failed");
      packed.resize(size);
    }
    const auto& payload = compress ? packed : data;
    Entry e{name, crc32(data), static_cast<std::uint32_t>(data.size()),
            static_cast<std::uint32_t>(position), static_cast<std::uint32_t>(payload.size()),
            static_cast<std::uint16_t>(compress ? 8 : 0)};
    Bytes h;
    put32(h, 0x04034b50);
    put16(h, 20);
    put16(h, 0x0800);
    put16(h, e.method);
    put16(h, 0);
    put16(h, 33);
    put32(h, e.crc);
    put32(h, e.packed);
    put32(h, e.size);
    put16(h, static_cast<std::uint16_t>(name.size()));
    put16(h, 0);
    h.insert(h.end(), name.begin(), name.end());
    write(h);
    write(payload);
    entries.push_back(std::move(e));
  }
  void finish() {
    auto start = static_cast<std::uint32_t>(position);
    for (const auto& e : entries) {
      Bytes h;
      put32(h, 0x02014b50);
      put16(h, 20);
      put16(h, 20);
      put16(h, 0x0800);
      put16(h, e.method);
      put16(h, 0);
      put16(h, 33);
      put32(h, e.crc);
      put32(h, e.packed);
      put32(h, e.size);
      put16(h, static_cast<std::uint16_t>(e.name.size()));
      put16(h, 0);
      put16(h, 0);
      put16(h, 0);
      put16(h, 0);
      put32(h, 0);
      put32(h, e.offset);
      h.insert(h.end(), e.name.begin(), e.name.end());
      write(h);
    }
    auto size = static_cast<std::uint32_t>(position - start);
    Bytes end;
    put32(end, 0x06054b50);
    put16(end, 0);
    put16(end, 0);
    put16(end, static_cast<std::uint16_t>(entries.size()));
    put16(end, static_cast<std::uint16_t>(entries.size()));
    put32(end, size);
    put32(end, start);
    put16(end, 0);
    write(end);
    stream.flush();
    stream.close();
  }
};
template <class Fill>
void atomic_write(const std::filesystem::path& path, Fill fill) {
  if (path.empty() || path.filename().empty())
    throw std::invalid_argument("NPZ destination must be filename");
  if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
  static std::atomic<std::uint64_t> sequence{0};
  auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
  auto temp = path;
  temp += ".tmp-" + std::to_string(stamp) + "-" + std::to_string(sequence++);
  try {
    {
      Zip zip(temp);
      fill(zip);
      zip.finish();
    }
    std::filesystem::rename(temp, path);
  } catch (...) {
    std::error_code ignored;
    std::filesystem::remove(temp, ignored);
    throw;
  }
}
}  // namespace
void write_npz(const std::filesystem::path& path, const std::map<std::string, NpyArray>& arrays, bool compress) {
  atomic_write(path, [&](Zip& z) {
    for (const auto& item : arrays) z.add(item.first, encode(item.second), compress);
  });
}
void write_npz(const std::filesystem::path& path, const Json& arrays, bool compress) {
  if (!arrays.is_object()) throw std::invalid_argument("write_npz expects object of named arrays");
  atomic_write(path, [&](Zip& z) {
    for (auto i = arrays.begin(); i != arrays.end(); ++i) z.add(i.key(), encode(i.value()), compress);
  });
}
NpyArray read_npz_array(const std::filesystem::path& path, const std::string& key) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("Cannot read profile archive: " + path.string());
  auto read = [&](std::size_t n) {
    Bytes b(n);
    if (!f.read(reinterpret_cast<char*>(b.data()), n)) throw std::runtime_error("Truncated NPZ archive");
    return b;
  };
  auto u = [](const Bytes& b, std::size_t p, int n) {
    std::uint64_t v=0; for (int j=0;j<n;++j) v |= std::uint64_t(b.at(p+j)) << (8*j); return v;
  };
  while (f.peek() != EOF) {
    auto h=read(4);
    if (u(h,0,4)!=0x04034b50) break;
    auto rest=read(26); h.insert(h.end(),rest.begin(),rest.end());
    const auto method=u(h,8,2), packed=u(h,18,4), size=u(h,22,4);
    if ((u(h,6,2)&9) || (method!=0 && method!=8)) throw std::runtime_error("Unsupported NPZ ZIP format");
    auto name=read(u(h,26,2)); f.seekg(u(h,28,2),std::ios::cur);
    if (std::string(name.begin(),name.end())!=key+".npy") { f.seekg(packed,std::ios::cur); continue; }
    auto bytes=read(packed);
    if (method==8) {
      Bytes unpacked(size); z_stream zs{};
      if (inflateInit2(&zs,-MAX_WBITS)!=Z_OK) throw std::runtime_error("NPZ inflate init failed");
      zs.next_in=bytes.data(); zs.avail_in=bytes.size(); zs.next_out=unpacked.data(); zs.avail_out=unpacked.size();
      auto status=inflate(&zs,Z_FINISH); auto total=zs.total_out; inflateEnd(&zs);
      if (status!=Z_STREAM_END || total!=size) throw std::runtime_error("NPZ decompression failed");
      bytes=std::move(unpacked);
    }
    if (bytes.size()!=size || crc32(bytes)!=u(h,14,4)) throw std::runtime_error("NPZ CRC mismatch");
    if (bytes.size()<10 || std::string(bytes.begin(),bytes.begin()+6)!=std::string("\x93NUMPY",6))
      throw std::runtime_error("Invalid NPY signature");
    const auto major=bytes[6];
    if (major!=1 && major!=2) throw std::runtime_error("Unsupported NPY version");
    std::size_t prefix=major==1?10:12, length=u(bytes,8,major==1?2:4), offset=prefix+length;
    if (offset>bytes.size()) throw std::runtime_error("Truncated NPY header");
    std::string header(bytes.begin()+prefix,bytes.begin()+offset);
    bool floats=header.find("'<f8'")!=std::string::npos;
    bool ints=header.find("'<i8'")!=std::string::npos;
    if ((!floats&&!ints) || header.find("'fortran_order': False")==std::string::npos)
      throw std::runtime_error("Expected native C-order float64/int64 NPY");
    std::smatch match;
    if (!std::regex_search(header,match,std::regex("'shape': \\(([^)]*)\\)"))) throw std::runtime_error("Invalid NPY shape");
    std::string dims=match[1]; std::replace(dims.begin(),dims.end(),',',' '); std::istringstream ss(dims);
    NpyArray a; std::size_t d; while(ss>>d) a.shape.push_back(d);
    const auto n=count(a.shape);
    if (n>(bytes.size()-offset)/8 || n*8!=bytes.size()-offset) throw std::runtime_error("NPY size mismatch");
    a.data.reserve(n);
    for(std::size_t i=0;i<n;++i) {
      auto bits=u(bytes,offset+8*i,8); double v;
      if(floats) std::memcpy(&v,&bits,8); else v=static_cast<double>(static_cast<std::int64_t>(bits));
      a.data.push_back(v);
    }
    return a;
  }
  throw std::runtime_error("Missing NPZ array: "+key);
}
}  // namespace cycle
