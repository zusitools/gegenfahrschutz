#ifndef GEGENFAHRSCHUTZ_FILECACHE_HPP_
#define GEGENFAHRSCHUTZ_FILECACHE_HPP_

#include <algorithm>
#include <string>
#include <unordered_map>

#include <boost/nowide/iostream.hpp>

#include "zusi_parser/utils.hpp"
#include "zusi_parser/zusi_parser.hpp"
#include "zusi_parser/zusi_types.hpp"

class FileCache {
public:
  Zusi* get_datei(const std::string_view zusi_pfad) {
    std::string key(zusi_pfad);
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);
    auto it = _cache.find(key);
    if (it != std::end(_cache)) {
      return it->second.get();
    }

    try {
      boost::nowide::cout << "Parsing " << zusi_pfad << "\n";
      return _cache.emplace(std::make_pair(key, zusixml::parseFile(zusixml::zusiPfadZuOsPfad(zusi_pfad, "")))).first->second.get();
    } catch (const std::exception& e) {
      boost::nowide::cerr << "Error parsing " << zusi_pfad << ": " << e.what() << "\n";
      _cache[key] = nullptr;
      return nullptr;
    }
  }

private:
  std::unordered_map<std::string, std::unique_ptr<Zusi>> _cache;
};

#endif  // GEGENFAHRSCHUTZ_FILECACHE_HPP_
