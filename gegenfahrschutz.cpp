#include "file_cache.hpp"

#include "zusi_parser/zusi_types.hpp"
#include "zusi_parser/zusi_parser.hpp"
#include "zusi_parser/utils.hpp"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/functional/hash.hpp>
#include <boost/program_options.hpp>
#include <boost/system/error_code.hpp>
#include <boost/nowide/args.hpp>
#include <boost/nowide/iostream.hpp>
#include <boost/nowide/fstream.hpp>

#include "rapidxml-1.13/rapidxml.hpp"
#include "rapidxml-1.13/rapidxml_print.hpp"

#include <algorithm>
#include <optional>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <set>
#include <string_view>
#include <vector>

#ifdef WIN32
#define UNICODE
#endif

namespace po = boost::program_options;

FileCache file_cache;

const ReferenzElement* get_refpunkt_by_nr(const Strecke& strecke, int nr) {
  if (nr < 0 || nr >= strecke.children_ReferenzElemente.size()) {
    return nullptr;
  }
  return strecke.children_ReferenzElemente[nr].get();
}

const StrElement* get_element_by_nr(const Strecke& strecke, int nr) {
  if (nr < 0 || nr >= strecke.children_StrElement.size()) {
    return nullptr;
  }
  return strecke.children_StrElement[nr].get();
}

struct ElementRichtungRef {
  const Strecke* st3;
  const StrElement* element;
  bool normrichtung;

  operator bool() const {
    return static_cast<bool>(element);
  }

  bool operator==(const ElementRichtungRef& other) {
    return (this->element == other.element) && (this->normrichtung == other.normrichtung);
  }

  const std::optional<StreckenelementRichtungsInfo>& richtungsInfo() const {
    return normrichtung ? element->InfoNormRichtung : element->InfoGegenRichtung;
  }

  ElementRichtungRef gegenrichtung() const {
    return { st3, element, !normrichtung };
  }

  const decltype(StrElement::children_NachNorm)& nachfolger_selbes_modul() const {
    return normrichtung ? element->children_NachNorm : element->children_NachGegen;
  }

  const decltype(StrElement::children_NachNormModul)& nachfolger_anderes_modul() const {
    return normrichtung ? element->children_NachNormModul : element->children_NachGegenModul;
  }
};

ElementRichtungRef get_element_by_ref_nr(const Strecke& strecke, int nr) {
  const auto& refpunkt = get_refpunkt_by_nr(strecke, nr);
  if (!refpunkt) {
    return { nullptr, nullptr, false };
  }
  return { &strecke, get_element_by_nr(strecke, refpunkt->StrElement), refpunkt->StrNorm };
}

// Eine geordnete Sequenz von Streckenelementen ohne Weichen.
// Damit ist es moeglich, zwei Streckenelemente bezueglich ihrer Position in der Sequenz zu vergleichen.
class ElementSeq {
 public:
  explicit ElementSeq(std::vector<ElementRichtungRef> elemente)
      : _elemente(std::move(elemente)) {
    for (size_t i = 0; i < _elemente.size(); ++i) {
      _indizes.emplace(_elemente[i].element, i);
    }
  }

  std::optional<size_t> index(const StrElement* element) const {
    const auto& it = _indizes.find(element);
    return (it == _indizes.end()) ? std::nullopt : std::optional(it->second);
  }

  const std::vector<ElementRichtungRef>& elemente() const {
    return _elemente;
  }

 private:
  std::vector<ElementRichtungRef> _elemente;
  std::unordered_map<const StrElement*, size_t> _indizes;
};

// Beschreibt eine Fahrstrasse, die ueber Elemente aus einer Elementsequenz verlaeuft.
// Da die Elementsequenz keine Weichen enthaelt, kann der Teil der Fahrstrasse, der innerhalb der Sequenz verlaeuft,
// durch das Intervall [`start_element_seq_idx`, `ende_element_seq_idx`] (jeweils einschliesslich) beschrieben werden.
// Register duerfen nur im Intervall [`min_register_element_seq_idx`, `ende_element_seq_idx`] (jeweils einschliesslich)
// gesetzt werden.
//
// Eine Fahrstrasse verlaeuft in *Normrichtung*, wenn start_element_seq_idx < ende_element_seq_idx,
// sonst in *Gegenrichtung*.
struct ElementSeqFahrstr {
  const Fahrstrasse* fahrstrasse;
  size_t start_element_seq_idx;
  size_t min_register_element_seq_idx;
  size_t ende_element_seq_idx;

  bool ist_normrichtung() const {
    return start_element_seq_idx < ende_element_seq_idx;
  }

  bool enthaelt_element_seq_idx(size_t element_seq_idx) const {
    const auto p = std::minmax(start_element_seq_idx, ende_element_seq_idx);
    return p.first <= element_seq_idx && element_seq_idx <= p.second;
  }
};

struct ModulInfo {
  std::string pfad_kurz;
  std::string pfad_os;
  std::string pfad_zusi;
};

std::unordered_map<const Strecke*, ModulInfo> modul_info;

std::string get_dateiname(const zusixml::ZusiPfad& zusi_pfad) {
  const auto& zusi_pfad_str = zusi_pfad.alsZusiPfad();
  const auto& pos = zusi_pfad_str.rfind(zusixml::zusiSep);
  return std::string((pos == std::string::npos || pos == zusi_pfad_str.size() - 1) ? zusi_pfad_str : zusi_pfad_str.substr(pos + 1));
}

const Strecke* get_strecke(const Dateiverknuepfung& datei) {
  // TODO: Wir gehen davon aus, dass die Dateiverknuepfung immer den kompletten Pfad enthaelt.
  // Das ist bei der in Zusi ueblichen Modulstruktur der Fall.
  const auto zusi_pfad = zusixml::ZusiPfad::vonZusiPfad(datei.Dateiname);
  const auto& zusi = file_cache.get_datei(zusi_pfad);
  if (!zusi || !zusi->Strecke) {
    return nullptr;
  }
  const auto* result = zusi->Strecke.get();
  const auto& it = modul_info.find(result);
  if (it == modul_info.end()) {
    modul_info.emplace(result, ModulInfo { get_dateiname(zusi_pfad),
        zusi_pfad.alsOsPfad(), std::string(zusi_pfad.alsZusiPfad()) });
  }
  return result;
}

const Strecke* get_strecke(const std::string& os_pfad) {
  const auto& zusi_pfad = zusixml::ZusiPfad::vonOsPfad(os_pfad);
  const auto& zusi = file_cache.get_datei(zusi_pfad);
  if (!zusi || !zusi->Strecke) {
    return nullptr;
  }
  const auto* result = zusi->Strecke.get();
  const auto& it = modul_info.find(result);
  if (it == modul_info.end()) {
    modul_info.emplace(result, ModulInfo { get_dateiname(zusi_pfad), os_pfad, std::string(zusi_pfad.alsZusiPfad()) });
  }
  return result;
}

struct NeuesRegister {
  decltype(StreckenelementRichtungsInfo::Reg) register_nr;
  decltype(ReferenzElement::ReferenzNr) verkn_refpunkt_nr;
  std::string verkn_refpunkt_modul;
};

struct NeuerRegisterRefpunkt {
  decltype(ReferenzElement::ReferenzNr) nr;
  decltype(StrElement::Nr) element_nr;
  bool element_normrichtung;
  int32_t register_nr;
};

struct StreckenAenderung {
  std::vector<NeuerRegisterRefpunkt> neue_refpunkte;
  std::unordered_map<decltype(StrElement::Nr), std::pair<NeuesRegister, NeuesRegister>> neue_register;
};

// Gibt zurueck, ob die gegebenen Fahrstrassen durch Register gegeneinander verriegelt sind.
// Es werden sowohl die in der Fahrstrasse verknuepften Register als auch die in `neue_register`
// angegebenen geplanten Registerpositionen beruecksichtigt.
bool gegeneinander_verriegelt(const ElementSeqFahrstr& fs1, const ElementSeqFahrstr& fs2, const std::vector<std::pair<size_t, size_t>>& neue_register)
{
  assert(fs1.ist_normrichtung() != fs2.ist_normrichtung());

  for (const auto& r : neue_register) {
    if (fs1.enthaelt_element_seq_idx(r.first) && fs2.enthaelt_element_seq_idx(r.second)) {
      return true;
    }
  }

  std::set<std::pair<const Strecke*, decltype(StreckenelementRichtungsInfo::Reg)>> register_fs1;
  for (const auto& registerVerkn : fs1.fahrstrasse->children_FahrstrRegister)
  {
    const auto st3 = get_strecke(registerVerkn->Datei);
    if (!st3) {
      continue;
    }
    const auto element = get_element_by_ref_nr(*st3, registerVerkn->Ref);
    if (!element || !element.richtungsInfo()) {
      continue;
    }
    register_fs1.insert(std::make_pair(st3, element.richtungsInfo()->Reg));
  }

  for (const auto& registerVerkn : fs2.fahrstrasse->children_FahrstrRegister)
  {
    const auto st3 = get_strecke(registerVerkn->Datei);
    if (!st3) {
      continue;
    }
    const auto element = get_element_by_ref_nr(*st3, registerVerkn->Ref);
    if (!element || !element.richtungsInfo()) {
      continue;
    }
    if (register_fs1.find(std::make_pair(st3, element.richtungsInfo()->Reg)) != register_fs1.end()) {
      return true;
    }
  }
  return false;
}

template<typename ValueType>
ValueType& get_pair_element(std::pair<ValueType, ValueType>& p, bool first_element) {
  return first_element ? p.first : p.second;
}

template<typename ValueType>
const ValueType& get_pair_element(const std::pair<ValueType, ValueType>& p, bool first_element) {
  return first_element ? p.first : p.second;
}

void print_visualisierung(
  const std::vector<ElementSeqFahrstr>& fahrstr_normrichtung,
  const std::vector<ElementSeqFahrstr>& fahrstr_gegenrichtung) {

  std::set<size_t> indices;
  for (const auto& fss : { fahrstr_normrichtung, fahrstr_gegenrichtung }) {
    for (const auto& fs : fss) {
      indices.emplace(fs.start_element_seq_idx);
      indices.emplace(fs.ende_element_seq_idx);
      indices.emplace(fs.min_register_element_seq_idx);
    }
  }

  std::unordered_map<size_t, size_t> start_pos;
  std::unordered_map<size_t, size_t> ende_pos;
  std::string index_str;
  std::optional<size_t> prev = std::nullopt;
  for (const auto& idx : indices) {
    if (prev.has_value()) {
      if (*prev == idx - 1) {
        index_str += '.';
      } else {
        index_str += "...";
      }
    }
    prev = idx;

    start_pos.emplace(idx, index_str.size());
    index_str += std::to_string(idx);
    ende_pos.emplace(idx, index_str.size() - 1);
  }

  for (const auto& fs : fahrstr_normrichtung) {
    size_t pos = 0;
    while (pos < start_pos[fs.start_element_seq_idx]) { boost::nowide::cerr << " "; ++pos; }
    boost::nowide::cerr << '|'; ++pos;
    while (pos < start_pos[fs.min_register_element_seq_idx]) { boost::nowide::cerr << "-"; ++pos; }
    while (pos + 1 < ende_pos[fs.ende_element_seq_idx]) { boost::nowide::cerr << "="; ++pos; }
    boost::nowide::cerr << ">| " << fs.fahrstrasse->FahrstrName << std::endl;
  }

  boost::nowide::cerr << index_str << std::endl;

  for (const auto& fs : fahrstr_gegenrichtung) {
    size_t pos = 0;
    while (pos < start_pos[fs.ende_element_seq_idx]) { boost::nowide::cerr << " "; ++pos; }
    boost::nowide::cerr << "|<"; pos += 2;
    while (pos < ende_pos[fs.min_register_element_seq_idx]) { boost::nowide::cerr << '='; ++pos; }
    while (pos < ende_pos[fs.start_element_seq_idx]) { boost::nowide::cerr << '-'; ++pos; }
    boost::nowide::cerr << "| " << fs.fahrstrasse->FahrstrName << std::endl;
  }
}

int main(int argc, char** argv) {
  boost::nowide::args a(argc, argv);

  std::string start_st3_pfad;
  int start_element_nr;
  std::string ziel_st3_pfad;
  int ziel_element_nr;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("debug", "Debug-Ausgaben aktivieren")
    ("start-st3", po::value<std::string>(&start_st3_pfad), "ST3-Datei des Start-Streckenelementes")
    ("start-element", po::value<int>(&start_element_nr), "Nummer des Start-Streckenelementes")
    ("ziel-st3", po::value<std::string>(&ziel_st3_pfad), "ST3-Datei des Ziel-Streckenelementes")
    ("ziel-element", po::value<int>(&ziel_element_nr), "Nummer des Ziel-Streckenelementes")
    ;

  po::positional_options_description p;
  p.add("start-st3", 1);
  p.add("start-element", 1);
  p.add("ziel-st3", 1);
  p.add("ziel-element", 1);

  po::variables_map vars;
  po::store(po::command_line_parser(argc, argv). options(desc).positional(p).run(), vars);
  po::notify(vars);

  if (vars.count("help")) {
    boost::nowide::cout << desc << std::endl;
    return 1;
  }

  const bool debug = vars.count("debug");

  const auto& start_st3 = get_strecke(start_st3_pfad);
  if (!start_st3)
  {
    boost::nowide::cerr << "Fehler beim Laden von " << start_st3_pfad << std::endl;
    return 1;
  }

  // Betrachte auch Fahrstrassen, die in Nachbarmodulen des Start- und Zielmoduls beginnen und enden
  for (const auto& modulDatei : start_st3->children_ModulDateien) {
    get_strecke(modulDatei->Datei);
  }

  const auto& ziel_st3 = get_strecke(ziel_st3_pfad);
  if (!ziel_st3)
  {
    boost::nowide::cerr << "Fehler beim Laden von " << ziel_st3_pfad << std::endl;
    return 1;
  }

  if (start_st3 != ziel_st3) {
    for (const auto& modulDatei : ziel_st3->children_ModulDateien) {
      get_strecke(modulDatei->Datei);
    }
  }

  const StrElement* start_element = get_element_by_nr(*start_st3, start_element_nr);
  if (!start_element)
  {
    boost::nowide::cerr << "Element " << start_element_nr << " existiert nicht in Datei " << start_st3_pfad << std::endl;
    return 1;
  }

  const StrElement* ziel_element = get_element_by_nr(*ziel_st3, ziel_element_nr);
  if (!ziel_element)
  {
    boost::nowide::cerr << "Element " << ziel_element_nr << " existiert nicht in Datei " << ziel_st3_pfad << std::endl;
    return 1;
  }

  std::vector<ElementRichtungRef> element_seq_elemente;

  // Finde Weg von Start nach Ziel ohne Weichen
  auto findeWegZuZiel = [&](bool start_normrichtung) -> bool {
    element_seq_elemente.clear();

    ElementRichtungRef cur_element { start_st3, start_element, start_normrichtung };
    constexpr size_t max_num_elemente = 10000;
    for (size_t idx = 0; idx < max_num_elemente; ++idx) {
      element_seq_elemente.emplace_back(cur_element);

      if (cur_element.element == ziel_element) {
        boost::nowide::cout << " - Fahrwegsuche in " << (start_normrichtung ? "Norm" : "Gegen") << "richtung erfolgreich" << std::endl;
        return true;
      }

      if (idx > 0 && (cur_element.gegenrichtung().nachfolger_selbes_modul().size() + cur_element.gegenrichtung().nachfolger_anderes_modul().size() > 1)) {
        // Weiche in der Gegenrichtung
        boost::nowide::cout << " - Fahrwegsuche in " << (start_normrichtung ? "Norm" : "Gegen") << "richtung: ignoriere stumpf befahrene Weiche an " << modul_info.at(cur_element.st3).pfad_kurz << ", Element " << cur_element.element->Nr << std::endl;
      }

      if (cur_element.nachfolger_selbes_modul().size() + cur_element.nachfolger_anderes_modul().size() == 0) {
        // Streckenende
        boost::nowide::cout << " - Fahrwegsuche in " << (start_normrichtung ? "Norm" : "Gegen") << "richtung abgebrochen an " << modul_info.at(cur_element.st3).pfad_kurz << ", Element " << cur_element.element->Nr << " (kein Nachfolger)" << std::endl;
        return false;
      } else if (cur_element.nachfolger_selbes_modul().size() + cur_element.nachfolger_anderes_modul().size() > 1) {
        boost::nowide::cout << " - Fahrwegsuche in " << (start_normrichtung ? "Norm" : "Gegen") << "richtung: ignoriere spitz befahrene Weiche an " << modul_info.at(cur_element.st3).pfad_kurz << ", Element " << cur_element.element->Nr << " (nimm ersten Nachfolger)" << std::endl;
      }

      const auto& next_element = [&]() {
        if (!cur_element.nachfolger_selbes_modul().empty()) {
          const int32_t anschluss_mask = (cur_element.normrichtung ? 0x1 : 0x100);
          return ElementRichtungRef { cur_element.st3, get_element_by_nr(*cur_element.st3, cur_element.nachfolger_selbes_modul()[0].Nr), (cur_element.element->Anschluss & anschluss_mask) == 0 };
        } else {
          assert(!cur_element.nachfolger_anderes_modul().empty());
          const auto& st3 = get_strecke(cur_element.nachfolger_anderes_modul()[0].Datei);
          if (!st3) {
            return ElementRichtungRef {};
          }
          return get_element_by_ref_nr(*st3, cur_element.nachfolger_anderes_modul()[0].Nr).gegenrichtung(); // Referenzpunkt zeigt in Richtung Modulgrenze
        }
      }();
      if (!next_element) {
        boost::nowide::cout << " - Fahrwegsuche in " << (start_normrichtung ? "Norm" : "Gegen") << "richtung abgebrochen an " << modul_info.at(cur_element.st3).pfad_kurz << ", Element " << cur_element.element->Nr << " (kein Nachfolger)\n";
        return false;
      }
      cur_element = next_element;
    }
    assert(cur_element);
    boost::nowide::cout << " - Fahrwegsuche in " << (start_normrichtung ? "Norm" : "Gegen") << "richtung abgebrochen an " << modul_info.at(cur_element.st3).pfad_kurz << ", Element " << cur_element.element->Nr << " (Zielelement nach " << max_num_elemente << " Schritten nicht erreicht)" << std::endl;

    return false;
  };

  if (!findeWegZuZiel(true) && !findeWegZuZiel(false)) {
    boost::nowide::cerr << "Keinen Weg zwischen den angegebenen Streckenelementen gefunden, der an Weichen immer der Vorzugslage folgt.\n";
    return 1;
  }

  const ElementSeq element_seq(std::move(element_seq_elemente));

  std::vector<ElementSeqFahrstr> fahrstr_normrichtung;
  std::vector<ElementSeqFahrstr> fahrstr_gegenrichtung;

  // Pruefe alle Fahrstrassen in allen Modulen auf Ueberschneidung mit `element_seq`
  // und speichere sie in `fahrstr_normrichtung` bzw. `fahrstr_gegenrichtung`.
  std::vector<const Strecke*> module;
  for (const auto& p : modul_info) {
    module.push_back(p.first);
  }
  for (const auto& modul : module) {
    for (const auto& fahrstrasse : modul->children_Fahrstrasse) {
      if (!fahrstrasse->FahrstrStart || !fahrstrasse->FahrstrZiel) {
        continue;
      }

      const auto& start_st3 = get_strecke(fahrstrasse->FahrstrStart->Datei);
      if (!start_st3) {
        continue;
      }
      const auto& start_element = get_element_by_ref_nr(*start_st3, fahrstrasse->FahrstrStart->Ref);
      const auto& ziel_st3 = get_strecke(fahrstrasse->FahrstrZiel->Datei);
      if (!ziel_st3) {
        continue;
      }
      const auto& ziel_element = get_element_by_ref_nr(*ziel_st3, fahrstrasse->FahrstrZiel->Ref);

      ElementRichtungRef cur_element = start_element;
      size_t idx = 0;
      size_t fahrstr_weichen_idx = 0;
      size_t fahrstr_aufloesepunkt_idx = 0;

      std::optional<size_t> start_element_seq_idx = std::nullopt;
      std::optional<size_t> letzter_element_seq_idx = std::nullopt;
      std::optional<size_t> letzter_aufloesepunkt_element_seq_idx = std::nullopt;

      while (cur_element && (idx <= 9999)) {
        const auto element_seq_idx = element_seq.index(cur_element.element);
        if ((idx != 0) && element_seq_idx.has_value()) {  // idx != 0: Startelement ist selbst nicht in der Fahrstrasse enthalten
          if (!start_element_seq_idx.has_value()) {
            start_element_seq_idx = element_seq_idx;
          }
          letzter_element_seq_idx = element_seq_idx;
        }

        idx++;
        if (cur_element == ziel_element) {
          if (start_element_seq_idx.has_value()) {
            assert(letzter_element_seq_idx.has_value());
            const auto min_register_element_seq_idx = [&]() {
              if (!letzter_aufloesepunkt_element_seq_idx.has_value()) {
                return *start_element_seq_idx;
              } else if (*letzter_aufloesepunkt_element_seq_idx < *letzter_element_seq_idx) {
                // Ein Register im Aufloeseelement selbst wird auch freigegeben,
                // deshalb darf das Register zum Verriegeln erst ein Element spaeter eingefuegt werden.
                return *letzter_aufloesepunkt_element_seq_idx + 1;
              } else {
                assert(*letzter_aufloesepunkt_element_seq_idx > *letzter_element_seq_idx);
                return *letzter_aufloesepunkt_element_seq_idx - 1;
              }
            }();
            ElementSeqFahrstr seq_fahrstr { fahrstrasse.get(), *start_element_seq_idx,
              min_register_element_seq_idx, *letzter_element_seq_idx };
            (seq_fahrstr.ist_normrichtung() ? fahrstr_normrichtung : fahrstr_gegenrichtung).push_back(std::move(seq_fahrstr));
          }
          break;
        }

        if (element_seq_idx.has_value()) {
          // Pruefe, ob am aktuellen Element ein Teilaufloesepunkt liegt.
          // (Ein Teilaufloesepunkt kann nicht im Zielelement liegen,
          // dort kann nur ein Aufloesepunkt liegen.)
          for (size_t i = 0; i < fahrstrasse->children_FahrstrTeilaufloesung.size(); ++i)
          {
            // Optimierung: Normalerweise sind die Teilaufloesepunkte in der korrekten Reihenfolge in der Fahrstrasse verknuepft.
            const auto aufloesepunkte_idx = (fahrstr_aufloesepunkt_idx + i) % fahrstrasse->children_FahrstrTeilaufloesung.size();
            const auto& aufloesepunkt = fahrstrasse->children_FahrstrTeilaufloesung[aufloesepunkte_idx];
            const auto st3 = get_strecke(aufloesepunkt->Datei);
            if (!st3) {
              break;
            }
            const auto refpunkt = get_refpunkt_by_nr(*st3, aufloesepunkt->Ref);
            if (!refpunkt) {
              break;
            }
            if ((refpunkt->StrElement == cur_element.element->Nr) && (refpunkt->StrNorm == cur_element.normrichtung)) {
              letzter_aufloesepunkt_element_seq_idx = element_seq_idx;
              fahrstr_aufloesepunkt_idx = aufloesepunkte_idx + 1;
              break;
            }
          }
        }

        if (cur_element.nachfolger_selbes_modul().size() + cur_element.nachfolger_anderes_modul().size() == 0) {
          break;
        }

        size_t nachfolgerIdx = 0;
        if (cur_element.nachfolger_selbes_modul().size() + cur_element.nachfolger_anderes_modul().size() > 1) {
          // Weiche: Bestimme Weichenstellung aus Fahrstrasse
          bool found = false;

          for (size_t i = 0; i < fahrstrasse->children_FahrstrWeiche.size(); ++i)
          {
            // Optimierung: Normalerweise sind die Weichen in der korrekten Reihenfolge in der Fahrstrasse verknuepft.
            const auto weichenliste_idx = (fahrstr_weichen_idx + i) % fahrstrasse->children_FahrstrWeiche.size();
            const auto& weiche = fahrstrasse->children_FahrstrWeiche[weichenliste_idx];
            const auto st3 = get_strecke(weiche->Datei);
            if (!st3) {
              break;
            }
            const auto refpunkt = get_refpunkt_by_nr(*st3, weiche->Ref);
            if (!refpunkt) {
              break;
            }
            if ((refpunkt->StrElement == cur_element.element->Nr) && (refpunkt->StrNorm == cur_element.normrichtung)) {
              nachfolgerIdx = weiche->FahrstrWeichenlage - 1;
              fahrstr_weichen_idx = weichenliste_idx + 1;
              found = true;
              break;
            }
          }

          if (!found) {
            break;
          }
        }

        // TODO: Momentan kann ein Streckenelement entweder Nachfolger im Modul oder ausserhalb des Moduls haben, nicht beides.
        if (nachfolgerIdx < cur_element.nachfolger_selbes_modul().size()) {
          const int32_t anschluss_mask = (cur_element.normrichtung ? 0x1 : 0x100) << nachfolgerIdx;
          cur_element = { cur_element.st3, get_element_by_nr(*cur_element.st3, cur_element.nachfolger_selbes_modul()[nachfolgerIdx].Nr), (cur_element.element->Anschluss & anschluss_mask) == 0 };
        } else {
          nachfolgerIdx -= cur_element.nachfolger_selbes_modul().size();
          if (nachfolgerIdx >= cur_element.nachfolger_anderes_modul().size()) {
            break;
          }
          const auto st3 = get_strecke(cur_element.nachfolger_anderes_modul()[nachfolgerIdx].Datei);
          if (!st3) {
            break;
          }
          cur_element = get_element_by_ref_nr(*st3, cur_element.nachfolger_anderes_modul()[nachfolgerIdx].Nr).gegenrichtung();
        }
      }
    }
  }

  // Gibt fuer jedes Element aus `element_seq` an, ob es fuer ein neues Register zur Verfuegung steht,
  // also kein Register und kein Ereignis "Register verknuepfen" enthaelt.
  // Einmal fuer die Laufrichtung von `element_seq`, einmal fuer die Gegenrichtung.
  std::vector<bool> freie_elemente_laufrichtung(element_seq.elemente().size()), freie_elemente_gegenrichtung(element_seq.elemente().size());
  for (size_t i = 0; i < element_seq.elemente().size(); ++i) {
    const auto ist_frei = [](const std::optional<StreckenelementRichtungsInfo>& richtungsInfo) -> bool {
      if (!richtungsInfo) {
        return true;
      }
      if (richtungsInfo->Reg != 0) {
        return false;
      }
      for (const auto& ereignis : richtungsInfo->children_Ereignis) {
        if (ereignis->Er == 34 || ereignis->Er == 35) {
          return false;
        }
      }
      return true;
    };
    const auto& element_ref = element_seq.elemente()[i];
    (element_ref.normrichtung ? freie_elemente_laufrichtung : freie_elemente_gegenrichtung)[i] = ist_frei(element_ref.element->InfoNormRichtung);
    (element_ref.normrichtung ? freie_elemente_gegenrichtung : freie_elemente_laufrichtung)[i] = ist_frei(element_ref.element->InfoGegenRichtung);
  }

  // Jedes neue Register verriegelt zwei Elemente aus `element_seq` gegeneinander.
  std::vector<std::pair<size_t, size_t>> neue_register;
  std::unordered_set<const Strecke*> geaenderte_module;

  // Sortiere Fahrstrassen so, dass kuerzere Fahrstrassen (Aufgleisfahrstrassen)
  // vor laengeren Fahrstrassen mit demselben Ziel bearbeitet werden.
  // So ergeben sich mehr gemeinsame Nutzungen von Registern.
  std::sort(fahrstr_normrichtung.begin(), fahrstr_normrichtung.end(), [](const auto& fs1, const auto& fs2) {
      return fs1.start_element_seq_idx > fs2.start_element_seq_idx;
    });
  std::sort(fahrstr_gegenrichtung.begin(), fahrstr_gegenrichtung.end(), [](const auto& fs1, const auto& fs2) {
      return fs1.start_element_seq_idx < fs2.start_element_seq_idx;
    });

  if (debug) {
    boost::nowide::cout << "Fahrstrassen in Normrichtung:" << std::endl;
    for (const auto& seq_fahrstr : fahrstr_normrichtung) {
      boost::nowide::cerr << " - " << seq_fahrstr.fahrstrasse->FahrstrName << " [" << seq_fahrstr.start_element_seq_idx << ", " << seq_fahrstr.ende_element_seq_idx << "] [" << seq_fahrstr.min_register_element_seq_idx << ", " << seq_fahrstr.ende_element_seq_idx << "]" << std::endl;
    }

    boost::nowide::cerr << "Fahrstrassen in Gegenrichtung:" << std::endl;
    for (const auto& seq_fahrstr : fahrstr_gegenrichtung) {
      boost::nowide::cerr << " - " << seq_fahrstr.fahrstrasse->FahrstrName << " [" << seq_fahrstr.start_element_seq_idx << ", " << seq_fahrstr.ende_element_seq_idx << "] [" << seq_fahrstr.min_register_element_seq_idx << ", " << seq_fahrstr.ende_element_seq_idx << "]" << std::endl;
    }

    print_visualisierung(fahrstr_normrichtung, fahrstr_gegenrichtung);
  }

  for (const auto& seq_fahrstr_norm : fahrstr_normrichtung) {
    for (const auto& seq_fahrstr_gegen : fahrstr_gegenrichtung) {
      // Normrichtung: start <= min_register < ende
      assert(seq_fahrstr_norm.start_element_seq_idx <= seq_fahrstr_norm.min_register_element_seq_idx);
      assert(seq_fahrstr_norm.min_register_element_seq_idx < seq_fahrstr_norm.ende_element_seq_idx);
      // Gegenrichtung: ende < min_register <= start
      assert(seq_fahrstr_gegen.ende_element_seq_idx < seq_fahrstr_gegen.min_register_element_seq_idx);
      assert(seq_fahrstr_gegen.min_register_element_seq_idx <= seq_fahrstr_gegen.start_element_seq_idx);

      if (seq_fahrstr_gegen.start_element_seq_idx < seq_fahrstr_norm.start_element_seq_idx) {
        // Die Fahrstrassen muessen nicht gegeneinander verriegelt werden,
        // weil eine Fahrstrasse "im Ruecken" der anderen liegt.
        if (debug) {
          boost::nowide::cerr << " + Die Fahrstraßen \"" << seq_fahrstr_norm.fahrstrasse->FahrstrName << "\" und \"" << seq_fahrstr_gegen.fahrstrasse->FahrstrName << "\" muessen nicht gegeneinander verriegelt werden." << std::endl;
        }
        continue;
      }

      if (gegeneinander_verriegelt(seq_fahrstr_norm, seq_fahrstr_gegen, neue_register)) {
        // Die Fahrstrassen sind bereits gegeneinander verriegelt -- pruefe aber noch,
        // ob das auch ohne die neuen Register der Fall waere. Falls nein, muessen die
        // Fahrstrassen in den betroffenen Modulen neu erzeugt werden.
        if (!gegeneinander_verriegelt(seq_fahrstr_norm, seq_fahrstr_gegen, {})) {
          if (debug) {
            boost::nowide::cerr << " + Die Fahrstraßen \"" << seq_fahrstr_norm.fahrstrasse->FahrstrName << "\" und \"" << seq_fahrstr_gegen.fahrstrasse->FahrstrName << "\" sind durch neue Register gegeneinander verriegelt" << std::endl;
          }
          geaenderte_module.insert(get_strecke(seq_fahrstr_norm.fahrstrasse->FahrstrStart->Datei));
          geaenderte_module.insert(get_strecke(seq_fahrstr_gegen.fahrstrasse->FahrstrStart->Datei));
        } else {
          if (debug) {
            boost::nowide::cerr << " + Die Fahrstraßen \"" << seq_fahrstr_norm.fahrstrasse->FahrstrName << "\" und \"" << seq_fahrstr_gegen.fahrstrasse->FahrstrName << "\" sind bereits gegeneinander verriegelt, auch ohne neue Register" << std::endl;
          }
        }
        continue;
      }

      if (debug) {
        boost::nowide::cerr << " + Verriegele die Fahrstraßen \"" << seq_fahrstr_norm.fahrstrasse->FahrstrName << "\" und \"" << seq_fahrstr_gegen.fahrstrasse->FahrstrName << "\" gegeneinander" << std::endl;
      }

      size_t idx_norm = seq_fahrstr_norm.min_register_element_seq_idx;
      while (!freie_elemente_laufrichtung[idx_norm] && (idx_norm <= seq_fahrstr_norm.ende_element_seq_idx)) {
        ++idx_norm;
      }

      size_t idx_gegen_plus_1 = seq_fahrstr_gegen.min_register_element_seq_idx + 1;  // + 1, damit es kleiner werden kann als seq_fahrstr_gegen.ende_element_seq_idx + 1, das hoechstens 1 ist.
      while (!freie_elemente_gegenrichtung[idx_gegen_plus_1 - 1] && (idx_gegen_plus_1 >= seq_fahrstr_gegen.ende_element_seq_idx + 1)) {
        --idx_gegen_plus_1;
      }

      if (idx_norm > seq_fahrstr_norm.ende_element_seq_idx) {
        const auto& el1 = element_seq.elemente().at(seq_fahrstr_norm.min_register_element_seq_idx);
        const auto& modul1 = modul_info.at(el1.st3);
        const auto& el2 = element_seq.elemente().at(seq_fahrstr_norm.ende_element_seq_idx);
        const auto& modul2 = modul_info.at(el2.st3);
        boost::nowide::cerr << "Kein Element ohne bestehendes oder neues Register zwischen Index "
           << seq_fahrstr_norm.min_register_element_seq_idx << " und " << seq_fahrstr_norm.ende_element_seq_idx << " - "
           << modul1.pfad_kurz << ", Element " << el1.element->Nr << " " << (el1.normrichtung ? "blau" : "grün") << " und "
           << modul2.pfad_kurz << ", Element " << el2.element->Nr << " " << (el2.normrichtung ? "blau" : "grün") << " gefunden, um die Fahrstraßen \""
           << seq_fahrstr_norm.fahrstrasse->FahrstrName << "\" und \"" << seq_fahrstr_gegen.fahrstrasse->FahrstrName << "\" gegeneinander zu verriegeln.\n";
        return 1;
      }
      if (idx_gegen_plus_1 < seq_fahrstr_gegen.ende_element_seq_idx + 1) {
        const auto& el1 = element_seq.elemente().at(seq_fahrstr_gegen.min_register_element_seq_idx);
        const auto& modul1 = modul_info.at(el1.st3);
        const auto& el2 = element_seq.elemente().at(seq_fahrstr_gegen.ende_element_seq_idx);
        const auto& modul2 = modul_info.at(el2.st3);
        boost::nowide::cerr << "Kein Element ohne bestehendes oder neues Register zwischen Index "
           << seq_fahrstr_gegen.min_register_element_seq_idx << " und " << seq_fahrstr_gegen.ende_element_seq_idx << " - "
           << modul1.pfad_kurz << ", Element " << el1.element->Nr << " " << (el1.normrichtung ? "blau" : "grün") << " und "
           << modul2.pfad_kurz << ", Element " << el2.element->Nr << " " << (el2.normrichtung ? "blau" : "grün") << " gefunden, um die Fahrstraßen \""
           << seq_fahrstr_norm.fahrstrasse->FahrstrName << "\" und \"" << seq_fahrstr_gegen.fahrstrasse->FahrstrName << "\" gegeneinander zu verriegeln\n";
        return 1;
      }

      if (debug) {
        boost::nowide::cerr << " + Neue Register an Index " << idx_norm << " und " << (idx_gegen_plus_1 - 1) << std::endl;
      }
      neue_register.emplace_back(idx_norm, idx_gegen_plus_1 - 1);
      geaenderte_module.insert(get_strecke(seq_fahrstr_norm.fahrstrasse->FahrstrStart->Datei));
      geaenderte_module.insert(get_strecke(seq_fahrstr_gegen.fahrstrasse->FahrstrStart->Datei));
      freie_elemente_laufrichtung[idx_norm] = false;
      freie_elemente_gegenrichtung[idx_gegen_plus_1 - 1] = false;
    }
  }

  // Umsetzung Verriegelung zweier Elemente:
  // Wenn in derselben Moduldatei: neues manuelles Register in beiden Elementen (jeweils eine Richtung)
  //   + entsprechende Referenzpunkte
  // Wenn nicht in derselben Moduldatei: neue manuelle Register in beiden Elementen (jeweils eine Richtung)
  //   + Referenzpunkt fuer das Register + Register in Fahrstr verknuepfen "A" des jeweils anderen Registers
  std::unordered_map<const Strecke*, decltype(StreckenelementRichtungsInfo::Reg)> maxManuellesRegister;
  const auto neuesManuellesRegister = [&maxManuellesRegister](const Strecke* strecke) -> decltype(StreckenelementRichtungsInfo::Reg) {
    const auto& it = maxManuellesRegister.find(strecke);
    const auto result = 1 + (it != maxManuellesRegister.end() ? it->second : [&strecke]() {
      int32_t result = 0;
      for (const auto& strElement : strecke->children_StrElement) {
        if (!strElement) {
          continue;
        }
        if (strElement->InfoNormRichtung.has_value() && (strElement->InfoNormRichtung->Reg < 5000)) {
          result = std::max(result, strElement->InfoNormRichtung->Reg);
        }
        if (strElement->InfoGegenRichtung.has_value() && (strElement->InfoGegenRichtung->Reg < 5000)) {
          result = std::max(result, strElement->InfoGegenRichtung->Reg);
        }
      }
      return result;
    }());
    maxManuellesRegister.insert_or_assign(strecke, result);
    return result;
  };

  std::unordered_map<const Strecke*, decltype(ReferenzElement::ReferenzNr)> maxRefpunkt;
  const auto neuerRefpunkt = [&maxRefpunkt](const Strecke* strecke) -> decltype(ReferenzElement::ReferenzNr) {
    const auto& it = maxRefpunkt.find(strecke);
    const auto result = it != maxRefpunkt.end() ? it->second + 1 : strecke->children_ReferenzElemente.size();
    maxRefpunkt.insert_or_assign(strecke, result);
    return result;
  };

  std::unordered_map<const Strecke*, StreckenAenderung> strecken_aenderungen;

  for (const auto& r : neue_register) {
    const auto& el1 = element_seq.elemente().at(r.first);
    const auto& el2 = element_seq.elemente().at(r.second);
    const auto& modul1 = modul_info.at(el1.st3);
    const auto& modul2 = modul_info.at(el2.st3);

    const auto neues_register = neuesManuellesRegister(el1.st3);
    const auto refpunkt_el1 = neuerRefpunkt(el1.st3);
    get_pair_element(strecken_aenderungen[el1.st3].neue_register[el1.element->Nr],  el1.normrichtung) = NeuesRegister { neues_register, 0, "" };
    strecken_aenderungen[el1.st3].neue_refpunkte.push_back(NeuerRegisterRefpunkt { refpunkt_el1, el1.element->Nr, el1.normrichtung, neues_register });

    if (el1.st3 == el2.st3) {
      const auto refpunkt_el2 = neuerRefpunkt(el2.st3);
      boost::nowide::cout << " - Neues manuelles Register " << neues_register << " an " << modul1.pfad_kurz << ", Element " << el1.element->Nr << " " << (el1.normrichtung ? "blau" : "grün") << " und " << el2.element->Nr << " " << (!el2.normrichtung ? "blau" : "grün") << std::endl;
      boost::nowide::cout << "   - Neuer Referenzpunkt " << refpunkt_el1 << " (Typ Register) an " << modul1.pfad_kurz << ", Element " << el1.element->Nr << " " << (el1.normrichtung ? "blau" : "grün") << std::endl;
      boost::nowide::cout << "   - Neuer Referenzpunkt " << refpunkt_el2 << " (Typ Register) an " << modul2.pfad_kurz << ", Element " << el2.element->Nr << " " << (!el2.normrichtung ? "blau" : "grün") << std::endl;
      get_pair_element(strecken_aenderungen[el2.st3].neue_register[el2.element->Nr], !el2.normrichtung) = NeuesRegister { neues_register, 0, "" };
      strecken_aenderungen[el2.st3].neue_refpunkte.push_back(NeuerRegisterRefpunkt { refpunkt_el2, el2.element->Nr, !el2.normrichtung, neues_register });
    } else {
      boost::nowide::cout << " - Neues manuelles Register " << neues_register << " an " << modul1.pfad_kurz << ", Element " << el1.element->Nr << " " << (el1.normrichtung ? "blau" : "grün") << std::endl;
      boost::nowide::cout << "   - Neuer Referenzpunkt " << refpunkt_el1 << " (Typ Register) an " << modul1.pfad_kurz << ", Element " << el1.element->Nr << " " << (el1.normrichtung ? "blau" : "grün") << std::endl;
      boost::nowide::cout << "   - Neues Ereignis \"Register in Fahrstrasse verknuepfen A\" mit Referenznummer " << refpunkt_el1 << " und Moduldatei " << modul1.pfad_zusi << " an " << modul2.pfad_kurz << ", Element " << el2.element->Nr << " " << (!el2.normrichtung ? "blau" : "grün") << std::endl;
      get_pair_element(strecken_aenderungen[el2.st3].neue_register[el2.element->Nr], !el2.normrichtung) = NeuesRegister { 0, refpunkt_el1, modul1.pfad_zusi };
    }
  }

  if (strecken_aenderungen.empty()) {
    return 0;
  }

  boost::nowide::cout << "Die folgenden Dateien wurden neu geschrieben:" << std::endl;
  for (const auto& [strecke, aenderung] : strecken_aenderungen) {
    const auto& dateiname = modul_info.at(strecke).pfad_os;
    zusixml::FileReader reader(dateiname);
    rapidxml::xml_document<> doc;
    doc.parse<rapidxml::parse_non_destructive | rapidxml::parse_comment_nodes | rapidxml::parse_declaration_node>(const_cast<char*>(reader.data()));

    auto* const zusi_node = doc.first_node("Zusi");
    auto* const strecke_node = zusi_node->first_node("Strecke");

    // Neue Register und Registerverknuepfungen
    for (auto* str_element_node = strecke_node->first_node("StrElement"); str_element_node; str_element_node = str_element_node->next_sibling("StrElement")) {
      const auto* const nr_attrib = str_element_node->first_attribute("Nr");
      if (!nr_attrib) {
        continue;
      }

      const std::string s { nr_attrib->value(), nr_attrib->value_size() };
      int nr = atoi(s.c_str());

      const auto& it = aenderung.neue_register.find(nr);
      if (it == aenderung.neue_register.end()) {
        continue;
      }

      for (bool normrichtung : { true, false }) {
        const auto& neues_register = get_pair_element(it->second, normrichtung);
        if ((neues_register.register_nr == 0) && (neues_register.verkn_refpunkt_nr == 0)) {
          continue;
        }
        const auto& richtungs_info_tag = (normrichtung ? "InfoNormRichtung" : "InfoGegenRichtung");
        auto* richtungs_info_node = str_element_node->first_node(richtungs_info_tag);
        if (!richtungs_info_node) {
          richtungs_info_node = doc.allocate_node(rapidxml::node_element, richtungs_info_tag);
          str_element_node->append_node(richtungs_info_node);
        }
        if (neues_register.register_nr != 0) {
          // Zusi: <InfoGegenRichtung vMax="single" km="single" pos="bool" Reg="integer"
          auto* insert_before = richtungs_info_node->first_attribute();
          while (insert_before) {
            const auto& attr_name = std::string_view(insert_before->name(), insert_before->name_size());
            if ((attr_name != "vMax") && (attr_name != "km") && (attr_name != "pos")) {
              break;
            }
            insert_before = insert_before->next_attribute();
          }
          richtungs_info_node->insert_attribute(insert_before, doc.allocate_attribute("Reg", doc.allocate_string(std::to_string(neues_register.register_nr).c_str())));
        }
        if (neues_register.verkn_refpunkt_nr != 0) {
          auto* ereignis_node = doc.allocate_node(rapidxml::node_element, "Ereignis");
          richtungs_info_node->insert_node(richtungs_info_node->first_node(), ereignis_node);
          ereignis_node->append_attribute(doc.allocate_attribute("Er", "34"));
          ereignis_node->append_attribute(doc.allocate_attribute("Beschr", neues_register.verkn_refpunkt_modul.c_str()));
          ereignis_node->append_attribute(doc.allocate_attribute("Wert", doc.allocate_string(std::to_string(neues_register.verkn_refpunkt_nr).c_str())));
        }
      }
    }

    // Neue Referenzelemente
    auto* str_element_node = strecke_node->first_node("StrElement");
    for (const auto& neuer_refpunkt : aenderung.neue_refpunkte) {
      auto* refpunkt_node = doc.allocate_node(rapidxml::node_element, "ReferenzElemente");
      if (str_element_node != nullptr) {
        strecke_node->insert_node(str_element_node, refpunkt_node);
      } else {
        strecke_node->append_node(refpunkt_node);
      }
      refpunkt_node->append_attribute(doc.allocate_attribute("ReferenzNr", doc.allocate_string(std::to_string(neuer_refpunkt.nr).c_str())));
      refpunkt_node->append_attribute(doc.allocate_attribute("StrElement", doc.allocate_string(std::to_string(neuer_refpunkt.element_nr).c_str())));
      if (neuer_refpunkt.element_normrichtung) {
        refpunkt_node->append_attribute(doc.allocate_attribute("StrNorm", "1"));
      }
      refpunkt_node->append_attribute(doc.allocate_attribute("RefTyp", "2"));
      refpunkt_node->append_attribute(doc.allocate_attribute("Info", doc.allocate_string((std::string("Register-Nr. ") + std::to_string(neuer_refpunkt.register_nr)).c_str())));
    }

    std::string out_string;
    rapidxml::print(std::back_inserter(out_string), doc, rapidxml::print_no_indenting);
    {
      boost::nowide::ofstream o(dateiname, std::ios::binary);
      o << out_string;
    }

    boost::nowide::cout << " - " << dateiname << std::endl;
  }

  boost::nowide::cout << "In den folgenden Modulen müssen die Fahrstrassen neu erstellt werden:" << std::endl;
  for (const auto& strecke : geaenderte_module) {
    boost::nowide::cout << " - " << modul_info.at(strecke).pfad_os << std::endl;
  }
}
