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
// durch das Intervall [`min_element_seq_idx`, `max_element_seq_idx`] (jeweils einschliesslich) beschrieben werden.
struct ElementSeqFahrstr {
  const Fahrstrasse* fahrstrasse;
  size_t min_element_seq_idx;
  size_t max_element_seq_idx;
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
  for (const auto& r : neue_register) {
    if ((r.first >= fs1.min_element_seq_idx) && (r.first <= fs1.max_element_seq_idx) &&
        (r.second >= fs2.min_element_seq_idx) && (r.second <= fs2.max_element_seq_idx)) {
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

int main(int argc, char** argv) {
  boost::nowide::args a(argc, argv);

  std::string start_st3_pfad;
  int start_element_nr;
  std::string ziel_st3_pfad;
  int ziel_element_nr;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
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
    boost::nowide::cout << desc << "\n";
    return 1;
  }

  const auto& start_st3 = get_strecke(start_st3_pfad);
  if (!start_st3)
  {
    boost::nowide::cerr << "Fehler beim Laden von " << start_st3_pfad << "\n";
    return 1;
  }

  // Betrachte auch Fahrstrassen, die in Nachbarmodulen des Start- und Zielmoduls beginnen und enden
  for (const auto& modulDatei : start_st3->children_ModulDateien) {
    get_strecke(modulDatei->Datei);
  }

  const auto& ziel_st3 = get_strecke(ziel_st3_pfad);
  if (!ziel_st3)
  {
    boost::nowide::cerr << "Fehler beim Laden von " << ziel_st3_pfad << "\n";
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
    boost::nowide::cerr << "Element " << start_element_nr << " existiert nicht in Datei " << start_st3_pfad << "\n";
    return 1;
  }

  const StrElement* ziel_element = get_element_by_nr(*ziel_st3, ziel_element_nr);
  if (!ziel_element)
  {
    boost::nowide::cerr << "Element " << ziel_element_nr << " existiert nicht in Datei " << ziel_st3_pfad << "\n";
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
        boost::nowide::cout << " - Fahrwegsuche in " << (start_normrichtung ? "Norm" : "Gegen") << "richtung erfolgreich\n";
        return true;
      }

      if (idx > 0 && (cur_element.gegenrichtung().nachfolger_selbes_modul().size() + cur_element.gegenrichtung().nachfolger_anderes_modul().size() > 1)) {
        // Weiche in der Gegenrichtung
        boost::nowide::cout << " - Fahrwegsuche in " << (start_normrichtung ? "Norm" : "Gegen") << "richtung: ignoriere stumpf befahrene Weiche an " << modul_info.at(cur_element.st3).pfad_kurz << ", Element " << cur_element.element->Nr << '\n';
      }

      if (cur_element.nachfolger_selbes_modul().size() + cur_element.nachfolger_anderes_modul().size() == 0) {
        // Streckenende
        boost::nowide::cout << " - Fahrwegsuche in " << (start_normrichtung ? "Norm" : "Gegen") << "richtung abgebrochen an " << modul_info.at(cur_element.st3).pfad_kurz << ", Element " << cur_element.element->Nr << " (kein Nachfolger)\n";
        return false;
      } else if (cur_element.nachfolger_selbes_modul().size() + cur_element.nachfolger_anderes_modul().size() > 1) {
        boost::nowide::cout << " - Fahrwegsuche in " << (start_normrichtung ? "Norm" : "Gegen") << "richtung: ignoriere spitz befahrene Weiche an " << modul_info.at(cur_element.st3).pfad_kurz << ", Element " << cur_element.element->Nr << " (nimm ersten Nachfolger)\n";
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
        break;
      }
      cur_element = next_element;
    }
    assert(cur_element);
    boost::nowide::cout << " - Fahrwegsuche in " << (start_normrichtung ? "Norm" : "Gegen") << "richtung abgebrochen an " << modul_info.at(cur_element.st3).pfad_kurz << ", Element " << cur_element.element->Nr << " (Zielelement nach nach " << max_num_elemente << " Schritten nicht erreicht)\n";

    return false;
  };

  if (!findeWegZuZiel(true) && !findeWegZuZiel(false)) {
    boost::nowide::cerr << "Keinen Weg zwischen den angegebenen Streckenelementen gefunden, der nicht über eine Weiche führt\n";
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

      bool in_element_seq = false;
      bool element_seq_in_normrichtung = false;
      size_t min_element_seq_idx = 0;
      size_t max_element_seq_idx = 0;

      while (cur_element && (idx <= 9999)) {
        const auto element_seq_idx = element_seq.index(cur_element.element);
        if ((idx != 0) && element_seq_idx.has_value()) {  // idx != 0: Startelement ist selbst nicht in der Fahrstrasse enthalten
          if (!in_element_seq) {
            in_element_seq = true;
            element_seq_in_normrichtung = (cur_element.normrichtung == element_seq.elemente().at(*element_seq_idx).normrichtung);
            min_element_seq_idx = *element_seq_idx;
            max_element_seq_idx = *element_seq_idx;
          } else {
            min_element_seq_idx = std::min(min_element_seq_idx, *element_seq_idx);
            max_element_seq_idx = std::max(max_element_seq_idx, *element_seq_idx);
          }
        }

        idx++;
        if (cur_element == ziel_element) {
          if (in_element_seq) {
            ElementSeqFahrstr seq_fahrstr { fahrstrasse.get(), min_element_seq_idx, max_element_seq_idx };
            if (element_seq_in_normrichtung) {
              fahrstr_normrichtung.push_back(std::move(seq_fahrstr));
            } else {
              fahrstr_gegenrichtung.push_back(std::move(seq_fahrstr));
            }
          }
          break;
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
      if ((!richtungsInfo) || (richtungsInfo->Reg != 0)) {
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

  std::sort(fahrstr_normrichtung.begin(), fahrstr_normrichtung.end(), [](const auto& fs1, const auto& fs2) {
      return fs1.min_element_seq_idx > fs2.min_element_seq_idx;
      });
  std::sort(fahrstr_gegenrichtung.begin(), fahrstr_gegenrichtung.end(), [](const auto& fs1, const auto& fs2) {
      return fs1.max_element_seq_idx < fs2.max_element_seq_idx;
      });

  for (const auto& seq_fahrstr_norm : fahrstr_normrichtung) {
    for (const auto& seq_fahrstr_gegen : fahrstr_gegenrichtung) {
      if ((seq_fahrstr_gegen.max_element_seq_idx >= seq_fahrstr_norm.min_element_seq_idx) && !gegeneinander_verriegelt(seq_fahrstr_norm, seq_fahrstr_gegen, neue_register)) {
        size_t idx_norm = seq_fahrstr_norm.min_element_seq_idx;
        while (!freie_elemente_laufrichtung[idx_norm] && (idx_norm <= seq_fahrstr_norm.max_element_seq_idx)) {
          ++idx_norm;
        }

        size_t idx_gegen_plus_1 = seq_fahrstr_gegen.max_element_seq_idx + 1;  // + 1, damit es kleiner werden kann als seq_fahrstr_gegen.min_element_seq_idx + 1, das hoechstens 1 ist.
        while (!freie_elemente_gegenrichtung[idx_gegen_plus_1 - 1] && (idx_gegen_plus_1 >= seq_fahrstr_gegen.min_element_seq_idx + 1)) {
          --idx_gegen_plus_1;
        }

        if ((idx_norm > seq_fahrstr_norm.max_element_seq_idx) || (idx_gegen_plus_1 < seq_fahrstr_gegen.min_element_seq_idx + 1)) {
          boost::nowide::cerr << "Kein freies Element gefunden, um die Fahrstraßen \"" << seq_fahrstr_norm.fahrstrasse->FahrstrName << "\" und \"" << seq_fahrstr_gegen.fahrstrasse->FahrstrName << "\" gegeneinander zu verriegeln\n";
          return 1;
        }

        neue_register.emplace_back(idx_norm, idx_gegen_plus_1 - 1);
        geaenderte_module.insert(get_strecke(seq_fahrstr_norm.fahrstrasse->FahrstrStart->Datei));
        geaenderte_module.insert(get_strecke(seq_fahrstr_gegen.fahrstrasse->FahrstrStart->Datei));
        freie_elemente_laufrichtung[idx_norm] = false;
        freie_elemente_gegenrichtung[idx_gegen_plus_1 - 1] = false;
      }
    }
  }

  // Umsetzung Verriegelung zweier Elemente:
  // Wenn in derselben Moduldatei: neues manuelles Register in beiden Elementen (jeweils eine Richtung)
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

    const auto register_el1 = neuesManuellesRegister(el1.st3);
    const auto refpunkt_el1 = neuerRefpunkt(el1.st3);
    get_pair_element(strecken_aenderungen[el1.st3].neue_register[el1.element->Nr],  el1.normrichtung) = NeuesRegister { register_el1, 0, "" };
    strecken_aenderungen[el1.st3].neue_refpunkte.push_back(NeuerRegisterRefpunkt { refpunkt_el1, el1.element->Nr, el1.normrichtung, register_el1 });

    if (el1.st3 == el2.st3) {
      boost::nowide::cout << " - Neues manuelles Register " << register_el1 << " an " << modul1.pfad_kurz << ", Element " << el1.element->Nr << " " << (el1.normrichtung ? "blau" : "grün") << " und " << el2.element->Nr << " " << (!el2.normrichtung ? "blau" : "grün") << "\n";
      get_pair_element(strecken_aenderungen[el2.st3].neue_register[el2.element->Nr], !el2.normrichtung) = NeuesRegister { register_el1, 0, "" };
    } else {
      boost::nowide::cout << " - Neues manuelles Register " << register_el1 << " an " << modul1.pfad_kurz << ", Element " << el1.element->Nr << " " << (el1.normrichtung ? "blau" : "grün") << "\n";
      boost::nowide::cout << " - Neuer Referenzpunkt " << refpunkt_el1 << " (Typ Register) an " << modul1.pfad_kurz << ", Element " << el1.element->Nr << " " << (el1.normrichtung ? "blau" : "grün") << "\n";
      boost::nowide::cout << " - Neues Ereignis \"Register in Fahrstrasse verknuepfen A\" mit Referenznummer " << refpunkt_el1 << " und Moduldatei " << modul1.pfad_zusi << " an " << modul2.pfad_kurz << ", Element " << el2.element->Nr << " " << (!el2.normrichtung ? "blau" : "grün") << "\n";
      get_pair_element(strecken_aenderungen[el2.st3].neue_register[el2.element->Nr], !el2.normrichtung) = NeuesRegister { 0, refpunkt_el1, modul1.pfad_zusi };
    }
  }

  if (strecken_aenderungen.empty()) {
    return 0;
  }

  boost::nowide::cout << "Die folgenden Dateien wurden neu geschrieben:\n";
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
      std::ofstream o(dateiname, std::ios::binary);
      o << out_string;
    }

    boost::nowide::cout << " - " << dateiname << "\n";
  }

  boost::nowide::cout << "In den folgenden Modulen müssen die Fahrstrassen neu erstellt werden:\n";
  for (const auto& strecke : geaenderte_module) {
    boost::nowide::cout << " - " << modul_info.at(strecke).pfad_os << "\n";
  }
}
