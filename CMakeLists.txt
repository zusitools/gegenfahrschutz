cmake_minimum_required (VERSION 3.1)
project (gegenfahrschutz)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(parser)
generate_zusi_parser(zusi_parser ${CMAKE_CURRENT_BINARY_DIR}/zusi_parser
  WHITELIST
    Zusi::Strecke

    Strecke::StrElement
    Strecke::ReferenzElemente
    Strecke::ModulDateien
    Strecke::Fahrstrasse

    StrElement::Nr
    StrElement::Anschluss
    StrElement::NachNorm
    StrElement::NachNormModul
    StrElement::NachGegen
    StrElement::NachGegenModul
    StrElement::InfoNormRichtung
    StrElement::InfoGegenRichtung

    StreckenelementRichtungsInfo::Reg
    StreckenelementRichtungsInfo::Ereignis

    Ereignis::Er

    NachfolgerAnderesModul::Datei
    NachfolgerAnderesModul::Nr
    NachfolgerSelbesModul::Nr

    Fahrstrasse::FahrstrName
    Fahrstrasse::FahrstrStart
    Fahrstrasse::FahrstrZiel
    Fahrstrasse::FahrstrWeiche
    Fahrstrasse::FahrstrRegister

    FahrstrStart::Datei
    FahrstrStart::Ref
    FahrstrZiel::Datei
    FahrstrZiel::Ref
    FahrstrWeiche::Datei
    FahrstrWeiche::Ref
    FahrstrWeiche::FahrstrWeichenlage
    FahrstrRegister::Datei
    FahrstrRegister::Ref

    Dateiverknuepfung::Dateiname

    ModulDateiVerknuepfung::Datei

    ReferenzElement::ReferenzNr
    ReferenzElement::StrElement
    ReferenzElement::StrNorm
  IGNORE_UNKNOWN)

find_package(Boost COMPONENTS filesystem program_options REQUIRED)

add_executable(gegenfahrschutz gegenfahrschutz.cpp nowide/src/iostream.cpp)
# TODO: target_link_libraries(Boost::filesystem) funktioniert derzeit unter MXE nicht
target_include_directories(gegenfahrschutz PRIVATE ${Boost_INCLUDE_DIRS} nowide/include)
target_link_libraries(gegenfahrschutz PRIVATE zusi_parser ${Boost_LIBRARIES})

install(TARGETS gegenfahrschutz RUNTIME DESTINATION bin)
