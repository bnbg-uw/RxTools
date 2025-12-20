#pragma once

#ifndef rxtools_rxunit_h
#define rxtools_rxunit_h

#include "taolist.hpp"
#include "structuresummary.hpp"
#include "models.hpp"

namespace rxtools {
    typedef boost::geometry::model::multi_point<StructureSummary> StructureMultiPoint_t;
    typedef boost::geometry::model::polygon<StructureSummary> StructurePolygon_t;

    // Columns and units for this csv should be:
    // SPCD     (species code as int)
    // Dia_cm   (diameter in centimeters)
    // HT_m     (height in meters)
    // STATUSCD (Live/Dead - 1 for live.  2 for dead I think? double check if mortality is introduced)
    // CBH_m    (Canopy base height in meters)
    // X_m, Y_m (X and Y locations in meters- idk if the units matter for this one?)
    inline void writeFastFuelsCsv(std::string path, TaoListMP taos, allometry::FastFuels ffa) {
        std::ofstream out;
        out.open(path);
        out << "SPCD,DIA_cm,HT_m,STATUSCD,CBH_m,X_m,Y_m\n";
        for (int i = 0; i < taos.size(); i++) {
            out << std::setprecision(std::numeric_limits<double>::max_digits10) << ffa.assignSpecies(taos.x(i), taos.y(i)) << "," << taos.dbh(i) << "," << taos.height(i) << ","
                << 1 << "," << ffa.predictCbh(taos.height(i)) << "," << taos.x(i) << "," << taos.y(i) << "\n";
        }
        out.close();
    }

class RxUnit {
public:
    bool paired = false;
    bool treated = false;
    TaoListMP taos;
    lapis::Raster<lapis::cell_t> unitMask;
    double areaHa = 0;

    double dbhMin = -1;
    double dbhMax = -1;

    StructureSummary currentStructure;
    StructureSummary targetStructure;

    TaoListMP treatedTaos;
    StructureSummary treatedStructure;

    RxUnit() = default;

    RxUnit(lapis::Raster<lapis::cell_t> mask, const TaoListMP& tl, double osi);

    RxUnit(std::string path, TaoGettersMP getters);

    std::pair<StructureSummary, StructureSummary> getVirtualMinMax(std::default_random_engine dre, double bbDbh);
    //StructurePolygon_t getTreatmentEnvelope(std::default_random_engine dre, std::string objective = "random", double bb_dbh = 21);
    //spatial::Raster<int> makeChm(spatial::Raster<int> chm, lico::TaoList tl);

    void write(std::string path, allometry::FastFuels ffa);

};
} // namespace rxtools

#endif // rxtools_rxunit_h