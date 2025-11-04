#pragma once

#ifndef licosim_rxunit_h
#define licosim_rxunit_h

#include "lico/LICO.hpp"
#include "landscapemetrics/landscapemetrics.hpp"
#include "shapefile/shapefile.h"
#include "licosim/structuresummary.hpp"
#include "licosim/allometry.hpp"
#include "lico/GraphLico.hpp"
#include "coregap/CoreGap.hpp"

namespace licosim {
    typedef boost::geometry::model::multi_point<StructureSummary> StructureMultiPoint_t;
    typedef boost::geometry::model::polygon<StructureSummary> StructurePolygon_t;

    using dbhFunction = std::function<std::vector<double>(lico::adapt_type<spatial::unit_t>)>;

    inline void writeFastFuelsCsv(std::string path, lico::TaoList tl, fastFuelsAllometry ffa, dbhFunction dbhFunc) {
        auto dbh = dbhFunc(tl.height());

        std::ofstream out;
        out.open(path);
        out << "SPCD,DIA_cm,HT_m,STATUSCD,CBH_m,X_m,Y_m\n";
        for (int i = 0; i < tl.size(); i++) {
            out << std::setprecision(std::numeric_limits<double>::max_digits10) << ffa.assignSpecies(tl.x()[i], tl.y()[i]) << "," << dbh[i] << "," << tl.height()[i] << ","
                << 1 << "," << ffa.predictCbh(tl.height()[i]) << "," << tl.x()[i] << "," << tl.y()[i] << "\n";
        }
        out.close();
    }

class RxUnit {
public:
    bool paired = false;
    bool treated = false;
    lico::TaoList taos{};
    spatial::Raster<int> unitMask;
    double areaHa = 0;
    dbhFunction dbhFunc;

    double canopycutoff = 2;
    double coregapdist = 6;
    double chmres = 0.75;

    double dbhMin = -1;
    double dbhMax = -1;

    StructureSummary currentStructure;
    StructureSummary targetStructure;

    lico::TaoList treatedTaos;
    StructureSummary treatedStructure;

    RxUnit() = default;

    template<class T>
    RxUnit(spatial::Raster<T> mask, lico::TaoList tl, double osi, double convFactor, dbhFunction dbhf) {
        dbhFunc = dbhf;

        unitMask = mask;
        for (spatial::cell_t i = 0; i < mask.ncell(); i++) {
            if (unitMask[i].has_value()) {
                areaHa += unitMask.xres() * unitMask.yres();
            }
        }

        areaHa *= convFactor * convFactor;
        areaHa /= 10000.0;

        for (int i = 0; i < tl.size(); i++) {
            if (unitMask.extract(tl.x()[i], tl.y()[i]).has_value())
                taos.addTAO(tl[i]);
        }

        coregapdist /= convFactor;
        currentStructure = summarizeStructure(taos, osi, dbhFunc);

    }

    RxUnit(std::string path, dbhFunction dbhf);

    std::pair<StructureSummary, StructureSummary> getVirtualMinMax(std::default_random_engine dre, double bbDbh);
    StructureSummary summarizeStructure(lico::TaoList& tl, double osi, dbhFunction dbhFunc) const;
    //StructurePolygon_t getTreatmentEnvelope(std::default_random_engine dre, std::string objective = "random", double bb_dbh = 21);
    //spatial::Raster<int> makeChm(spatial::Raster<int> chm, lico::TaoList tl);

    void write(std::string path, fastFuelsAllometry ffa);

private:
    double calcOsi(spatial::Raster<int> chm) const;

};
} // namespace licosim

#endif // licosim_rxunit_h