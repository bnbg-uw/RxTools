#pragma once

#ifndef licosim_lmu_h
#define licosim_lmu_h

#include "shapefile/shapefile.h"
#include "licosim/treatment.hpp"
#include "licosim/utilities.hpp"
#include "ProcessedFolder/ProcessedFolder.hpp"

namespace licosim {

enum class LmuType { valleyBottom, ridgeTop, neFacing, swFacing, all };

class Lmu {
public:
    std::vector<RxUnit> units;
    spatial::Raster<int> mask;
    LmuType type = LmuType::all;
    
    //assigned reference area observed structures, and associated weights
    std::vector<std::string> targetLmuNames;
    std::vector<StructureSummary> structures;

    Lmu() {};
    Lmu(spatial::Raster<int> mask, LmuType t) : mask(mask), type(t) {};
    Lmu(std::string path, dbhFunction dbhf);

    void makeUnits(spatial::SpVectorDataset<spatial::SpMultiPolygon> unitsPoly, lico::TaoList tl, spatial::Raster<int> osiNum, spatial::Raster<int> osiDen, double convFactor, dbhFunction dbhFunc, bool overrideTargets);

    template<class T>
    void makeUnits(spatial::Raster<T> unitsRaster, lico::TaoList tl, spatial::Raster<int> osiNum, spatial::Raster<int> osiDen, double convFactor, dbhFunction dbhFunc) {
        if (units.size()) throw std::runtime_error("Units have already been calculated");

        if (!unitsRaster.overlaps(mask)) throw spatial::OutsideExtentException();
        unitsRaster.crop(mask);
        unitsRaster.mask(mask);
        std::unordered_set<T> ids;
        for (spatial::cell_t i = 0; i < unitsRaster.ncell(); i++) {
            if (unitsRaster[i].has_value())
                ids.emplace(unitsRaster[i].value());
        }

        for (T id : ids) {
            auto unitMask = unitsRaster;
            xt::filtration(unitMask.values().has_value(), xt::not_equal(unitMask.values().value(), id)) = false;
            unitMask.trim();

            osiNum.crop(unitMask);
            osiDen.crop(unitMask);
            osiNum.mask(unitMask);
            osiDen.mask(unitMask);
            double num = 0;
            double den = 0;
            for (spatial::cell_t i = 0; i < osiNum.ncell(); ++i) {
                if (osiNum[i].has_value()) {
                    num += osiNum[i].value();
                    den += osiDen[i].value();
                }
            }
            double osi = num / den * 100;

            if (unitMask.anyHasValue()) {
                auto rx = RxUnit(unitMask, tl, osi, convFactor, dbhFunc);
                if(rx.areaHa > 0.5)
                    units.push_back(rx);
            }
        }
    }

    void assignUnitTargets(std::default_random_engine dre, double bb_dbh, bool overrideTargets);

    void write(std::string path, fastFuelsAllometry ffa);
};
} //namespace licosim 

#endif //!licosim_lmu_h