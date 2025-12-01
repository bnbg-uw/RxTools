#pragma once

#ifndef rxtools_lmu_h
#define rxtools_lmu_h

#include "rxtools_pch.hpp"
#include "treatment.hpp"
#include "utilities.hpp"
#include "rxunit.hpp"
#include "ProcessedFolder/src/ProcessedFolder.hpp"

namespace rxtools {

enum class LmuType { valleyBottom, ridgeTop, neFacing, swFacing, all };

class Lmu {
public:
    std::vector<RxUnit> units;
    lapis::Raster<lapis::cell_t> mask;
    LmuType type = LmuType::all;
    
    //assigned reference area observed structures, and associated weights
    std::vector<std::string> targetLmuNames;
    std::vector<StructureSummary> structures;

    Lmu() {};
    Lmu(lapis::Raster<lapis::cell_t> mask, LmuType t) : mask(mask), type(t) {};
    Lmu(std::string path, TaoGetters<lapis::VectorDataset<lapis::Point>> getters);

    void makeUnits(lapis::VectorDataset<lapis::MultiPolygon> unitsPoly, TaoList<lapis::VectorDataset<lapis::Point>> tl, lapis::Raster<int> osiNum, lapis::Raster<int> osiDen, double convFactor, bool overrideTargets);

    template<class T>
    void makeUnits(lapis::Raster<T> unitsRaster, TaoList<lapis::VectorDataset<lapis::Point>> taos, lapis::Raster<int> osiNum, lapis::Raster<int> osiDen, double convFactor) {
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

    void write(std::string path, allometry::FastFuels ffa);
};
} //namespace rxtools 

#endif //!rxtools_lmu_h