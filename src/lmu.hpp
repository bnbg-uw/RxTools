#pragma once

#ifndef rxtools_lmu_h
#define rxtools_lmu_h

#include "rxtools_pch.hpp"
#include "treatment.hpp"
#include "utilities.hpp"
#include "rxunit.hpp"
#include "ProcessedFolder.hpp"

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
    Lmu(std::string path, TaoGettersMP getters);

    void makeUnits(const lapis::VectorDataset<lapis::MultiPolygon>& unitsPoly, const TaoListMP& tl, const lapis::Raster<int>& osiNum, const lapis::Raster<int>& osiDen, const bool& overrideTargets);

    template<class T>
    void makeUnits(const lapis::Raster<T>& unitsRaster, const TaoListMP& tl, const lapis::Raster<int>& osiNum, const lapis::Raster<int>& osiDen) {
        if (units.size()) {
            throw std::runtime_error("Units have already been calculated");
        }
        if (!unitsRaster.overlaps(mask)) {
            throw lapis::OutsideExtentException();
        }

        auto thisUnits = lapis::cropRaster(unitsRaster, mask, lapis::SnapType::out);
        thisUnits.mask(mask);
        std::unordered_set<T> ids;
        for (lapis::cell_t i = 0; i < thisUnits.ncell(); i++) {
            if (thisUnits[i].has_value())
                ids.emplace(thisUnits[i].value());
        }

        for (T id : ids) {
            lapis::Raster<lapis::cell_t> unitMask{ (lapis::Alignment)thisUnits };
            for (lapis::cell_t c = 0; c < thisUnits.ncell(); ++c) {
                if (thisUnits[c].has_value() && thisUnits[c].value() == id) {
                    unitMask[c].has_value() = true;
                    unitMask[c].value() = static_cast<lapis::cell_t>(thisUnits[c].value());
                }
            }
            unitMask = lapis::trimRaster(unitMask);
            
            if (unitMask.hasAnyValue()) {
                auto thisNum = lapis::cropRaster(osiNum, unitMask, lapis::SnapType::out);
                auto thisDen = lapis::cropRaster(osiDen, unitMask, lapis::SnapType::out);
                thisNum.mask(unitMask);
                thisDen.mask(unitMask);
                double num = 0;
                double den = 0;
                for (lapis::cell_t i = 0; i < thisNum.ncell(); ++i) {
                    if (thisNum[i].has_value()) {
                        num += thisNum[i].value();
                        den += thisDen[i].value();
                    }
                }
                double osi = num / den * 100;
                
                auto rx = RxUnit(unitMask, tl, osi);
                if (rx.areaHa > 0.5) {
                    units.push_back(rx);
                }
            }
        }
    }

    void assignUnitTargets(std::default_random_engine dre, double bb_dbh, bool overrideTargets);

    void write(std::string path, allometry::FastFuels ffa);
};
} //namespace rxtools 

#endif //!rxtools_lmu_h