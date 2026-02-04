#include "lmu.hpp"


namespace rxtools {

    Lmu::Lmu(std::string path, TaoGettersMP getters) {
        try {
            mask = lapis::Raster<lapis::cell_t>(path + "ridgeTop.tif");
            type = LmuType::ridgeTop;
        }
        catch (lapis::InvalidRasterFileException e) {
            try {
                mask = lapis::Raster<lapis::cell_t>(path + "valleyBottom.tif");
                type = LmuType::valleyBottom;
            }
            catch (lapis::InvalidRasterFileException e) {
                try {
                    mask = lapis::Raster<lapis::cell_t>(path + "swFacing.tif");
                    type = LmuType::swFacing;
                }
                catch (lapis::InvalidRasterFileException e) {
                    mask = lapis::Raster<lapis::cell_t>(path + "nwFacing.tif");
                    type = LmuType::neFacing;
                }
            }
        }
        auto p = std::filesystem::path(path);
        p /= "units";
        for (auto& entry : std::filesystem::directory_iterator(p)) {
            auto rx = RxUnit(entry.path().string(), getters);
            units.push_back(rx);
        }
    }

    void Lmu::makeUnits(const lapis::VectorDataset<lapis::MultiPolygon>& unitsPoly, const TaoListMP& tl, const bool& overrideTargets) {
        if (units.size()) {
            throw std::runtime_error("Units have already been calculated");
        }
        
        for (int i = 0; i < unitsPoly.nFeature(); i++) {
            if (unitsPoly.getFeature(i).getGeometry().boundingBox().overlaps(mask)) {
                auto unitMask = mask;
                unitMask = lapis::extendRaster(unitMask, unitsPoly.getFeature(i).getGeometry().boundingBox(), lapis::SnapType::out);
                unitMask = lapis::cropRaster(unitMask, unitsPoly.getFeature(i).getGeometry().boundingBox(), lapis::SnapType::out);
                unitMask.maskByMultiPolygon(unitsPoly.getFeature(i).getGeometry());
                if (!unitMask.hasAnyValue())
                    continue;
                unitMask = lapis::trimRaster(unitMask);

                try {
                    auto rx = RxUnit(unitMask, tl);
                    if (rx.areaHa > 0.5) {
                        units.push_back(rx);
                        if (overrideTargets) {
                            try {
                                units[units.size() - 1].targetStructure.ba = unitsPoly.getNumericField<double>(i, "ba");
                            }
                            catch (std::out_of_range) {}
                            try {
                                units[units.size() - 1].dbhMax = unitsPoly.getNumericField<double>(i, "dbhMax");
                            }
                            catch (std::out_of_range) {}
                            try {
                                units[units.size() - 1].dbhMin = unitsPoly.getNumericField<double>(i, "dbhMin");
                            }
                            catch (std::out_of_range) {}
                        }
                    }
                }
                catch (processedfolder::FileNotFoundException e) {
                    continue;
                }
            }
        }
    }

    void Lmu::assignUnitTargets(std::default_random_engine dre, double bbDbh, bool overrideTargets) {
        for (int i = 0; i < units.size(); i++) {
            auto minMax = units[i].getVirtualMinMax(dre, bbDbh);
            //std::cout << "Min " << minMax.first.ba << " " << minMax.first.tph << " " << minMax.first.mcs << "\n";
            //std::cout << "Max " << minMax.second.ba << " " << minMax.second.tph << " " << minMax.second.mcs << "\n";

            std::vector<int> outerIdx;
            double minDist = std::numeric_limits<double>::max();
            int minIdx = 0;
            for(int j = 0; j < structures.size(); ++j) {
                //std::cout << structures[j].ba << " " << structures[j].tph << " " << structures[j].mcs << << "\n";
                //Check for containment
                if (minMax.first.ba <= structures[j].ba && minMax.second.ba >= structures[j].ba &&
                    minMax.first.tph <= structures[j].tph && minMax.second.tph >= structures[j].tph &&
                    minMax.first.mcs <= structures[j].mcs && minMax.second.mcs >= structures[j].mcs &&
                    minMax.first.cc <= structures[j].cc && minMax.second.cc >= structures[j].cc)
                {
                    outerIdx.push_back(j); 
                }
                else if(outerIdx.size() == 0)
                { //Otherwise calc the distance if there are no contained points yet.
                    double lowerDist = 0;
                    double upperDist = 0;
                    for (int k = 0; k < 4; ++k) {
                        if (k == 2) continue; come back and fix this?
                        lowerDist += (structures[j][k] - minMax.first[k]) * (structures[j][k] - minMax.first[k]);
                        upperDist += (structures[j][k] - minMax.second[k]) * (structures[j][k] - minMax.second[k]);
                    }
                    for (int k = 0; k < minMax.first.csd.size(); ++k) {
                        lowerDist += (structures[j].csd[k] - minMax.first.csd[k]) * (structures[j].csd[k] - minMax.first.csd[k]);
                        upperDist += (structures[j].csd[k] - minMax.second.csd[k]) * (structures[j].csd[k] - minMax.second.csd[k]);
                    }
                    double thisDist = std::min(lowerDist, upperDist);
                    if (thisDist < minDist) {
                        minDist = thisDist;
                        minIdx = j;
                    }
                }
            }

            int j;
            if (outerIdx.size()) {
                std::cout << "Choosing random target from options\n";
                std::vector<double> w(outerIdx.size(), 1);
                j = outerIdx[utilities::randomIdxFromWeights(w, dre)];
                units[i].paired = true;
            }
            else {
                std::cout << "Choosing closest target.\n";
                j = minIdx;
            }

            double existingBa;
            if (overrideTargets && units[i].targetStructure.ba != 0)
                existingBa = units[i].targetStructure.ba;

            units[i].targetStructure = structures[j];
            if(overrideTargets)
                units[i].targetStructure.ba = existingBa;
        }
    }

    void Lmu::write(std::string path, allometry::FastFuels ffa) {
        std::string t;
        if (type == LmuType::all)
            t = "all";
        else if (type == LmuType::ridgeTop)
            t = "ridgeTop";
        else if (type == LmuType::valleyBottom)
            t = "valleyBottom";
        else if (type == LmuType::swFacing)
            t = "swFacing";
        else
            t = "nwFacing";
        mask.writeRaster(path + "/" + t + ".tif");
        
        std::ofstream out;
        out.open(path + "/structures.csv");
        for (int i = 0; i < structures.size(); ++i) {
            out << structures[i].ba << "," << structures[i].tph << "," << structures[i].mcs;
            for (int j = 0; j < structures[i].csd.size(); ++j) {
                out << "," << structures[i].binMins[j] << "," << structures[i].csd[j] << "," << structures[i].binMaxs[j];
            }
            out << "\n";
        }
        out.close();

        for (int i = 0; i < units.size(); ++i) {
            auto p = path + "/units/" + std::to_string(i);
            std::filesystem::create_directories(p);
            units[i].write(p, ffa);
        }
    }
}