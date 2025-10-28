#include "lmu.hpp"


namespace licosim {

    Lmu::Lmu(std::string path, dbhFunction dbhf) {
        try {
            mask = spatial::Raster<int>(path + "ridgeTop.img");
            type = LmuType::ridgeTop;
        }
        catch (spatial::InvalidRasterFileException e) {
            try {
                mask = spatial::Raster<int>(path + "valleyBottom.img");
                type = LmuType::valleyBottom;
            }
            catch (spatial::InvalidRasterFileException e) {
                try {
                    mask = spatial::Raster<int>(path + "swFacing.img");
                    type = LmuType::swFacing;
                }
                catch (spatial::InvalidRasterFileException e) {
                    mask = spatial::Raster<int>(path + "nwFacing.img");
                    type = LmuType::neFacing;
                }
            }
        }
        auto p = boost::filesystem::path(path);
        p /= "units";
        for (auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(p), {})) {
            auto rx = RxUnit(entry.path().string(), dbhf);
            units.push_back(rx);
        }
    }

    void Lmu::makeUnits(spatial::SpVectorDataset<spatial::SpMultiPolygon> unitsPoly, lico::TaoList tl, spatial::Raster<int> osiNum, spatial::Raster<int> osiDen, double convFactor, dbhFunction dbhFunc, bool overrideTargets) {
        if (units.size()) throw std::runtime_error("Units have already been calculated");
        units.reserve(unitsPoly.nFeatures());
        int baCol = -1;
        int dbhMinCol = -1;
        int dbhMaxCol = -1;
        for (int i = 0; i < unitsPoly.colNames().size(); i++) {
            if (unitsPoly.colNames()[i] == "ba")
                baCol = i;
            if (unitsPoly.colNames()[i] == "dbhMin")
                dbhMinCol = i;
            if (unitsPoly.colNames()[i] == "dbhMax")
                dbhMaxCol = i;
        }
        
        for (int i = 0; i < unitsPoly.nFeatures(); i++) {
            if (unitsPoly[i].geom.overlaps(mask)) {
                auto unitMask = mask;
                unitMask.extend(unitsPoly[i].geom);
                auto x = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                unitMask.crop(unitsPoly[i].geom);
                unitMask.mask(unitsPoly[i].geom);
                if (!unitMask.anyHasValue())
                    continue;

                unitMask.trim();

                auto thisNum = spatial::crop(osiNum, unitMask);
                auto thisDen = spatial::crop(osiDen, unitMask);
                thisNum.mask(unitMask);
                thisDen.mask(unitMask);
                double num = 0;
                double den = 0;


                for (spatial::cell_t x = 0; x < thisNum.ncell(); ++x) {
                    if (thisNum[x].has_value()) {
                        num += thisNum[x].value();
                        den += thisDen[x].value();
                    }
                }
                double osi = num / den * 100;

                try {
                    auto rx = RxUnit(unitMask, tl, osi, convFactor, dbhFunc);
                    if (rx.areaHa > 0.5) {
                        units.push_back(rx);
                        if (overrideTargets) {
                            if (baCol >= 0)
                                units[units.size() - 1].targetStructure.ba = unitsPoly[i].getAttributeAsDouble(baCol);
                            if (dbhMaxCol >= 0)
                                units[units.size() - 1].dbhMax = unitsPoly[i].getAttributeAsDouble(dbhMaxCol);
                            if (dbhMinCol >= 0)
                                units[units.size() - 1].dbhMin = unitsPoly[i].getAttributeAsDouble(dbhMinCol);
                        }
                    }
                }
                catch (lidar::FileNotFoundException e) {
                    continue;
                }
            }
        }
    }

    void Lmu::assignUnitTargets(std::default_random_engine dre, double bbDbh, bool overrideTargets) {
        for (int i = 0; i < units.size(); i++) {
            auto minMax = units[i].getVirtualMinMax(dre, bbDbh);
            //std::cout << "Min " << minMax.first.ba << " " << minMax.first.tph << " " << minMax.first.mcs << " " << minMax.first.osi << "\n";
            //std::cout << "Max " << minMax.second.ba << " " << minMax.second.tph << " " << minMax.second.mcs << " " << minMax.second.osi << "\n";

            std::vector<int> outerIdx;
            double minDist = std::numeric_limits<double>::max();
            int minIdx = 0;
            for(int j = 0; j < structures.size(); ++j) {
                //std::cout << structures[j].ba << " " << structures[j].tph << " " << structures[j].mcs << " " << structures[j].osi << "\n";
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
                    for (int k = 0; k < 5; ++k) {
                        if (k == 2 || k == 3) continue; //Skip mcs and osi.
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
                std::vector<double> w(1, outerIdx.size());
                j = outerIdx[randomIdxFromWeights(w, dre)];
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

    void Lmu::write(std::string path, fastFuelsAllometry ffa) {
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
        mask.writeRaster(path + "/" + t + ".img");
        
        std::ofstream out;
        out.open(path + "/structures.csv");
        for (int i = 0; i < structures.size(); ++i) {
            out << structures[i].ba << "," << structures[i].tph << "," << structures[i].mcs << "," << structures[i].osi;
            for (int j = 0; j < structures[i].csd.size(); ++j) {
                out << "," << structures[i].binMins[j] << "," << structures[i].csd[j] << "," << structures[i].binMaxs[j];
            }
            out << "\n";
        }
        out.close();

        for (int i = 0; i < units.size(); ++i) {
            auto p = path + "/units/" + std::to_string(i);
            boost::filesystem::create_directories(p);
            units[i].write(p, ffa);
        }
    }
}