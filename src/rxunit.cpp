#include "licosim/rxunit.hpp"

namespace licosim {

    StructureSummary RxUnit::summarizeStructure(lico::TaoList& tl, double osi, dbhFunction dbhFunc) const {
        auto sd = lico::SparseDistMatrix();
        sd.nearDist(tl.x(), tl.y(), tl.crown(), false);

        lico::ClumpInfo clumps = sd.getClumps(tl.size());

        double ba;
        double mcs;
        double tph;
        double cc;
        try {
            mcs = 0;
            cc = 0;
            for (int i = 0; i < tl.size(); ++i) {
                mcs += clumps.clumpSize[clumps.clumpID[i]];
                cc += tl.area()[i];
            }
            mcs /= static_cast<double>(tl.size());

            cc = cc / 10000. / areaHa;

            auto dbhList = dbhFunc(tl.height()); // cm

            ba = 0;
            for (int i = 0; i < dbhList.size(); i++) {
                ba += M_PI * (dbhList[i] / 2.) * (dbhList[i] / 2.); // square cm
            }
            ba = ba / 10000. / areaHa; // to square meters then to sqm/ha
        }
        catch(std::bad_alloc e) {
            std::cout << "ba/mcs\n";
            throw e;
        }
        catch (std::exception e) {
            std::cout << e.what();
            throw e;
        }

        tph = tl.size() / areaHa;

        std::vector<double> csd(6, 0);

        try {

            for (int i = 0; i < clumps.clumpSize.size(); i++) { //TODO: don't hardcode these bins
                if (clumps.clumpSize[i] == 1)
                    csd[0]++;
                else if (clumps.clumpSize[i] < 5)
                    csd[1]++;
                else if (clumps.clumpSize[i] < 10)
                    csd[2]++;
                else if (clumps.clumpSize[i] < 15)
                    csd[3]++;
                else if (clumps.clumpSize[i] < 31)
                    csd[4]++;
                else
                    csd[5]++;
            }

            for (int i = 0; i < csd.size(); i++)
                csd[i] /= static_cast<double>(clumps.clumpSize.size());
        }
        catch (std::bad_alloc e) {
            std::cout << "csd\n";
            throw e;
        }

        auto out = StructureSummary(ba, tph, mcs, osi, cc, csd);
        return out;
    }

    double RxUnit::calcOsi(spatial::Raster<int> chm) const {
        spatial::Alignment thisalign{ unitMask };
        thisalign.crop(chm, spatial::SnapType::out);
        thisalign.extend(chm, spatial::SnapType::out);
        lsmetrics::crop_mw_function<double, int> numFunc = [&](const lsmetrics::crop_view<int>& e)->xtl::xoptional<double> {return lsmetrics::OSInumerator(e, chmres, canopycutoff, coregapdist); };
        spatial::Raster<double> coreNum = lsmetrics::movingWindowByRaster(chm, thisalign, numFunc, coregapdist);
        lsmetrics::crop_mw_function<double, int> totalFunc = [&](const lsmetrics::crop_view<int>& e)->xtl::xoptional<double> {return lsmetrics::totalAreaForOSI(e, chmres, coregapdist); };
        spatial::Raster<double> thistotal = lsmetrics::movingWindowByRaster(chm, thisalign, totalFunc, coregapdist);

        int osinum = 0;
        int osiden = 0;
        for (spatial::cell_t i = 0; i < coreNum.ncell(); ++i) {
            if (coreNum[i].has_value()) {
                osinum += coreNum[i].value();
                osiden += thistotal[i].value();
            }
        }
        return (static_cast<double>(osinum) / static_cast<double>(osiden)) * 100.;
    }

    std::pair<StructureSummary, StructureSummary> RxUnit::getVirtualMinMax(std::default_random_engine dre, double bbDbh) {
        StructureSummary min;
        StructureSummary max;
        auto dbh = dbhFunc(taos.height());

        max.ba = currentStructure.ba;
        max.tph = currentStructure.tph;
        max.mcs = currentStructure.mcs;
        max.cc = currentStructure.cc;
        
        //min.osi = calcOsi(thisChm);
        min.cc = 0;
        min.ba = 0;
        min.tph = 0;
        min.mcs = std::numeric_limits<double>::max();

        //calc min and max csd distributions
        spatial::Alignment a{ unitMask.xmin(), unitMask.ymin(),(spatial::rowcol_t)std::ceil((unitMask.ymax() - unitMask.ymin()) / 3),
                    (spatial::rowcol_t)std::ceil((unitMask.xmax() - unitMask.xmin()) / 3), 3, 3 };
        lico::GraphLico g{ a };

        //add all bb trees.
        std::vector<lico::index_t> idxs;
        for (lico::index_t i = 0; i < taos.size(); ++i) {
            if (dbh[i] > bbDbh) {
                g.addTAO(taos[i], dbh[i], lico::NodeStatus::on);
                min.ba += g.nodes[g.nodes.size() - 1].ba;
                min.tph++;
                min.cc += g.nodes[g.nodes.size() - 1].area;
                continue;
            }
            g.addTAO(taos[i], dbh[i]);
            idxs.push_back(i);
        }
        min.ba /= areaHa;
        min.tph /= areaHa;
        min.cc = min.cc / 10000. / areaHa;

        for (int i = 0; i < min.csd.size(); ++i) {
            min.csd[i] = std::numeric_limits<double>::max();
            max.csd[i] = std::numeric_limits<double>::min();
        }

        //get current structure and update mins and maxs.
        std::vector<double> csdCurrent(min.csd.size());
        double numCurrent = 0;
        int denCurrent = 0;
        for (lico::index_t i = 0; i < g.nodes.size(); ++i) {
            if (g.nodes[i].status != lico::NodeStatus::on) {
                continue;
            }
            int clumpSize = g.nodes[i].clumpSize();
            numCurrent += clumpSize;
            denCurrent++;
            for (int j = 0; j < csdCurrent.size(); ++j) {
                if (clumpSize >= min.binMins[j] && (clumpSize <= min.binMaxs[j] || j == csdCurrent.size() - 1)) {
                    csdCurrent[j] += g.nodes[i].ba;
                }
            }
        }
        for (int i = 0; i < csdCurrent.size(); ++i) {
            min.csd[i] = std::min(csdCurrent[i], min.csd[i]);
            max.csd[i] = std::max(csdCurrent[i], max.csd[i]);
        }
        double mcsCurrent = numCurrent / denCurrent;
        min.mcs = std::min(mcsCurrent, min.mcs);
        max.mcs = std::max(mcsCurrent, max.mcs);

        std::cout << numCurrent << " : " << denCurrent << "\n";
        for (int i = 0; i < 5; ++i) {
            lico::GraphLico g2{ a };
            for (lico::index_t i = 0; i < taos.size(); ++i) {
                if (dbh[i] > bbDbh) {
                    g2.addTAO(taos[i], dbh[i], lico::NodeStatus::on);
                    continue;
                }
                g2.addTAO(taos[i], dbh[i]);
            }

            double thisNum = numCurrent;
            double thisDen = denCurrent;
            auto thisCsd = csdCurrent;

            std::shuffle(std::begin(idxs), std::end(idxs), dre);
            for (auto idx : idxs) {
                std::unordered_set<lico::index_t> adjClumps;
                for (lico::index_t adj : g2.nodes[idx].adjList) {
                    if (adjClumps.count(g2.nodes[adj].findAncestor().index)) {
                        continue;
                    }
                    adjClumps.insert(g2.nodes[adj].findAncestor().index);
                    int size = g2.nodes[adj].clumpSize();
                    thisNum -= size * size;
                    double ba = g2.nodes[adj].getClumpBA();
                    for (int i = 0; i < thisCsd.size(); ++i) {
                        if (size >= min.binMins[i] && (size <= min.binMaxs[i] || i == thisCsd.size() - 1)) {
                            thisCsd[i] -= ba;
                        }
                    }
                }

                //turn on addTAO and update current ba and csd
                int size = g2.nodes[idx].turnOn();
                thisNum += size * size;
                thisDen++;
                double ba = g2.nodes[idx].getClumpBA();
                for (int i = 0; i < thisCsd.size(); ++i) {
                    if (size >= min.binMins[i] && (size <= min.binMaxs[i] || i == thisCsd.size() - 1)) {
                        thisCsd[i] += ba;
                    }
                }

                for (int i = 0; i < csdCurrent.size(); ++i) {
                    min.csd[i] = std::min(thisCsd[i], min.csd[i]);
                    max.csd[i] = std::max(thisCsd[i], max.csd[i]);
                }
                double thisMcs = thisNum / thisDen;
                min.mcs = std::min(thisMcs, min.mcs);
                max.mcs = std::max(thisMcs, max.mcs);
            }
        }
        return std::pair(min, max);
    }

    void RxUnit::write(std::string path, fastFuelsAllometry ffa) {
        taos.writeCsv(path + "/taos.csv");
        unitMask.writeRaster(path + "/unitMask.img");
        //chm.writeRaster(path + "/chm.img");
        //basinMap.writeRaster(path + "/basinmap.img");
        treatedTaos.writeCsv(path + "/treatedTaos.csv");
        //treatedChm.writeRaster(path + "/treatedChm.img");

        if (ffa.init) {
            writeFastFuelsCsv(path + "/taos_fastFuels.csv", taos, ffa, dbhFunc);
            writeFastFuelsCsv(path + "/treatedTaos_fastFuels.csv", treatedTaos, ffa, dbhFunc);
        }

        std::ofstream out;
        out.open(path + "/metadata.csv");
        out << "treated," << treated << "\n";
        out << "areaHa," << areaHa << "\n";
        out << "canopycutoff," << canopycutoff << "\n";
        out << "coregapdist," << coregapdist << "\n";
        out << "dbhMin," << dbhMin << "\n";
        out << "dbhMax," << dbhMax << "\n";

        out << "currentStructure," << currentStructure.ba << "," << currentStructure.tph << "," << currentStructure.mcs << "," << currentStructure.osi << "," << currentStructure.cc;
        for (int i = 0; i < currentStructure.csd.size(); ++i) {
            out << "," << currentStructure.binMins[i] << "," << currentStructure.csd[i] << "," << currentStructure.binMaxs[i];
        }
        out << "\n";

        out << "targetStructure," << targetStructure.ba << "," << targetStructure.tph << "," << targetStructure.mcs << "," << targetStructure.osi << "," << targetStructure.cc;
        for (int i = 0; i < targetStructure.csd.size(); ++i) {
            out << "," << targetStructure.binMins[i] << "," << targetStructure.csd[i] << "," << targetStructure.binMaxs[i];
        }
        out << "\n";

        out << "treatedStructure," << treatedStructure.ba << "," << treatedStructure.tph << "," << treatedStructure.mcs << "," << treatedStructure.osi << "," << treatedStructure.cc;
        for (int i = 0; i < treatedStructure.csd.size(); ++i) {
            out << "," << treatedStructure.binMins[i] << "," << treatedStructure.csd[i] << "," << treatedStructure.binMaxs[i];
        }
        out << "\n";

        out.close();
    }

    RxUnit::RxUnit(std::string path, dbhFunction dbhf) {
        dbhFunc = dbhf;
        //std::cout << "dbhf load\n";
        taos = lico::TaoList(path + "/taos.csv");
       // std::cout << "taos load\n";
        unitMask = spatial::Raster<int>(path + "/unitMask.img");
        //chm = spatial::Raster<double>(path + "/chm.img");
        //unitMask = spatial::Raster<int>(path + "/basinmap.img");
        treatedTaos = lico::TaoList(path + "/treatedTaos.csv");
        //treatedChm = spatial::Raster<double>(path + "/treatedChm.img");

        //std::cout << "metaddata\n";
        std::filebuf fb;
        if (!fb.open(path + "/metadata.csv", std::ios::in)) throw std::runtime_error("Cannot open metdata file.");
        std::istream is{ &fb };

        treated = readCSVLine(is)[1] == "1";
        areaHa = std::stod(readCSVLine(is)[1]);
        canopycutoff = std::stod(readCSVLine(is)[1]);
        coregapdist = std::stod(readCSVLine(is)[1]);
        dbhMin = std::stod(readCSVLine(is)[1]);
        dbhMax = std::stod(readCSVLine(is)[1]);

        auto row = readCSVLine(is);
        std::vector<int> binMins;
        std::vector<double> csd;
        std::vector<int> binMaxs;
        for (int i = 6; i < row.size(); ++i) {
            if (i % 3 == 1)
                csd.push_back(std::stod(row[i]));
            else if (i % 3 == 2)
                binMaxs.push_back(std::stoi(row[i]));
            else
                binMins.push_back(std::stoi(row[i]));
        }
        std::cout << csd.size() << " " << binMaxs.size() << " " << binMins.size() << "\n";
        currentStructure = StructureSummary(std::stod(row[1]), std::stod(row[2]), std::stod(row[3]), std::stod(row[4]), std::stod(row[5]), csd, binMins, binMaxs);
        
        row = readCSVLine(is);
        binMins.clear();
        csd.clear();
        binMaxs.clear();
        for (int i = 5; i < row.size(); ++i) {
            if (i % 3 == 1)
                csd.push_back(std::stod(row[i]));
            else if (i % 3 == 2)
                binMaxs.push_back(std::stoi(row[i]));
            else
                binMins.push_back(std::stoi(row[i]));
        }
        std::cout << csd.size() << " " << binMaxs.size() << " " << binMins.size() << "\n";

        targetStructure = StructureSummary(std::stod(row[1]), std::stod(row[2]), std::stod(row[3]), std::stod(row[4]), std::stod(row[5]), csd, binMins, binMaxs);

        row = readCSVLine(is);
        binMins.clear();
        csd.clear();
        binMaxs.clear();
        for (int i = 5; i < row.size(); ++i) {
            if (i % 3 == 1)
                csd.push_back(std::stod(row[i]));
            else if (i % 3 == 2)
                binMaxs.push_back(std::stoi(row[i]));
            else
                binMins.push_back(std::stoi(row[i]));
        }
        std::cout << csd.size() << " " << binMaxs.size() << " " << binMins.size() << "\n";

        treatedStructure = StructureSummary(std::stod(row[1]), std::stod(row[2]), std::stod(row[3]), std::stod(row[4]), std::stod(row[5]), csd, binMins, binMaxs);

        fb.close();
    }
}