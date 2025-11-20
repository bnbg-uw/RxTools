#include "rxunit.hpp"

namespace rxtools {
    std::pair<StructureSummary, StructureSummary> RxUnit::getVirtualMinMax(std::default_random_engine dre, double bbDbh) {
        StructureSummary min;
        StructureSummary max;

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
        lapis::Alignment a{ unitMask.xmin(), unitMask.ymin(),(lapis::rowcol_t)std::ceil((unitMask.ymax() - unitMask.ymin()) / 3),
                    (lapis::rowcol_t)std::ceil((unitMask.xmax() - unitMask.xmin()) / 3), 3, 3 };
        lapis::lico::GraphLico g{ a };

        //add all bb trees.
        std::vector<size_t> idxs;
        for (size_t i = 0; i < taos.size(); ++i) {
            if (taos.dbh(i) > bbDbh) {
                g.addTAO(taos[i], dbh[i], lapis::lico::NodeStatus::on);
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
        for (size_t i = 0; i < g.nodes.size(); ++i) {
            if (g.nodes[i].status != lapis::lico::NodeStatus::on) {
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
            lapis::lico::GraphLico g2{ a };
            for (size_t i = 0; i < taos.size(); ++i) {
                if (taos.dbh(i) > bbDbh) {
                    g2.addTAO(taos[i], dbh[i], lapis::lico::NodeStatus::on);
                    continue;
                }
                g2.addTAO(taos[i], dbh[i]);
            }

            double thisNum = numCurrent;
            double thisDen = denCurrent;
            auto thisCsd = csdCurrent;

            std::shuffle(std::begin(idxs), std::end(idxs), dre);
            for (auto idx : idxs) {
                std::unordered_set<size_t> adjClumps;
                for (size_t adj : g2.nodes[idx].adjList) {
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

    void RxUnit::write(std::string path, allometry::FastFuels ffa) {
        taos.writeCsv(path + "/taos.csv");
        unitMask.writeRaster(path + "/unitMask.img");
        //chm.writeRaster(path + "/chm.img");
        //basinMap.writeRaster(path + "/basinmap.img");
        treatedTaos.writeCsv(path + "/treatedTaos.csv");
        //treatedChm.writeRaster(path + "/treatedChm.img");

        if (ffa.init) {
            writeFastFuelsCsv(path + "/taos_fastFuels.csv", taos, ffa);
            writeFastFuelsCsv(path + "/treatedTaos_fastFuels.csv", treatedTaos, ffa);
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

    RxUnit::RxUnit(std::string path, TaoGetters<lapis::VectorDataset<lapis::Point>> getters) {
        //std::cout << "dbhf load\n";
        taos = TaoList(path + "/taos.csv");
       // std::cout << "taos load\n";
        unitMask = lapis::Raster<int>(path + "/unitMask.img");
        //chm = spatial::Raster<double>(path + "/chm.img");
        //unitMask = spatial::Raster<int>(path + "/basinmap.img");
        treatedTaos = TaoList(path + "/treatedTaos.csv");
        //treatedChm = spatial::Raster<double>(path + "/treatedChm.img");

        //std::cout << "metaddata\n";
        std::filebuf fb;
        if (!fb.open(path + "/metadata.csv", std::ios::in)) throw std::runtime_error("Cannot open metdata file.");
        std::istream is{ &fb };

        treated = utilities::readCSVLine(is)[1] == "1";
        areaHa = std::stod(utilities::readCSVLine(is)[1]);
        canopycutoff = std::stod(utilities::readCSVLine(is)[1]);
        coregapdist = std::stod(utilities::readCSVLine(is)[1]);
        dbhMin = std::stod(utilities::readCSVLine(is)[1]);
        dbhMax = std::stod(utilities::readCSVLine(is)[1]);

        auto row = utilities::readCSVLine(is);
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
        
        row = utilities::readCSVLine(is);
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

        row = utilities::readCSVLine(is);
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