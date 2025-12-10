#pragma once

#include "raster.hpp"
#include "rxunit.hpp"

namespace rxtools {
    class Output {
    public:  
        lapis::Raster<int> ids;
        lapis::AttributeTable atts;
        std::vector<lapis::Raster<double>> pre;
        std::vector<lapis::Raster<double>> post;
        std::vector<lapis::Raster<double>> target;
        lapis::Raster<lapis::cell_t> lmus;
        lapis::Raster<lapis::cell_t> lmuIds;
        std::string commandLine;

        Output() {
            auto colnames = std::vector<std::string>{
                "ID",
                "curBaHa", "curBaAc", "curTpHa", "curTpAc", "curMCS", "curOSI", "curCC",
                "refBaHa", "refBaAc", "refTpHa", "refTpAc", "refMCS", "refOSI", "refCC",
                "treatedBaHa", "treatedBaAc", "treatedTpHa", "treatedTpAc", "treatedMCS", "treatedOSI", "treatedCC",
                "dbhMin", "dbhMax", "lmuCode"
            };
            auto types = std::vector<OGRFieldType>{
                OFTString,
                OFTReal, OFTReal, OFTReal, OFTReal, OFTReal, OFTReal, OFTReal,
                OFTReal, OFTReal, OFTReal, OFTReal, OFTReal, OFTReal, OFTReal,
                OFTReal, OFTReal, OFTReal, OFTReal, OFTReal, OFTReal, OFTReal,
                OFTReal, OFTReal, OFTInteger
            };
            for (size_t i = 0; i < colnames.size(); ++i) {
                if (types.at(i) == OFTString) {
                    atts.addStringField(colnames.at(i), 64);
                }
                else if (types.at(i) == OFTReal) {
                    atts.addRealField(colnames.at(i));
                }
                else if (types.at(i) == OFTInteger) {
                    atts.addIntegerField(colnames.at(i));
                }
                else {
                    throw std::invalid_argument("unsupported OFT Type");
                }
            }
        }

        void addRxUnit(RxUnit& rx, int lmuCode) {
            std::vector<lapis::Raster<int>*> poIdsToMerge;
            poIdsToMerge.push_back(&ids);
            poIdsToMerge.push_back(&rx.unitMask);
            ids = lapis::mosaic(poIdsToMerge, lapis::firstCombiner<int>);
            int thisId = -1;
            for (lapis::cell_t i = 0; i < rx.unitMask.ncell(); ++i) {
                if (rx.unitMask[i].has_value()) {
                    thisId = rx.unitMask[i].value();
                }
            }
            if (thisId == -1) {
                throw std::invalid_argument("unitmask either has no value or has a negative value");
            }

            for (int i = 0; i < names.size(); ++i) {
                std::cout << i << "\n";
                lapis::Raster<double> rPre{ (lapis::Alignment)rx.unitMask };
                lapis::Raster<double> rPost{ (lapis::Alignment)rx.unitMask };
                lapis::Raster<double> rTarg{ (lapis::Alignment)rx.unitMask };
                for (lapis::cell_t j = 0; j < rPre.ncell(); ++j) {
                    if (rx.unitMask[j].has_value()) {
                        rPre[j].value() = rx.currentStructure[i];
                        rPost[j].value() = rx.treatedStructure[i];
                        rTarg[j].value() = rx.targetStructure[i];

                        rPre[j].has_value() = true;
                        rPost[j].has_value() = true;
                        rTarg[j].has_value() = true;
                    }
                }
                std::cout << "a\n";

                if (pre.size() == names.size()) {
                    std::vector<lapis::Raster<double>*> poPreToMerge;
                    std::vector<lapis::Raster<double>*> poPostToMerge;
                    std::vector<lapis::Raster<double>*> poTargToMerge;

                    //Add existing output layers
                    poPreToMerge.push_back(&pre.at(i));
                    poPostToMerge.push_back(&post.at(i));
                    poTargToMerge.push_back(&target.at(i));

                    //Add new layers
                    poPreToMerge.push_back(&rPre);
                    poPostToMerge.push_back(&rPost);
                    poTargToMerge.push_back(&rTarg);

                    pre[i] = lapis::mosaic(poPreToMerge, lapis::firstCombiner<double>);
                    post[i] = lapis::mosaic(poPostToMerge, lapis::firstCombiner<double>);
                    target[i] = lapis::mosaic(poTargToMerge, lapis::firstCombiner<double>);
                }
                else if(pre.size() < names.size()) {
                    pre.push_back(rPre);
                    post.push_back(rPost);
                    target.push_back(rTarg);
                }
                else {
                    std::cout << "Failed at merge " + std::to_string(pre.size()) + " " + std::to_string(names.size()) + "\n";
                    throw std::range_error("big ouch");
                }
            }

            atts.addRow();
            atts.setStringField(atts.nFeature()-1, "ID", std::to_string(thisId));
            atts.setRealField(atts.nFeature() - 1, "curBaHa", rx.currentStructure.ba);
            atts.setRealField(atts.nFeature() - 1, "curBaAc", rx.currentStructure.ba * 4.356);
            atts.setRealField(atts.nFeature() - 1, "curTpHa", rx.currentStructure.tph);
            atts.setRealField(atts.nFeature() - 1, "curTpAc", rx.currentStructure.tph / 2.47105);
            atts.setRealField(atts.nFeature() - 1, "curMCS", rx.currentStructure.mcs);
            atts.setRealField(atts.nFeature() - 1, "curOSI", rx.currentStructure.osi);
            atts.setRealField(atts.nFeature() - 1, "curCC", rx.currentStructure.cc);

            atts.setRealField(atts.nFeature() - 1, "refBaHa", rx.targetStructure.ba);
            atts.setRealField(atts.nFeature() - 1, "refBaAc", rx.targetStructure.ba * 4.356);
            atts.setRealField(atts.nFeature() - 1, "refTpHa", rx.targetStructure.tph);
            atts.setRealField(atts.nFeature() - 1, "refTpAc", rx.targetStructure.tph / 2.47105);
            atts.setRealField(atts.nFeature() - 1, "refMCS", rx.targetStructure.mcs);
            atts.setRealField(atts.nFeature() - 1, "refOSI", rx.targetStructure.osi);
            atts.setRealField(atts.nFeature() - 1, "refCC", rx.targetStructure.cc);

            atts.setRealField(atts.nFeature() - 1, "treatedBaHa", rx.treatedStructure.ba);
            atts.setRealField(atts.nFeature() - 1, "treatedBaAc", rx.treatedStructure.ba * 4.356);
            atts.setRealField(atts.nFeature() - 1, "treatedTpHa", rx.treatedStructure.tph);
            atts.setRealField(atts.nFeature() - 1, "treatedTpAc", rx.treatedStructure.tph / 2.47105);
            atts.setRealField(atts.nFeature() - 1, "treatedMCS", rx.treatedStructure.mcs);
            atts.setRealField(atts.nFeature() - 1, "treatedOSI", rx.treatedStructure.osi);
            atts.setRealField(atts.nFeature() - 1, "treatedCC", rx.treatedStructure.cc);

            atts.setRealField(atts.nFeature() - 1, "dbhMin", rx.dbhMin);
            atts.setRealField(atts.nFeature() - 1, "dbhMax", rx.dbhMax);
            atts.setIntegerField(atts.nFeature() - 1, "lmuCode", lmuCode);
            std::cout << "d\n";

        }

        void write(std::string path) {
            auto p = std::filesystem::path(path);
            std::ofstream cmdLine;
            cmdLine.open((p / "commandLine.txt").string());
            cmdLine << commandLine;
            cmdLine.close();

            lmus.writeRaster((p / "lmus.img").string());
            lmuIds.writeRaster((p / "lmuIds.img").string());
            lapis::Raster<double> lmuStats{ (lapis::Alignment)lmuIds };
            lapis::Raster<double> lmuDelta{ (lapis::Alignment)lmuIds };
            for (lapis::cell_t i = 0; i < lmuStats.ncell(); ++i) {
                if (lmuIds[i].has_value()) {
                    lmuStats[i].has_value() = true;
                    lmuDelta[i].has_value() = true;
                }
            }

            for (int i = 0; i < names.size(); ++i) {
                if (i == 0 || i == 1) {
                    pre[i].writeRaster((p / ("Rx_CurrentStructure_" + names[i] + "Ha.img")).string());
                    post[i].writeRaster((p / ("Rx_TreatedStructure_" + names[i] + "Ha.img")).string());
                    target[i].writeRaster((p / ("Rx_TargetStructure_" + names[i] + "Ha.img")).string());

                    auto delta = post[i] - pre[i];
                    delta.writeRaster((p / ("Rx_DeltaStructure_" + names[i] + "Ha.img")).string());

                    if (i == 0) {
                        (pre[i] * 4.356).writeRaster((p / ("Rx_CurrentStructure_" + names[i] + "Ac.img")).string());
                        (post[i] * 4.356).writeRaster((p / ("Rx_TreatedStructure_" + names[i] + "Ac.img")).string());
                        (target[i] * 4.356).writeRaster((p / ("Rx_TargetStructure_" + names[i] + "Ac.img")).string());
                        (delta * 4.356).writeRaster((p / ("Rx_DeltaStructure_" + names[i] + "Ac.img")).string());
                    }
                    else {
                        (pre[i] / 2.47105).writeRaster((p / ("Rx_CurrentStructure_" + names[i] + "Ac.img")).string());
                        (post[i] / 2.47105).writeRaster((p / ("Rx_TreatedStructure_" + names[i] + "Ac.img")).string());
                        (target[i] / 2.47105).writeRaster((p / ("Rx_TargetStructure_" + names[i] + "Ac.img")).string());
                        (delta / 2.47105).writeRaster((p / ("Rx_DeltaStructure_" + names[i] + "Ac.img")).string());
                    }
                }
                else {
                    pre[i].writeRaster((p / ("Rx_CurrentStructure_" + names[i] + ".img")).string());
                    post[i].writeRaster((p / ("Rx_TreatedStructure_" + names[i] + ".img")).string());
                    target[i].writeRaster((p / ("Rx_TargetStructure_" + names[i] + ".img")).string());
                    auto delta = post[i] - pre[i];
                    delta.writeRaster((p / ("Rx_DeltaStructure_" + names[i] + ".img")).string());

                }

                auto preZonal = lapis::zonalStatisticsByRaster(pre[i], lmuIds, zMean);
                for (auto z : preZonal) {
                    for (int j = 0; j < lmuStats.ncell(); ++j) {
                        if (lmuIds[j].value() == z.first) {
                            lmuStats[j].value() = z.second.value();
                        }
                    }
                }

                if (i == 0 || i == 1) {
                    lmuStats.writeRaster((p / ("LMU_CurrentStructure_" + names[i] + "Ha.img")).string());

                    if (i == 0)
                        (lmuStats * 4.356).writeRaster((p / ("LMU_CurrentStructure_" + names[i] + "Ac.img")).string());
                    else
                        (lmuStats / 2.47105).writeRaster((p / ("LMU_CurrentStructure_" + names[i] + "Ac.img")).string());
                }
                else
                    lmuStats.writeRaster((p / ("LMU_CurrentStructure_" + names[i] + ".img")).string());
                lmuDelta = lmuStats;

                auto postZonal = lsmetrics::zonalStatisticsByRaster(post[i], lmuIds, zMean);
                for (auto z : postZonal) {
                    for (int j = 0; j < lmuStats.ncell(); ++j) {
                        if (lmuIds[j].value() == z.first) {
                            lmuStats[j].value() = z.second.value();
                        }
                    }
                }
                lmuDelta = lmuStats - lmuDelta;

                if (i == 0 || i == 1) {
                    lmuStats.writeRaster((p / ("LMU_TreatedStructure_" + names[i] + "Ha.img")).string());
                    lmuDelta.writeRaster((p / ("LMU_DeltaStructure_" + names[i] + "Ha.img")).string());

                    if (i == 0) {
                        (lmuStats * 4.356).writeRaster((p / ("LMU_TreatedStructure_" + names[i] + "Ac.img")).string());
                        (lmuDelta * 4.356).writeRaster((p / ("LMU_DeltaStructure_" + names[i] + "Ac.img")).string());
                    }
                    else {
                        (lmuStats / 2.47105).writeRaster((p / ("LMU_TreatedStructure_" + names[i] + "Ac.img")).string());
                        (lmuDelta / 2.47105).writeRaster((p / ("LMU_DeltaStructure_" + names[i] + "Ac.img")).string());
                    }
                }
                else {
                    lmuStats.writeRaster((p / ("LMU_TreatedStructure_" + names[i] + ".img")).string());
                    lmuDelta.writeRaster((p / ("LMU_DeltaStructure_" + names[i] + ".img")).string());
                }

                auto targZonal = lsmetrics::zonalStatisticsByRaster(target[i], lmuIds, zMean);
                for (auto z : targZonal) {
                    for (int j = 0; j < lmuStats.ncell(); ++j) {
                        if (lmuIds[j].value() == z.first) {
                            lmuStats[j].value() = z.second.value();
                        }
                    }
                }

                if (i == 0 || i == 1) {
                    lmuStats.writeRaster((p / ("LMU_TargetStructure_" + names[i] + "Ha.img")).string());

                    if (i == 0)
                        (lmuStats * 4.356).writeRaster((p / ("LMU_TargetStructure_" + names[i] + "Ac.img")).string());
                    else
                        (lmuStats / 2.47105).writeRaster((p / ("LMU_TargetStructure_" + names[i] + "Ac.img")).string());
                }
                else
                    lmuStats.writeRaster((p / ("LMU_TargetStructure_" + names[i] + ".img")).string());
            }

            auto shp = lapis::rasterToMultiPolygonForTaos(ids, &atts);
            shp.writeShapefile((p / "licosim_units.shp").string());
        }

    private:
        lsmetrics::zonal_function<xtl::xoptional<double>, double> zMean = lsmetrics::zonalNoDataMean<double>;
        std::vector<std::string> names = { "Ba", "Tp", "Mcs", "Osi", "Cc"};
    };
} //namespace rxtools