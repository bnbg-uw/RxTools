#pragma once

#include "raster.hpp"
#include "rxunit.hpp"

namespace rxtools {
    class Output {
    public:  
        lapis::Raster<lapis::cell_t> ids;
        lapis::AttributeTable atts;
        std::vector<lapis::Raster<double>> pre;
        std::vector<lapis::Raster<double>> post;
        std::vector<lapis::Raster<double>> target;
        lapis::Raster<lapis::cell_t> lmus;
        lapis::Raster<lapis::cell_t> lmuIds;
        std::string commandLine;
        
        Output() {}
        Output(lapis::Alignment outAlign) : ids(outAlign), lmus(outAlign), lmuIds(outAlign){
            auto colnames = std::vector<std::string>{
                "ID",
                "curBaHa", "curBaAc", "curTpHa", "curTpAc", "curMCS", "curCC",
                "refBaHa", "refBaAc", "refTpHa", "refTpAc", "refMCS", "refCC",
                "trtBaHa", "trtBaAc", "trtTpHa", "trtTpAc", "trtMCS", "trtCC",
                "dbhMin", "dbhMax", "lmuCode"
            };
            auto types = std::vector<OGRFieldType>{
                OFTInteger,
                OFTReal, OFTReal, OFTReal, OFTReal, OFTReal, OFTReal,
                OFTReal, OFTReal, OFTReal, OFTReal, OFTReal, OFTReal,
                OFTReal, OFTReal, OFTReal, OFTReal, OFTReal, OFTReal,
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
            ids.overlayInside(rx.unitMask);
            lapis::cell_t thisId = -1;
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

                if (pre.size() < names.size()) {
                    pre.push_back(lapis::Raster<double>{ (lapis::Alignment)ids });
                    post.push_back(lapis::Raster<double>{ (lapis::Alignment)ids });
                    target.push_back(lapis::Raster<double>{ (lapis::Alignment)ids });
                }
                pre[i].overlayInside(rPre);
                post[i].overlayInside(rPost);
                target[i].overlayInside(rTarg);
            }

            atts.addRow();
            atts.setIntegerField(atts.nFeature()-1, "ID", thisId);
            atts.setRealField(atts.nFeature() - 1, "curBaHa", rx.currentStructure.ba);
            atts.setRealField(atts.nFeature() - 1, "curBaAc", rx.currentStructure.ba * 4.356);
            atts.setRealField(atts.nFeature() - 1, "curTpHa", rx.currentStructure.tph);
            atts.setRealField(atts.nFeature() - 1, "curTpAc", rx.currentStructure.tph / 2.47105);
            atts.setRealField(atts.nFeature() - 1, "curMCS", rx.currentStructure.mcs);
            atts.setRealField(atts.nFeature() - 1, "curCC", rx.currentStructure.cc);

            atts.setRealField(atts.nFeature() - 1, "refBaHa", rx.targetStructure.ba);
            atts.setRealField(atts.nFeature() - 1, "refBaAc", rx.targetStructure.ba * 4.356);
            atts.setRealField(atts.nFeature() - 1, "refTpHa", rx.targetStructure.tph);
            atts.setRealField(atts.nFeature() - 1, "refTpAc", rx.targetStructure.tph / 2.47105);
            atts.setRealField(atts.nFeature() - 1, "refMCS", rx.targetStructure.mcs);
            atts.setRealField(atts.nFeature() - 1, "refCC", rx.targetStructure.cc);

            atts.setRealField(atts.nFeature() - 1, "trtBaHa", rx.treatedStructure.ba);
            atts.setRealField(atts.nFeature() - 1, "trtBaAc", rx.treatedStructure.ba * 4.356);
            atts.setRealField(atts.nFeature() - 1, "trtTpHa", rx.treatedStructure.tph);
            atts.setRealField(atts.nFeature() - 1, "trtTpAc", rx.treatedStructure.tph / 2.47105);
            atts.setRealField(atts.nFeature() - 1, "trtMCS", rx.treatedStructure.mcs);
            atts.setRealField(atts.nFeature() - 1, "trtCC", rx.treatedStructure.cc);

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

            lmus.writeRaster((p / "lmus.tif").string());
            lmuIds.writeRaster((p / "lmuIds.tif").string());

            for (int i = 0; i < names.size(); ++i) {
                if (i == 0 || i == 1) {
                    pre[i].writeRaster((p / ("Rx_CurrentStructure_" + names[i] + "Ha.tif")).string());
                    post[i].writeRaster((p / ("Rx_TreatedStructure_" + names[i] + "Ha.tif")).string());
                    target[i].writeRaster((p / ("Rx_TargetStructure_" + names[i] + "Ha.tif")).string());

                    auto delta = post[i] - pre[i];
                    delta.writeRaster((p / ("Rx_DeltaStructure_" + names[i] + "Ha.tif")).string());

                    if (i == 0) {
                        (pre[i] * 4.356).writeRaster((p / ("Rx_CurrentStructure_" + names[i] + "Ac.tif")).string());
                        (post[i] * 4.356).writeRaster((p / ("Rx_TreatedStructure_" + names[i] + "Ac.tif")).string());
                        (target[i] * 4.356).writeRaster((p / ("Rx_TargetStructure_" + names[i] + "Ac.tif")).string());
                        (delta * 4.356).writeRaster((p / ("Rx_DeltaStructure_" + names[i] + "Ac.tif")).string());
                    }
                    else {
                        (pre[i] / 2.47105).writeRaster((p / ("Rx_CurrentStructure_" + names[i] + "Ac.tif")).string());
                        (post[i] / 2.47105).writeRaster((p / ("Rx_TreatedStructure_" + names[i] + "Ac.tif")).string());
                        (target[i] / 2.47105).writeRaster((p / ("Rx_TargetStructure_" + names[i] + "Ac.tif")).string());
                        (delta / 2.47105).writeRaster((p / ("Rx_DeltaStructure_" + names[i] + "Ac.tif")).string());
                    }
                }
                else {
                    pre[i].writeRaster((p / ("Rx_CurrentStructure_" + names[i] + ".tif")).string());
                    post[i].writeRaster((p / ("Rx_TreatedStructure_" + names[i] + ".tif")).string());
                    target[i].writeRaster((p / ("Rx_TargetStructure_" + names[i] + ".tif")).string());
                    auto delta = post[i] - pre[i];
                    delta.writeRaster((p / ("Rx_DeltaStructure_" + names[i] + ".tif")).string());

                }

                lapis::Raster<double> lmuStats{ (lapis::Alignment)lmuIds };
                lapis::Raster<double> lmuDelta{ (lapis::Alignment)lmuIds };

                auto preZonal = lapis::zonalMean(pre[i], lmuIds);
                for (auto z : preZonal) {
                    for (int j = 0; j < lmuStats.ncell(); ++j) {
                        if (lmuIds[j].value() == z.first) {
                            lmuStats[j].value() = z.second;
                            lmuStats[j].has_value() = true;
                        }
                    }
                }

                if (i == 0 || i == 1) {
                    lmuStats.writeRaster((p / ("LMU_CurrentStructure_" + names[i] + "Ha.tif")).string());

                    if (i == 0)
                        (lmuStats * 4.356).writeRaster((p / ("LMU_CurrentStructure_" + names[i] + "Ac.tif")).string());
                    else
                        (lmuStats / 2.47105).writeRaster((p / ("LMU_CurrentStructure_" + names[i] + "Ac.tif")).string());
                }
                else
                    lmuStats.writeRaster((p / ("LMU_CurrentStructure_" + names[i] + ".tif")).string());
                lmuDelta = lmuStats;

                auto postZonal = lapis::zonalMean(post[i], lmuIds);
                for (auto z : postZonal) {
                    for (int j = 0; j < lmuStats.ncell(); ++j) {
                        if (lmuIds[j].value() == z.first) {
                            lmuStats[j].value() = z.second;
                        }
                    }
                }
                lmuDelta = lmuStats - lmuDelta;

                if (i == 0 || i == 1) {
                    lmuStats.writeRaster((p / ("LMU_TreatedStructure_" + names[i] + "Ha.tif")).string());
                    lmuDelta.writeRaster((p / ("LMU_DeltaStructure_" + names[i] + "Ha.tif")).string());

                    if (i == 0) {
                        (lmuStats * 4.356).writeRaster((p / ("LMU_TreatedStructure_" + names[i] + "Ac.tif")).string());
                        (lmuDelta * 4.356).writeRaster((p / ("LMU_DeltaStructure_" + names[i] + "Ac.tif")).string());
                    }
                    else {
                        (lmuStats / 2.47105).writeRaster((p / ("LMU_TreatedStructure_" + names[i] + "Ac.tif")).string());
                        (lmuDelta / 2.47105).writeRaster((p / ("LMU_DeltaStructure_" + names[i] + "Ac.tif")).string());
                    }
                }
                else {
                    lmuStats.writeRaster((p / ("LMU_TreatedStructure_" + names[i] + ".tif")).string());
                    lmuDelta.writeRaster((p / ("LMU_DeltaStructure_" + names[i] + ".tif")).string());
                }

                auto targZonal = lapis::zonalMean(target[i], lmuIds);
                for (auto z : targZonal) {
                    for (int j = 0; j < lmuStats.ncell(); ++j) {
                        if (lmuIds[j].value() == z.first) {
                            lmuStats[j].value() = z.second;
                        }
                    }
                }

                if (i == 0 || i == 1) {
                    lmuStats.writeRaster((p / ("LMU_TargetStructure_" + names[i] + "Ha.tif")).string());

                    if (i == 0)
                        (lmuStats * 4.356).writeRaster((p / ("LMU_TargetStructure_" + names[i] + "Ac.tif")).string());
                    else
                        (lmuStats / 2.47105).writeRaster((p / ("LMU_TargetStructure_" + names[i] + "Ac.tif")).string());
                }
                else
                    lmuStats.writeRaster((p / ("LMU_TargetStructure_" + names[i] + ".tif")).string());
            }

            auto shp = lapis::rasterToMultiPolygonForTaos(ids, &atts);
            shp.writeShapefile((p / "licosim_units.shp").string());
        }

    private:
        std::vector<std::string> names = { "Ba", "Tp", "Mcs", "Cc"};
    };
} //namespace rxtools