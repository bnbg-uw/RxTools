#pragma once

#include "raster/raster.hpp"
#include "landscapemetrics/landscapemetrics.hpp"
#include "licosim/rxunit.hpp"

namespace licosim {
    class Output {
    public:  
        std::vector<spatial::Raster<double>> pre;
        std::vector<spatial::Raster<double>> post;
        std::vector<spatial::Raster<double>> target;
        spatial::Raster<int> lmus;
        spatial::Raster<int> lmuIds;
        std::string commandLine;

        spatial::SpVectorDataset<spatial::SpMultiPolygon> shp;

        Output() {
            auto colnames = std::vector<std::string>{
                "id",
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
            shp.setColnames(colnames);
            shp.setDtypes(types);
        }

        void addRxUnit(RxUnit& rx, int lmuCode) {
            for (int i = 0; i < names.size(); ++i) {
                std::cout << i << "\n";
                spatial::Raster<double> rPre{ (spatial::Alignment)rx.unitMask };
                xt::filtration(rPre.values().value(), rx.unitMask.values().has_value()) = rx.currentStructure[i];
                rPre.values().has_value() = rx.unitMask.values().has_value();
                spatial::Raster<double> rPost{ (spatial::Alignment)rx.unitMask };
                xt::filtration(rPost.values().value(), rx.unitMask.values().has_value()) = rx.treatedStructure[i];
                rPost.values().has_value() = rx.unitMask.values().has_value();
                spatial::Raster<double> rTarg{ (spatial::Alignment)rx.unitMask };
                xt::filtration(rTarg.values().value(), rx.unitMask.values().has_value()) = rx.targetStructure[i];
                rTarg.values().has_value() = rx.unitMask.values().has_value();
                std::cout << "a\n";

                if (pre.size() == names.size()) {
                    std::vector<spatial::Raster<double>*> poPreToMerge;
                    std::vector<spatial::Raster<double>*> poPostToMerge;
                    std::vector<spatial::Raster<double>*> poTargToMerge;

                    //Add existing output layers
                    poPreToMerge.push_back(&pre.at(i));
                    poPostToMerge.push_back(&post.at(i));
                    poTargToMerge.push_back(&target.at(i));

                    //Add new layers
                    poPreToMerge.push_back(&rPre);
                    poPostToMerge.push_back(&rPost);
                    poTargToMerge.push_back(&rTarg);

                    pre[i] = spatial::rasterMerge(poPreToMerge, mergeFunc);
                    post[i] = spatial::rasterMerge(poPostToMerge, mergeFunc);
                    target[i] = spatial::rasterMerge(poTargToMerge, mergeFunc);
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
            std::cout << "b\n";

            spatial::SpFeature<spatial::SpMultiPolygon> ft;
            std::cout << "b1\n";
            spatial::SpMultiPolygon mp;
            std::cout << "b2\n";
            rx.unitMask.writeRaster("F:/unitmask.img");
            auto geom = spatial::polygonize(rx.unitMask, true, true);
            std::cout << "b3\n";
            for (int i = 0; i < geom.nrow(); ++i) {
                mp.addPolygon(geom.getFeaturesPtr()->at(i).geom);
            }
            std::cout << "b4\n";
            ft.geom = mp;
            std::cout << "b5\n";
            ft.geom.updateExtent();
            std::cout << "c\n";

            try {
                ft.addAttribute(shp.colNames()[0], shp.dTypes()[0], geom[0].getAttribute("value"));
            }
            catch (std::exception e) {
                std::cout << "Failed at addAttribute()\n";
                rx.unitMask.writeRaster("E:/r_tmp.img");
                geom.write("E:/shp_tmp.shp");
                throw std::runtime_error("polygonize");
            }
            ft.addAttribute(shp.colNames()[1], shp.dTypes()[1], std::to_string(rx.currentStructure.ba));
            ft.addAttribute(shp.colNames()[2], shp.dTypes()[2], std::to_string(rx.currentStructure.ba * 4.356));
            ft.addAttribute(shp.colNames()[3], shp.dTypes()[3], std::to_string(rx.currentStructure.tph));
            ft.addAttribute(shp.colNames()[4], shp.dTypes()[4], std::to_string(rx.currentStructure.tph / 2.47105));
            ft.addAttribute(shp.colNames()[5], shp.dTypes()[5], std::to_string(rx.currentStructure.mcs));
            ft.addAttribute(shp.colNames()[6], shp.dTypes()[6], std::to_string(rx.currentStructure.osi));
            ft.addAttribute(shp.colNames()[7], shp.dTypes()[7], std::to_string(rx.currentStructure.cc));

            ft.addAttribute(shp.colNames()[8], shp.dTypes()[8], std::to_string(rx.targetStructure.ba));
            ft.addAttribute(shp.colNames()[9], shp.dTypes()[9], std::to_string(rx.targetStructure.ba * 4.356));
            ft.addAttribute(shp.colNames()[10], shp.dTypes()[10], std::to_string(rx.targetStructure.tph));
            ft.addAttribute(shp.colNames()[11], shp.dTypes()[11], std::to_string(rx.targetStructure.tph / 2.47105));
            ft.addAttribute(shp.colNames()[12], shp.dTypes()[12], std::to_string(rx.targetStructure.mcs));
            ft.addAttribute(shp.colNames()[13], shp.dTypes()[13], std::to_string(rx.targetStructure.osi));
            ft.addAttribute(shp.colNames()[14], shp.dTypes()[14], std::to_string(rx.targetStructure.cc));

            ft.addAttribute(shp.colNames()[15], shp.dTypes()[15], std::to_string(rx.treatedStructure.ba));
            ft.addAttribute(shp.colNames()[16], shp.dTypes()[16], std::to_string(rx.treatedStructure.ba * 4.356));
            ft.addAttribute(shp.colNames()[17], shp.dTypes()[17], std::to_string(rx.treatedStructure.tph));
            ft.addAttribute(shp.colNames()[18], shp.dTypes()[18], std::to_string(rx.treatedStructure.tph / 2.47105));
            ft.addAttribute(shp.colNames()[19], shp.dTypes()[19], std::to_string(rx.treatedStructure.mcs));
            ft.addAttribute(shp.colNames()[20], shp.dTypes()[20], std::to_string(rx.treatedStructure.osi));
            ft.addAttribute(shp.colNames()[21], shp.dTypes()[21], std::to_string(rx.treatedStructure.cc));

            ft.addAttribute(shp.colNames()[22], shp.dTypes()[22], std::to_string(rx.dbhMin));
            ft.addAttribute(shp.colNames()[23], shp.dTypes()[23], std::to_string(rx.dbhMax));
            ft.addAttribute(shp.colNames()[24], shp.dTypes()[24], std::to_string(lmuCode));
            shp.addFeature(ft);
            std::cout << "d\n";

        }

        void write(std::string path) {
            auto p = boost::filesystem::path(path);
            std::ofstream cmdLine;
            cmdLine.open((p / "commandLine.txt").string());
            cmdLine << commandLine;
            cmdLine.close();

            lmus.writeRaster((p / "lmus.img").string());
            lmuIds.writeRaster((p / "lmuIds.img").string());
            spatial::Raster<double> lmuStats{ (spatial::Alignment)lmuIds };
            lmuStats.values().has_value() = lmuIds.values().has_value();
            spatial::Raster<double> lmuDelta{ (spatial::Alignment)lmuIds };
            lmuDelta.values().has_value() = lmuIds.values().has_value();

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

                auto preZonal = lsmetrics::zonalStatisticsByRaster(pre[i], lmuIds, zMean);
                for (auto z : preZonal)
                    xt::filtration(lmuStats.values().value(), xt::equal(lmuIds.values().value(), z.first)) = z.second.value();

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
                for (auto z : postZonal)
                    xt::filtration(lmuStats.values().value(), xt::equal(lmuIds.values().value(), z.first)) = z.second.value();

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
                for (auto z : targZonal)
                    xt::filtration(lmuStats.values().value(), xt::equal(lmuIds.values().value(), z.first)) = z.second.value();

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

            shp.projection(pre[0].projection());
            shp.write((p / "licosim_units.shp").string());
        }

    private:
        spatial::rasterMergeFunction<double, double> mergeFunc =spatial::rasterMergeFunction<double, double>(spatial::mergeByFirst<double>);
        lsmetrics::zonal_function<xtl::xoptional<double>, double> zMean = lsmetrics::zonalNoDataMean<double>;
        std::vector<std::string> names = { "Ba", "Tp", "Mcs", "Osi", "Cc"};
    };
} //namespace licosim