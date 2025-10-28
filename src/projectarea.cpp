#include "licosim/projectarea.hpp"

namespace licosim {
    ProjectArea::ProjectArea(std::string lidarDatasetPath, std::string projectPolygonPath, std::string aetPath, std::string cwdPath, std::string tmnPath, int nThread, std::string lmuPath, std::string lmuParam) {
        lidarDataset = lidar::getProcessedFolderReader(lidarDatasetPath);
        
        projectPoly = spatial::SpVectorDataset<spatial::SpMultiPolygon>(projectPolygonPath);
        if (!spatial::consistentProjection(projectPoly.projection(), lidarDataset->getProjection()))
            projectPoly.project(lidarDataset->getProjection());


        auto mask = spatial::Raster<int>(lidar::stringOrThrow(lidarDataset->getMaskRaster()));
        mask.crop((spatial::Extent)projectPoly);

        std::cout << "\t\tReading and resampling climate layers...";
        auto r = spatial::Raster<double>(aetPath);
        aet = spatial::resampleBilinear(spatial::Raster<double>(aetPath), mask);
        cwd = spatial::resampleBilinear(spatial::Raster<double>(cwdPath), mask);
        tmn = spatial::resampleBilinear(spatial::Raster<double>(tmnPath), mask);
        std::cout << tmnPath << "\n";

        if (lmuPath.empty()) {
            std::cout << "\t\tCreating LMU's:\n";
            lmuRaster = createLmuRasterFromTpiAndAsp(lidarDataset, lmuParam, projectPoly);
            lmuIds = lsmetrics::connCompLabel(lmuRaster);
            std::cout << "\t\tLMU creation done!\n";
            //lmuRaster.writeRaster("E:/yubalmus.img");
        }
        else {
            // check extents etc...
            std::cout << "\t\tReading LMU raster...";
            lmuRaster = spatial::Raster<int>(lmuPath);
            if (!spatial::consistentProjection(lmuRaster.projection(), mask.projection()))
                lmuRaster = spatial::resampleNGB(lmuRaster, mask);
            else {
                lmuRaster.extend(mask);
                lmuRaster.crop(mask);
                lmuRaster.mask(mask);
            }
            lmuRaster.trim();
            lmuIds = lsmetrics::connCompLabel(lmuRaster);
            std::cout << " Done!\n";
        }

        osiNum = spatial::Raster<int>{ (spatial::Alignment)lmuIds };
        osiDen = spatial::Raster<int>{ (spatial::Alignment)lmuIds };
        bbOsiNum = spatial::Raster<int>{ (spatial::Alignment)lmuIds };
        bbOsiDen = spatial::Raster<int>{ (spatial::Alignment)lmuIds };

        std::cout << "\t\tCreating regiontype map..\n";
        for (spatial::cell_t c = 0; c < lmuIds.ncell(); ++c) {
            if (lmuIds[c].has_value()) {
                regionType.emplace(lmuIds[c].value(), lmuRaster[c].value());
            }
        }

        std::cout << "\t\tProject Area Construction done.\n";
    }

    void ProjectArea::createCoreGapAndReadTaos(int nThread, double bbDbh, dbhFunction dbhFunc) {
        std::cout << "\t\tcalculating pretreat OSI and reading taos...\n";
        std::pair<spatial::coord_t, spatial::coord_t> expectedRes{};
        if (lidarDataset->getConvFactor() == 1) {
            expectedRes.first = 0.75;
        }
        else {
            expectedRes.first = 2.4606;
        }
        expectedRes.second = 0;

        std::mutex mut{};
        std::vector<std::thread> threads{};
        int sofar = -1;
        auto mwThreadFunc = [&](int i) {
            try {
                coreGapThread(osiNum, osiDen, bbOsiNum, bbOsiDen, nThread, i, mut, sofar, lmuRaster, 2, 6, expectedRes, dbhFunc, bbDbh);
            }
            catch (lidar::FileNotFoundException e) {
                std::cerr << e.what() << '\n';
                std::cerr << "Aborting\n";
                throw e;
            }
            catch (spatial::InvalidRasterFileException e) {
                std::cerr << "Error reading file:\n";
                std::cerr << e.what() << '\n';
                std::cerr << "Aborting\b";
                throw e;
            }
        };
        for (int i = 0; i < nThread; ++i) {
            threads.push_back(std::thread(mwThreadFunc, i));
        }
        for (int i = 0; i < nThread; ++i) {
            threads[i].join();
        }

        osiNum.extend(lmuRaster);
        osiDen.extend(lmuRaster);
        osiNum.crop(lmuRaster);
        osiDen.crop(lmuRaster);
        osiNum.mask(lmuRaster);
        osiDen.mask(lmuRaster);

        bbOsiNum.extend(lmuRaster);
        bbOsiDen.extend(lmuRaster);
        bbOsiNum.crop(lmuRaster);
        bbOsiDen.crop(lmuRaster);
        bbOsiNum.mask(lmuRaster);
        bbOsiDen.mask(lmuRaster);
    }
    
    Lmu ProjectArea::createLmuThread(int& sofar, const int thisThread) {
        int id = std::next(regionType.begin(), sofar)->first;
        int type = std::next(regionType.begin(), sofar)->second;
        
        std::cout << "\t Creating lmu " + std::to_string(sofar) + "/" + std::to_string(regionType.size()) + " " + std::to_string(id) + " on thread " + std::to_string(thisThread) + "\n";
        auto r = lmuIds;
        xt::filter(r.values(), xt::not_equal(r.values().value(), id)).has_value() = false;
        r.trim();
        return Lmu(r, static_cast<LmuType>(type));
    }

    spatial::Raster<int> ProjectArea::createLmuRasterFromTpiAndAsp(std::unique_ptr<lidar::ProcessedFolder>& lds,
                                                                   std::string terrain,
                                                                   spatial::SpVectorDataset<spatial::SpMultiPolygon> projectPoly) {
        double ridgeSep, canyonSep;
        double ridgeArea = 20000;
        double canyonArea = 40000;
        int expSize = 0;
        if (terrain == "steep") {
            ridgeSep = 25; //ridges
            canyonSep = 20; //canyons
        }
        else {
            ridgeSep = 15;
            canyonSep = 14;
        }

        int tpiDist = 500;
        int aspDist = 135;
        if (lds->getConvFactor() != 1.) {
            ridgeArea /= (lds->getConvFactor() * lds->getConvFactor());
            canyonArea /= (lds->getConvFactor() * lds->getConvFactor());
        }

        std::cout << " TPI ";
        // create ridge groups
        std::cout << lds->getTPI(tpiDist).value() << "\n";
        std::cout << lds->getTPI(tpiDist).has_value() << "\n";
        auto tpi = spatial::Raster<double>(lidar::stringOrThrow(lds->getTPI(tpiDist))) / 100. * lds->getConvFactor();
        if (!projectPoly.isEmpty()) {
            tpi.crop((spatial::Extent)projectPoly);
            //tpi.mask(projectPoly);
        }

        std::cout << lds->getAspect(aspDist).has_value() << "\n";
        std::cout << lds->getAspect(aspDist).value() << "\n";
        auto aspect = spatial::Raster<double>(lidar::stringOrThrow(lds->getAspect(aspDist)));
        if (!projectPoly.isEmpty()) {
            aspect.crop((spatial::Extent)projectPoly);
            //aspect.mask(projectPoly);
        }
        aspect.mask(tpi);
        tpi.mask(aspect);

        auto ridge = tpi > ridgeSep;
        xt::filter(ridge.values(), !ridge.values().value()).has_value() = false;
        auto ridgeGroup = lsmetrics::connCompLabel(ridge);
        // Get unique names of ridges, then remove ridges smaller than min patch size.
        double cellArea = tpi.xres() * tpi.yres();
        int nCellArea = ridgeArea / cellArea;
        std::unordered_map<int, int> regionArea;
        for (auto v : ridgeGroup.values()) {
            if (v.has_value()) {
                regionArea.emplace(v.value(), 0);
                ++regionArea[v.value()];
            }
        }
        for (spatial::cell_t c = 0; c < ridgeGroup.ncell(); ++c) {
            if (regionArea[ridgeGroup[c].value()] < nCellArea) {
                ridgeGroup[c].has_value() = false;
            }
        }
        std::cout << "\t\t\tRidges identified.\n";

        //create canyon groups
        auto canyon = tpi < -canyonSep;
        xt::filter(canyon.values(), !canyon.values().value()).has_value() = false;
        auto canyonGroup = lsmetrics::connCompLabel(canyon);
        //get unique names of ridge the remove ridge
        cellArea = canyon.xres() * canyon.yres();
        nCellArea = canyonArea / cellArea;
        regionArea.clear();
        for (auto v : canyonGroup.values()) {
            if (v.has_value()) {
                regionArea.emplace(v.value(), 0);
                ++regionArea[v.value()];
            }
        }
        for (spatial::cell_t c = 0; c < canyonGroup.ncell(); ++c) {
            if (regionArea[canyonGroup[c].value()] < nCellArea) {
                canyonGroup[c].has_value() = false;
            }
        }
        std::cout << "\t\t\tValley bottoms identified.\n";


        // Calculate slope position
        spatial::Raster<int> spos{ (spatial::Alignment)tpi };
        for (spatial::cell_t i = 0; i < spos.ncell(); ++i) {
            //canyon = 0
            //ridge = 1
            //slope = 2
            if (canyon[i].has_value()) {
                spos[i].value() = 0; // canyon
            }
            else if (ridge[i].has_value()) {
                spos[i].value() = 1; // ridge
            }
            else {
                spos[i].value() = 2; // slope
            }
            spos[i].has_value() = true;
        }
        std::cout << "\t\t\tSlope position identified.\n";

        // Code aspect, ready to be put into output lmu raster.
        spatial::Raster<int> aspectClassified{ (spatial::Alignment)aspect };
        for (spatial::cell_t i = 0; i < aspect.ncell(); ++i) {
            if (aspect[i].has_value()) {
                aspectClassified[i].has_value() = true;
                if ((aspect[i].value() >= 315 && aspect[i].value() <= 360) ||
                    (aspect[i].value() >= 0 && aspect[i].value() < 135)) {
                    aspectClassified[i].value() = 2; // Northeasterly.
                }
                else {
                    aspectClassified[i].value() = 3; // Southwesterly.
                }
            }
        }
        std::cout << "\t\t\tAspect identified.\n";


        //propagate LMUs
        //canyon = 0
        //ridge = 1
        //NE facing = 2
        //SW facing = 3
        spatial::Raster<int> lmu{ (spatial::Alignment)aspect };
        for (spatial::cell_t i = 0; i < lmu.ncell(); ++i) {
            if (aspect[i].has_value()) {
                lmu[i].has_value() = true;
                if (spos[i].value() == 0 || spos[i].value() == 1) {
                    lmu[i].value() = spos[i].value();
                }
                else {
                    lmu[i].value() = aspectClassified[i].value();
                }
            }
        }
        std::cout << "\t\t\tFirst pass LMU's identified. Beginning post-processing.\n";

        ///------------------
        // LMU post processing.
        //--------------------
        // select large lmu's as seeds.
        auto lmuGrp = lsmetrics::connCompLabel(lmu);
        cellArea = canyon.xres() * canyon.yres();
        nCellArea = canyonArea / cellArea; // canyon area cutoff ~10acres, which is our lmu cutoff.
        regionArea.clear();
        for (auto v : lmuGrp.values()) {
            if (v.has_value()) {
                regionArea.emplace(v.value(), 0);
                ++regionArea[v.value()];
            }
        }
        for (spatial::cell_t c = 0; c < lmuGrp.ncell(); ++c) {
            if (regionArea[lmuGrp[c].value()] < nCellArea) {
                lmuGrp[c].has_value() = false;
            }
        }
        std::cout << "\t\t\tLarge lmu's identified (connected components created). Nibbling small lmu's into larger lmu's\n";

        //max dist needed to look would be radius of the circle representing the largest removable LMU.
        //a = pi*r^2, I'm not dividing by pi to leave buffer.
        auto lmuNibble = lsmetrics::nibble(lmu, lmuGrp, std::ceil(std::sqrt(nCellArea)));
        std::cout << "\t\t\tNibbling done.\n";

        return lmuNibble;
    }

    void ProjectArea::subdivideLmus(std::string climateClassPath, int nThread) {
        auto cc = spatial::Raster<int>(climateClassPath);
        cc = spatial::resampleNGB(cc, lmuIds);
        lmuIds = lmuIds * 100;
        lmuIds = lmuIds + cc;
        std::cout << "climate done \n";

        std::unordered_map<int, int> regionArea;
        for (auto v : lmuIds.values()) {
            if (v.has_value()) {
                regionArea.emplace(v.value(), 0);
                ++regionArea[v.value()];
            }
        }

        auto newlmus = lmuIds;
        std::cout << "threading\n";
        std::mutex mut{};
        int sofar = 0;
        std::vector<std::thread> threads{};
        auto threadFunc = [&](int i) { divideLmusThread(sofar, mut, regionArea, lmuIds, newlmus, i); };
        for (int i = 0; i < nThread; ++i) {
            threads.push_back(std::thread(threadFunc, i));
        }
        for (int i = 0; i < nThread; ++i) {
            threads[i].join();
        }
        std::cout << "threading done\n";
        regionType.clear();
        for (spatial::cell_t c = 0; c < lmuIds.ncell(); ++c) {
            if (lmuIds[c].has_value()) {
                regionType.emplace(lmuIds[c].value(), lmuRaster[c].value());
            }
        }
        auto nCellArea = 300; // canyon area cutoff ~10acres, which is our lmu cutoff.
        regionArea.clear();
        auto mask = newlmus;
        for (auto v : newlmus.values()) {
            if (v.has_value()) {
                regionArea.emplace(v.value(), 0);
                ++regionArea[v.value()];
            }
        }
        for (spatial::cell_t c = 0; c < newlmus.ncell(); ++c) {
            if (regionArea[newlmus[c].value()] < nCellArea) {
                mask[c].has_value() = false;
            }
        }
        auto lmuNibble = lsmetrics::nibble(lmuIds, mask, std::ceil(std::sqrt(nCellArea)));
        lmuIds = lmuNibble;
    }

    int ProjectArea::getIndex(int n) {
        static int i = 0;
        static std::mutex mut;
        std::scoped_lock lock{ mut };
        auto out = i;
        i = i + n;
        return out;
    }

    void ProjectArea::divideLmusThread(int& sofar, std::mutex& mut, std::unordered_map<int, int>& regionArea, spatial::Raster<int>& lmus, spatial::Raster<int>& newlmus, const int thisThread) {
        int nLmu = regionArea.size();
        while (true) {
            mut.lock();
            int i = sofar;
            ++sofar;
            mut.unlock();
            if (i >= nLmu)
                break;
            int id = std::next(regionArea.begin(), i)->first;
            int area = std::next(regionArea.begin(), i)->second;
            int k = std::ceil((double)area / 650); // 650 cells = 150 acre unit which is max operational size according to jacob baker at stanislaus nf.
            if (k > 1) {
                std::cout << "Thread " + std::to_string(thisThread) + " starting lmu subdivision " + std::to_string(i) + "/" + std::to_string(nLmu) + "\n";
                auto kMeans = stats::Kmeans(k, 100);
                std::vector<stats::Kpoint> allPoints;
                for (spatial::cell_t j = 0; j < lmus.ncell(); ++j) {
                    if (lmus[j].has_value() && lmus[j].value() == id) {
                        allPoints.push_back(stats::Kpoint(lmus.xFromCell(j), lmus.yFromCell(j), j));
                    }
                }
                kMeans.run(allPoints);
                auto idx = getIndex(k);
                for (spatial::cell_t j = 0; j < allPoints.size(); ++j) {
                    newlmus[allPoints[j].cell].value() = idx + allPoints[j].clusterid;
                }
            }
            else {
                xt::filter(newlmus.values(), xt::equal(newlmus.values().value(), id)).value() = getIndex(1);
            }
        }
    }

    void ProjectArea::coreGapThread(spatial::Raster<int>& osiNum, spatial::Raster<int>& osiDen, spatial::Raster<int>& bbOsiNum, spatial::Raster<int>& bbOsiDen, const int nThread, const int thisThread, std::mutex& mut, int& sofar,
        const spatial::Raster<int>& maskr, const double canopycutoff, const double coregapdist, std::pair<spatial::coord_t, spatial::coord_t> expectedRes,
        dbhFunction dbhFunc, double bbDbh) {
        int ntile = lidarDataset->nTiles();

        spatial::coord_t chmres = spatial::Raster<int>(lidar::stringOrThrow(lidarDataset->getMaxHeightRaster(0))).xres() * lidarDataset->getConvFactor();
        while (true) {

            mut.lock();
            ++sofar;
            if (sofar >= ntile) {
                mut.unlock();
                break;
            }
            int i = sofar;
            mut.unlock();

            spatial::Raster<int> chm;
            spatial::Raster<int> basinMap;
            lico::TaoList taos;
            auto e = lidarDataset->extentByTile(i);
            if (!e.has_value()) {
                std::cerr << "Extent you tried to create from tile does not exist\n";
                throw std::runtime_error("No extent created.");
            }
            try {
                if (e.value().overlaps(maskr)) {
                    std::cout << "Tile " + std::to_string(i) + "/" + std::to_string(ntile) + " on thread " + std::to_string(thisThread) + "\n";
                    taos = lidarDataset->getHighPointsByTile(i, 3*lidarDataset->getConvFactor());
                    chm = spatial::Raster<int>(lidar::stringOrThrow(lidarDataset->getMaxHeightRaster(i)));
                    basinMap = spatial::Raster<int>(lidar::stringOrThrow(lidarDataset->getSegmentRaster(i)));
                    chm.repairResolution(expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);
                    basinMap.repairResolution(expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);
                }
                else {
                    continue;
                }
            }
            catch (lidar::FileNotFoundException e) {
                continue;
            }

            auto bbChm = spatial::Raster<int>((spatial::Alignment)chm);
            auto dbh = dbhFunc(taos.height());
            for (lico::index_t i = 0; i < taos.size(); ++i) {
                if (dbh[i] >= bbDbh) {
                    auto v = basinMap.extract(taos.x()[i], taos.y()[i]).value();
                    for (spatial::cell_t j = 0; j < basinMap.ncell(); j++) {
                        if (basinMap[j].has_value() && basinMap[j].value() == v) {
                            if (v == 1) {
                                if (e.value().contains(taos.x()[i], taos.y()[i]))
                                    throw std::runtime_error("yuck");
                                else
                                    continue;
                            }
                            bbChm[j].value() = 999;
                        }
                    }
                }
            }

            bbChm.values().has_value() = chm.values().has_value();

            chm.projection(maskr.projection());
            bbChm.projection(maskr.projection());
            spatial::Alignment thisalign{ maskr };
            thisalign.crop(chm, spatial::SnapType::out);
            thisalign.extend(chm, spatial::SnapType::out);
                   
            lsmetrics::crop_mw_function<int, int> numFunc = [&](const lsmetrics::crop_view<int>& e)->xtl::xoptional<int> {return lsmetrics::OSInumerator(e, chmres, canopycutoff, coregapdist); };
            spatial::Raster<int> thiscorenum = lsmetrics::movingWindowByRaster(chm, thisalign, numFunc, coregapdist);
            spatial::Raster<int> thisbbnum = lsmetrics::movingWindowByRaster(bbChm, thisalign, numFunc, coregapdist);

            lsmetrics::crop_mw_function<int, int> totalFunc = [&](const lsmetrics::crop_view<int>& e)->xtl::xoptional<int> {return lsmetrics::totalAreaForOSI(e, chmres, coregapdist); };
            spatial::Raster<int> thistotal = lsmetrics::movingWindowByRaster(chm, thisalign, totalFunc, coregapdist);
            spatial::Raster<int> thisbbtotal = lsmetrics::movingWindowByRaster(bbChm, thisalign, totalFunc, coregapdist);

            thiscorenum.crop(e.value(), spatial::SnapType::out);
            thistotal.crop(e.value(), spatial::SnapType::out);
            thisbbnum.crop(e.value(), spatial::SnapType::out);
            thisbbtotal.crop(e.value(), spatial::SnapType::out);
            
            std::vector<spatial::Raster<int>*> v = { &thiscorenum,&osiNum };
            std::vector<spatial::Raster<int>*> vtotal = { &thistotal,&osiDen };
            std::vector<spatial::Raster<int>*> vbb = { &thisbbnum,&bbOsiNum };
            std::vector<spatial::Raster<int>*> vbbtotal = { &thisbbtotal,&bbOsiDen };

            mut.lock();
            osiNum = spatial::rasterMergeInterior(v);
            osiDen = spatial::rasterMergeInterior(vtotal);
            bbOsiNum = spatial::rasterMergeInterior(vbb);
            bbOsiDen = spatial::rasterMergeInterior(vbbtotal);

            for (lico::index_t i = 0; i < taos.size(); ++i) {
                if (e.value().contains(taos.x()[i], taos.y()[i]))
                    allTaos.addTAO(taos[i]);
            }
            mut.unlock();
        }
    }

    void ProjectArea::postGapThread(spatial::Raster<int>& osiNum, spatial::Raster<int>& osiDen, const lico::TaoList& taos, const int nThread, const int thisThread, std::mutex& mut, int& sofar,
        const spatial::Raster<int>& maskr, const double canopycutoff, const double coregapdist, std::pair<spatial::coord_t, spatial::coord_t> expectedRes) {
        int ntile = lidarDataset->nTiles();

        auto l = spatial::SpVectorDataset<spatial::SpPolygon>(lidar::stringOrThrow(lidarDataset->getTileLayoutVector()));
        auto poFt = l.getFeaturesPtr();

        spatial::coord_t chmres = spatial::Raster<int>(lidar::stringOrThrow(lidarDataset->getMaxHeightRaster(0))).xres() * lidarDataset->getConvFactor();
        while (true) {
            mut.lock();
            ++sofar;
            if (sofar >= ntile) {
                mut.unlock();
                break;
            }
            int i = sofar;
            mut.unlock();

            spatial::Raster<int> basinMap;
            auto e = lidarDataset->extentByTile(i);
            if (!e.has_value()) {
                std::cerr << "Extent you tried to create from tile does not exist\n";
                throw std::runtime_error("No extent created.");
            }
            try {
                if (e.value().overlaps(maskr)) {
                    std::cout << "Tile " + std::to_string(i) + "/" + std::to_string(ntile) + " on thread " + std::to_string(thisThread) + "\n";
                    basinMap = spatial::Raster<int>(lidar::stringOrThrow(lidarDataset->getSegmentRaster(i)));
                    basinMap.repairResolution(expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);

                }
                else {
                    continue;
                }
            }
            catch (lidar::FileNotFoundException e) {
                continue;
            }

            spatial::Raster<int> chm{ (spatial::Alignment) basinMap };
            for (lico::index_t i = 0; i < taos.size(); ++i) {
                if (basinMap.contains(taos.x()[i], taos.y()[i])) {
                    auto v = basinMap.extract(taos.x()[i], taos.y()[i]);
                    if (v.has_value()) {
                        for (spatial::cell_t j = 0; j < basinMap.ncell(); j++) {
                            if (basinMap[j].has_value() && basinMap[j].value() == v.value()) {
                                if (v.value() == 1) {
                                    /*if (e.contains(taos.x()[i], taos.y()[i]))
                                        throw std::runtime_error("yuck");
                                    else*/
                                        continue;
                                }
                                chm[j].value() = 999;
                            }
                        }
                    }
                }
            }
            chm.values().has_value() = basinMap.values().has_value();
            chm.projection(maskr.projection());

            spatial::Alignment thisalign{ maskr };
            thisalign.crop(chm, spatial::SnapType::out);
            thisalign.extend(chm, spatial::SnapType::out);
            lsmetrics::crop_mw_function<int, int> numFunc = [&](const lsmetrics::crop_view<int>& e)->xtl::xoptional<int> {return lsmetrics::OSInumerator(e, chmres, canopycutoff, coregapdist); };
            spatial::Raster<int> thiscorenum = lsmetrics::movingWindowByRaster(chm, thisalign, numFunc, coregapdist);

            lsmetrics::crop_mw_function<int, int> totalFunc = [&](const lsmetrics::crop_view<int>& e)->xtl::xoptional<int> {return lsmetrics::totalAreaForOSI(e, chmres, coregapdist); };
            spatial::Raster<int> thistotal = lsmetrics::movingWindowByRaster(chm, thisalign, totalFunc, coregapdist);

            thiscorenum.crop(e.value(), spatial::SnapType::out);
            thistotal.crop(e.value(), spatial::SnapType::out);

            std::vector<spatial::Raster<int>*> v = { &thiscorenum,&osiNum };
            std::vector<spatial::Raster<int>*> vtotal = { &thistotal,&osiDen };

            mut.lock();
            osiNum = spatial::rasterMergeInterior(v);
            osiDen = spatial::rasterMergeInterior(vtotal);
            mut.unlock();
        }
    }
}