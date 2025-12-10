#include "projectarea.hpp"

namespace rxtools {
    ProjectArea::ProjectArea(std::string lidarDatasetPath, std::string projectPolygonPath, std::string aetPath, std::string cwdPath, std::string tmnPath, int nThread, std::string lmuPath, std::string lmuParam) {
        lidarDataset = processedfolder::readProcessedFolder(lidarDatasetPath);

        projectPoly = lapis::VectorDataset<lapis::MultiPolygon>(projectPolygonPath);
        if (!projectPoly.crs().isConsistent(lidarDataset->crs()))
            projectPoly.projectInPlacePreciseExtent(lidarDataset->crs());


        auto mask = lapis::Raster<int>(processedfolder::stringOrThrow(lidarDataset->maskRaster()));
        mask = lapis::cropRaster(mask, projectPoly.extent(), lapis::SnapType::out);

        std::cout << "\t\tReading and resampling climate layers...";
        aet = lapis::resampleRaster(lapis::Raster<double>(aetPath), mask, lapis::ExtractMethod::bilinear);
        cwd = lapis::resampleRaster(lapis::Raster<double>(cwdPath), mask, lapis::ExtractMethod::bilinear);
        tmn = lapis::resampleRaster(lapis::Raster<double>(tmnPath), mask, lapis::ExtractMethod::bilinear);
        std::cout << tmnPath << "\n";

        if (lmuPath.empty()) {
            std::cout << "\t\tCreating LMU's:\n";
            lmuRaster = createLmuRasterFromTpiAndAsp(lidarDataset, lmuParam, projectPoly);
            lmuIds = lapis::connectedComponents(lmuRaster, false);
            std::cout << "\t\tLMU creation done!\n";
            //lmuRaster.writeRaster("E:/yubalmus.img");
        }
        else {
            // check extents etc...
            std::cout << "\t\tReading LMU raster...";
            lmuRaster = lapis::Raster<lapis::cell_t>(lmuPath);
            if (!lmuRaster.crs().isConsistent(mask.crs()))
                lmuRaster = lapis::resampleRaster(lmuRaster, mask, lapis::ExtractMethod::near);
            else {
                lmuRaster = lapis::extendRaster(lmuRaster, mask, lapis::SnapType::out);
                lmuRaster = lapis::cropRaster(lmuRaster, mask, lapis::SnapType::out);
            }
            lmuRaster.mask(mask);
            lmuRaster = lapis::trimRaster(lmuRaster);
            lmuIds = lapis::connectedComponents(lmuRaster, false);
            std::cout << " Done!\n";
        }

        osiNum = lapis::Raster<int>{ (lapis::Alignment)lmuIds };
        osiDen = lapis::Raster<int>{ (lapis::Alignment)lmuIds };
        bbOsiNum = lapis::Raster<int>{ (lapis::Alignment)lmuIds };
        bbOsiDen = lapis::Raster<int>{ (lapis::Alignment)lmuIds };

        std::cout << "\t\tCreating regiontype map..\n";
        for (lapis::cell_t c = 0; c < lmuIds.ncell(); ++c) {
            if (lmuIds[c].has_value()) {
                regionType.emplace(lmuIds[c].value(), lmuRaster[c].value());
            }
        }

        std::cout << "\t\tProject Area Construction done.\n";
    }

    void ProjectArea::createCoreGapAndReadTaos(int nThread, double bbDbh, TaoGettersMP getters) {
        std::cout << "\t\tcalculating pretreat OSI and reading taos...\n";
        std::pair<lapis::coord_t, lapis::coord_t> expectedRes{};
        if (lidarDataset->units() == lapis::linearUnitPresets::meter) {
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
                coreGapThread(osiNum, osiDen, bbOsiNum, bbOsiDen, nThread, i, mut, sofar, lmuRaster, 2, 6, expectedRes, getters, bbDbh);
            }
            catch (processedfolder::FileNotFoundException e) {
                std::cerr << e.what() << '\n';
                std::cerr << "Aborting\n";
                throw e;
            }
            catch (lapis::InvalidRasterFileException e) {
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

        osiNum = lapis::extendRaster(osiNum, lmuRaster, lapis::SnapType::out);
        osiDen = lapis::extendRaster(osiDen, lmuRaster, lapis::SnapType::out);
        osiNum = lapis::cropRaster(osiNum, lmuRaster, lapis::SnapType::out);
        osiDen = lapis::cropRaster(osiDen, lmuRaster, lapis::SnapType::out);
        osiNum.mask(lmuRaster);
        osiDen.mask(lmuRaster);

        bbOsiNum = lapis::extendRaster(bbOsiNum, lmuRaster, lapis::SnapType::out);
        bbOsiDen = lapis::extendRaster(bbOsiDen, lmuRaster, lapis::SnapType::out);
        bbOsiNum = lapis::cropRaster(bbOsiNum, lmuRaster, lapis::SnapType::out);
        bbOsiDen = lapis::cropRaster(bbOsiDen, lmuRaster, lapis::SnapType::out);
        bbOsiNum.mask(lmuRaster);
        bbOsiDen.mask(lmuRaster);
    }
    
    Lmu ProjectArea::createLmuThread(int& sofar, const int thisThread) {
        int id = std::next(regionType.begin(), sofar)->first;
        int type = std::next(regionType.begin(), sofar)->second;
        
        std::cout << "\t Creating lmu " + std::to_string(sofar) + "/" + std::to_string(regionType.size()) + " " + std::to_string(id) + " on thread " + std::to_string(thisThread) + "\n";
        auto r = lmuIds;
        for (lapis::cell_t c = 0; c < r.ncell(); ++c) {
            if (r[c].value() != id) {
                r[c].has_value() = false;
            }
        }
        r = lapis::trimRaster(r);
        return Lmu(r, static_cast<LmuType>(type));
    }

    lapis::Raster<lapis::cell_t> ProjectArea::createLmuRasterFromTpiAndAsp(std::unique_ptr<processedfolder::ProcessedFolder>& lds,
                                                                   std::string terrain,
                                                                   lapis::VectorDataset<lapis::MultiPolygon> projectPoly) {
        double ridgeSep, canyonSep;

        int expSize = 0;
        if (terrain == "steep") {
            ridgeSep = 25; //ridges
            canyonSep = 20; //canyons
        }
        else {
            ridgeSep = 15;
            canyonSep = 14;
        }

        double ridgeArea = 20000;
        double canyonArea = 40000;
        auto convFactor = lapis::linearUnitPresets::meter.convertOneToThis(1, lds->units());
        if (lds->units() != lapis::linearUnitPresets::meter) {
            ridgeArea /= (convFactor * convFactor);
            canyonArea /= (convFactor * convFactor);
        }

        std::cout << " TPI ";
        int tpiDist = 500;
        int aspDist = 135;
        // create ridge groups
        auto tpi = lapis::Raster<double>(processedfolder::stringOrThrow(lds->tpi(tpiDist, lapis::linearUnitPresets::meter)));
        if (lds->type() == processedfolder::RunType::fusion) {
            tpi = tpi / 100; do i need convfactor
        }
        if (projectPoly.nFeature()) {
            tpi = lapis::cropRaster(tpi, projectPoly.extent(), lapis::SnapType::out);
            //tpi.mask(projectPoly);
        }

        auto aspect = lapis::Raster<double>(processedfolder::stringOrThrow(lds->aspect(aspDist, lapis::linearUnitPresets::meter)));
        if (projectPoly.nFeature()) {
            aspect = lapis::cropRaster(aspect, projectPoly.extent(), lapis::SnapType::out);
            //aspect.mask(projectPoly);
        }
        aspect.mask(tpi);
        tpi.mask(aspect);

        auto ridge = tpi > ridgeSep;
        for (lapis::cell_t c = 0; c < ridge.ncell(); ++c) {
            if (!ridge[c].value()) {
                ridge[c].has_value() = false;
            }
        }
        auto ridgeGroup = lapis::connectedComponents(ridge, false);
        // Get unique names of ridges, then remove ridges smaller than min patch size.
        double cellArea = tpi.xres() * tpi.yres();
        int nCellArea = ridgeArea / cellArea;
        std::unordered_map<int, int> regionArea;
        for (lapis::cell_t c = 0; c < ridgeGroup.ncell(); ++c) {
            if (ridgeGroup[c].has_value()) {
                regionArea.emplace(ridgeGroup[c].value(), 0);
                ++regionArea[ridgeGroup[c].value()];
            }
        }
        for (lapis::cell_t c = 0; c < ridgeGroup.ncell(); ++c) {
            if (regionArea[ridgeGroup[c].value()] < nCellArea) {
                ridgeGroup[c].has_value() = false;
            }
        }
        std::cout << "\t\t\tRidges identified.\n";

        //create canyon groups
        auto canyon = tpi < -canyonSep;
        for (lapis::cell_t c = 0; c < canyon.ncell(); ++c) {
            if (!ridge[c].value()) {
                ridge[c].has_value() = false;
            }
        }
        auto canyonGroup = lapis::connectedComponents(canyon, false);
        //get unique names of ridge the remove ridge
        cellArea = canyon.xres() * canyon.yres();
        nCellArea = canyonArea / cellArea;
        regionArea.clear();
        for (auto c = 0; c < canyonGroup.ncell(); ++c) {
            if (canyonGroup[c].has_value()) {
                regionArea.emplace(canyonGroup[c].value(), 0);
                ++regionArea[canyonGroup[c].value()];
            }
        }
        for (lapis::cell_t c = 0; c < canyonGroup.ncell(); ++c) {
            if (regionArea[canyonGroup[c].value()] < nCellArea) {
                canyonGroup[c].has_value() = false;
            }
        }
        std::cout << "\t\t\tValley bottoms identified.\n";


        // Calculate slope position
        lapis::Raster<int> spos{ (lapis::Alignment)tpi };
        for (lapis::cell_t i = 0; i < spos.ncell(); ++i) {
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
        lapis::Raster<int> aspectClassified{ (lapis::Alignment)aspect };
        for (lapis::cell_t i = 0; i < aspect.ncell(); ++i) {
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
        lapis::Raster<lapis::cell_t> lmu{ (lapis::Alignment)aspect };
        for (lapis::cell_t i = 0; i < lmu.ncell(); ++i) {
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
        auto lmuGrp = lapis::connectedComponents(lmu, false);
        cellArea = canyon.xres() * canyon.yres();
        nCellArea = canyonArea / cellArea; // canyon area cutoff ~10acres, which is our lmu cutoff.
        regionArea.clear();
        for (lapis::cell_t c = 0; c < lmuGrp.ncell(); ++c) {
            if (lmuGrp[c].has_value()) {
                regionArea.emplace(lmuGrp[c].value(), 0);
                ++regionArea[lmuGrp[c].value()];
            }
        }
        for (lapis::cell_t c = 0; c < lmuGrp.ncell(); ++c) {
            if (regionArea[lmuGrp[c].value()] < nCellArea) {
                lmuGrp[c].has_value() = false;
            }
        }
        std::cout << "\t\t\tLarge lmu's identified (connected components created). Nibbling small lmu's into larger lmu's\n";

        //max dist needed to look would be radius of the circle representing the largest removable LMU.
        //a = pi*r^2, I'm not dividing by pi to leave buffer.
        auto lmuNibble = lapis::nibble(lmu, lmuGrp);
        std::cout << "\t\t\tNibbling done.\n";

        return lmuNibble;
    }

    void ProjectArea::subdivideLmus(std::string climateClassPath, int nThread) {
        auto cc = lapis::Raster<int>(climateClassPath);
        cc = lapis::resampleRaster(cc, lmuIds, lapis::ExtractMethod::near);
        lmuIds = lmuIds * 100;
        lmuIds = lmuIds + cc;
        std::cout << "climate done \n";

        std::unordered_map<int, int> regionArea;
        for (lapis::cell_t c = 0; c < lmuIds.ncell(); ++c) {
            if (lmuIds[c].has_value()) {
                regionArea.emplace(lmuIds[c].value(), 0);
                ++regionArea[lmuIds[c].value()];
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
        for (lapis::cell_t c = 0; c < lmuIds.ncell(); ++c) {
            if (lmuIds[c].has_value()) {
                regionType.emplace(lmuIds[c].value(), lmuRaster[c].value());
            }
        }
        auto nCellArea = 300; // canyon area cutoff ~10acres, which is our lmu cutoff.
        regionArea.clear();
        auto mask = newlmus;
        for (lapis::cell_t c = 0; c < newlmus.ncell(); ++c) {
            if (newlmus[c].has_value()) {
                regionArea.emplace(newlmus[c].value(), 0);
                ++regionArea[newlmus[c].value()];
            }
        }
        for (lapis::cell_t c = 0; c < newlmus.ncell(); ++c) {
            if (regionArea[newlmus[c].value()] < nCellArea) {
                mask[c].has_value() = false;
            }
        }
        auto lmuNibble = lapis::nibble(lmuIds, mask);
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

    void ProjectArea::divideLmusThread(int& sofar, std::mutex& mut, std::unordered_map<int, int>& regionArea, lapis::Raster<lapis::cell_t>& lmus, lapis::Raster<lapis::cell_t>& newlmus, const int thisThread) {
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
                auto kMeans = utilities::Kmeans(k, 100);
                std::vector<utilities::Kpoint> allPoints;
                for (lapis::cell_t j = 0; j < lmus.ncell(); ++j) {
                    if (lmus[j].has_value() && lmus[j].value() == id) {
                        allPoints.push_back(utilities::Kpoint(lmus.xFromCell(j), lmus.yFromCell(j), j));
                    }
                }
                kMeans.run(allPoints);
                auto idx = getIndex(k);
                for (lapis::cell_t j = 0; j < allPoints.size(); ++j) {
                    newlmus[allPoints[j].cell].value() = idx + allPoints[j].clusterid;
                }
            }
            else {
                for (lapis::cell_t c = 0; c < newlmus.ncell(); ++c) {
                    if (newlmus[c].value() == id) {
                        newlmus[c].value() = getIndex(1);
                    }
                }
            }
        }
    }

    void ProjectArea::coreGapThread(lapis::Raster<int>& osiNum, lapis::Raster<int>& osiDen, lapis::Raster<int>& bbOsiNum, lapis::Raster<int>& bbOsiDen, const int nThread, const int thisThread, std::mutex& mut, int& sofar,
        const lapis::Raster<lapis::cell_t>& maskr, const double canopycutoff, double coregapdist, std::pair<lapis::coord_t, lapis::coord_t> expectedRes,
        TaoGettersMP getters, double bbDbh) {

        int ntile = lidarDataset->nTiles();
        coregapdist = lidarDataset->units()->convertOneToThis(coregapdist, lapis::linearUnitPresets::meter);
        while (true) {

            mut.lock();
            ++sofar;
            if (sofar >= ntile) {
                mut.unlock();
                break;
            }
            int i = sofar;
            mut.unlock();

            lapis::Raster<lapis::coord_t> chm;
            lapis::Raster<int> basinMap;
            TaoListMP taos;
            auto e = lidarDataset->extentByTile(i);
            if (!e.has_value()) {
                std::cerr << "Extent you tried to create from tile does not exist\n";
                throw std::runtime_error("No extent created.");
            }
            try {
                if (e.value().overlaps(maskr)) {
                    auto eBuffer = lapis::bufferExtent(e.value(), coregapdist);
                    std::cout << "Tile " + std::to_string(i) + "/" + std::to_string(ntile) + " on thread " + std::to_string(thisThread) + "\n";
                    taos = TaoListMP(lidarDataset->polygons(eBuffer), getters);
                    chm = lidarDataset->csmRaster(eBuffer).value();
                    basinMap = lidarDataset->watershedSegmentRaster(eBuffer).value();

                    chm = processedfolder::repairResolution(chm, expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);
                    basinMap = processedfolder::repairResolution(basinMap, expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);
                }
                else {
                    continue;
                }
            }
            catch (processedfolder::FileNotFoundException e) {
                continue;
            }

            auto bbChm = lapis::Raster<int>((lapis::Alignment)chm);
            for (size_t i = 0; i < taos.size(); ++i) {
                if (taos.dbh(i) >= bbDbh) {
                    auto v = basinMap.extract(taos.x(i), taos.y(i), lapis::ExtractMethod::near).value();
                    for (lapis::cell_t j = 0; j < basinMap.ncell(); j++) {
                        if (basinMap[j].has_value() && basinMap[j].value() == v) {
                            if (v == 1) {
                                if (e.value().contains(taos.x(i), taos.y(i)))
                                    throw std::runtime_error("yuck");
                                else
                                    continue;
                            }
                            bbChm[j].value() = 999;
                        }
                    }
                }
            }

            for (lapis::cell_t c = 0; c < bbChm.ncell(); ++c) {
                bbChm[c].has_value() = chm[c].has_value();
            }

            chm.defineCRS(maskr.crs());
            bbChm.defineCRS(maskr.crs());
            lapis::Alignment thisalign{ maskr };
            thisalign = lapis::cropAlignment(thisalign, chm, lapis::SnapType::out);
            thisalign = lapis::extendAlignment(thisalign, chm, lapis::SnapType::out);
            
            lapis::Raster<lapis::coord_t> edt = lapis::euclideanDistanceTransform(chm);
            lapis::Raster<char> isCoreGap = edt >= 6;
            lapis::Raster<char> thiscorenum = lapis::aggregateSum(isCoreGap, thisalign);
            lapis::Raster<char> thisden = lapis::aggregateCount(isCoreGap, thisalign);

            edt = lapis::euclideanDistanceTransform(bbChm);
            isCoreGap = edt >= 6;
            lapis::Raster<char> thisbbnum = lapis::aggregateSum(isCoreGap, thisalign);
            lapis::Raster<char> thisbbden = lapis::aggregateCount(isCoreGap, thisalign);

            thiscorenum = lapis::cropRaster(thiscorenum, e.value(), lapis::SnapType::out);
            thisden = lapis::cropRaster(thisden, e.value(), lapis::SnapType::out);
            thisbbnum = lapis::cropRaster(thisbbnum, e.value(), lapis::SnapType::out);
            thisbbden = lapis::cropRaster(thisbbden, e.value(), lapis::SnapType::out);
            
            std::vector<lapis::Raster<int>*> v = { &thiscorenum,&osiNum };
            std::vector<lapis::Raster<int>*> vtotal = { &thisden,&osiDen };
            std::vector<lapis::Raster<int>*> vbb = { &thisbbnum,&bbOsiNum };
            std::vector<lapis::Raster<int>*> vbbtotal = { &thisbbden,&bbOsiDen };
            mut.lock();

            osiNum = lapis::mosaicInside(v);
            osiDen = lapis::mosaicInside(vtotal);
            bbOsiNum = lapis::mosaicInside(vbb);
            bbOsiDen = lapis::mosaicInside(vbbtotal);

            for (size_t i = 0; i < taos.size(); ++i) {
                if (e.value().contains(taos.x(i), taos.y(i)))
                    allTaos.taoVector.addFeature(taos.taoVector.getFeature(i));
            }
            mut.unlock();
        }
    }



    void ProjectArea::postGapThread(lapis::Raster<int>& osiNum, lapis::Raster<int>& osiDen, const TaoListMP& taos, const int nThread, const int thisThread, std::mutex& mut, int& sofar,
        const lapis::Raster<int>& maskr, const double canopycutoff, const double coregapdist, std::pair<lapis::coord_t, lapis::coord_t> expectedRes) {
        int ntile = lidarDataset->nTiles();

        auto l = lapis::VectorDataset<lapis::Polygon>(processedfolder::stringOrThrow(lidarDataset->tileLayoutVector()));

        while (true) {
            mut.lock();
            ++sofar;
            if (sofar >= ntile) {
                mut.unlock();
                break;
            }
            int i = sofar;
            mut.unlock();

            lapis::Raster<int> basinMap;
            auto e = lidarDataset->extentByTile(i);
            if (!e.has_value()) {
                std::cerr << "Extent you tried to create from tile does not exist\n";
                throw std::runtime_error("No extent created.");
            }
            try {
                if (e.value().overlaps(maskr)) {
                    auto eBuffer = lapis::bufferExtent(e.value(), coregapdist);
                    std::cout << "Tile " + std::to_string(i) + "/" + std::to_string(ntile) + " on thread " + std::to_string(thisThread) + "\n";
                    basinMap = lidarDataset->watershedSegmentRaster(eBuffer).value();
                    basinMap = processedfolder::repairResolution(basinMap, expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);
                }
                else {
                    continue;
                }
            }
            catch (processedfolder::FileNotFoundException e) {
                continue;
            }

            lapis::Raster<int> chm{ (lapis::Alignment) basinMap };
            for (size_t i = 0; i < taos.size(); ++i) {
                if (basinMap.contains(taos.x(i), taos.y(i))) {
                    auto v = basinMap.extract(taos.x(i), taos.y(i), lapis::ExtractMethod::near);
                    if (v.has_value()) {
                        for (lapis::cell_t j = 0; j < basinMap.ncell(); j++) {
                            if (basinMap[j].has_value() && basinMap[j].value() == v.value()) {
                                if (v.value() == 1) {
                                    continue;
                                }
                                chm[j].value() = 999;
                            }
                        }
                    }
                }
            }
            for (lapis::cell_t c = 0; c < chm.ncell(); ++c) {
                chm[c].has_value() = basinMap[c].has_value();
            }
            chm.defineCRS(maskr.crs());

            //lapis::coord_t chmres = lapis::Raster<int>(processedfolder::stringOrThrow(lidarDataset->maxHeightRaster(0))).xres() * lidarDataset->getConvFactor();

            lapis::Alignment thisalign{ maskr };
            thisalign = lapis::cropAlignment(thisalign, chm, lapis::SnapType::out);
            thisalign = lapis::extendAlignment(thisalign, chm, lapis::SnapType::out);
            
            lapis::Raster<lapis::coord_t> edt = lapis::euclideanDistanceTransform(chm);
            lapis::Raster<char> isCoreGap = edt >= 6;
            lapis::Raster<char> thiscorenum = lapis::aggregateSum(isCoreGap, thisalign);
            lapis::Raster<char> thisden = lapis::aggregateCount(isCoreGap, thisalign);

            thiscorenum = lapis::cropRaster(thiscorenum, e.value(), lapis::SnapType::out);
            thisden = lapis::cropRaster(thisden, e.value(), lapis::SnapType::out);

            std::vector<lapis::Raster<int>*> v = { &thiscorenum,&osiNum };
            std::vector<lapis::Raster<int>*> vtotal = { &thisden,&osiDen };

            mut.lock();
            osiNum = lapis::mosaicInside(v);
            osiDen = lapis::mosaicInside(vtotal);
            mut.unlock();
        }
    }
}