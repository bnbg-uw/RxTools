#pragma once

#ifndef rxtools_projectarea_h
#define rxtools_projectarea_h

#include "utilities.hpp"
#include "lmu.hpp"
#include "ProcessedFolder/src/readProcessedFolder.hpp"
#include "RasterAlgos.hpp"

namespace rxtools {
    class ProjectArea {
    public:
        std::unique_ptr<processedfolder::ProcessedFolder> lidarDataset;
        TaoListMP allTaos{};
        lapis::VectorDataset<lapis::MultiPolygon> projectPoly;

        lapis::Raster<double> aet;
        lapis::Raster<double> cwd;
        lapis::Raster<double> tmn;

        lapis::Raster<lapis::cell_t> lmuRaster;
        lapis::Raster<lapis::cell_t> lmuIds;
        lapis::Raster<int> osiDen;
        lapis::Raster<int> osiNum;
        lapis::Raster<int> bbOsiDen;
        lapis::Raster<int> bbOsiNum;

        std::unordered_map<lapis::cell_t, lapis::cell_t> regionType;

        ProjectArea() = default;
        //lmu param can be moderate or steep
        ProjectArea(std::string lidarDatasetPath, std::string projectPolygonPath, std::string aetPath, std::string cwdPath, std::string tmnPath, int nThread, std::string lmuPath = "", std::string lmuParam = "moderate");

        //terrain can be moderate or steep.
        //units in lidar units.
        lapis::Raster<lapis::cell_t> createLmuRasterFromTpiAndAsp(std::unique_ptr<processedfolder::ProcessedFolder>& lds,
                                                            std::string terrain = "moderate",
                                                            lapis::VectorDataset<lapis::MultiPolygon> projectPoly = lapis::VectorDataset<lapis::MultiPolygon>());
        void subdivideLmus(std::string climateClassPath, int nThread);
        int getIndex(int n);
        void divideLmusThread(int& sofar, std::mutex& mut, std::unordered_map<lapis::cell_t, int>& regionArea, lapis::Raster<lapis::cell_t>& lmus, lapis::Raster<lapis::cell_t>& newlmus, const int thisThread);
        void createCoreGapAndReadTaos(int nThread, double bbDbh, TaoGettersMP getters);

        Lmu createLmuThread(int& sofar, const int thisThread);
        void coreGapThread(lapis::Raster<int>& osiNum, lapis::Raster<int>& osiDen, lapis::Raster<int>& bbOsiNum, lapis::Raster<int>& bbOsiDen, const int nThread, const int thisThread, std::mutex& mut, int& sofar,
            const lapis::Raster<lapis::cell_t>& maskr, const double canopycutoff, double coregapdist, std::pair<lapis::coord_t, lapis::coord_t> expectedRes,
            TaoGettersMP getters, double bbDbh);
        void postGapThread(lapis::Raster<int>& osiNum, lapis::Raster<int>& osiDen, const TaoListMP& taos, const int nThread, const int thisThread, std::mutex& mut, int& sofar,
            const lapis::Raster<int>& maskr, const double canopycutoff, const double coregapdist, std::pair<lapis::coord_t, lapis::coord_t> expectedRes);

    };
}  // namespace rxtools

#endif  // !rxtools_projectarea_h