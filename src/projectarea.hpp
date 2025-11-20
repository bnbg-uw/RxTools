#pragma once

#ifndef rxtools_projectarea_h
#define rxtools_projectarea_h

#include "lmu.hpp"
#include "ProcessedFolder/src/readProcessedFolder.hpp"

namespace rxtools {
    class ProjectArea {
    public:
        std::unique_ptr<processedfolder::ProcessedFolder> lidarDataset;
        TaoPointList allTaos{};
        lapis::VectorDataset<lapis::MultiPolygon> projectPoly;

        lapis::Raster<double> aet;
        lapis::Raster<double> cwd;
        lapis::Raster<double> tmn;

        lapis::Raster<int> lmuRaster;
        lapis::Raster<int> lmuIds;
        lapis::Raster<int> osiDen;
        lapis::Raster<int> osiNum;
        lapis::Raster<int> bbOsiDen;
        lapis::Raster<int> bbOsiNum;

        std::unordered_map<int, int> regionType;

        ProjectArea() = default;
        //lmu param can be moderate or steep
        ProjectArea(std::string lidarDatasetPath, std::string projectPolygonPath, std::string aetPath, std::string cwdPath, std::string tmnPath, int nThread, std::string lmuPath = "", std::string lmuParam = "moderate");

        //terrain can be moderate or steep.
        //units in lidar units.
        lapis::Raster<int> createLmuRasterFromTpiAndAsp(std::unique_ptr<processedfolder::ProcessedFolder>& lds,
                                                            std::string terrain = "moderate",
                                                            lapis::VectorDataset<lapis::MultiPolygon> projectPoly = lapis::VectorDataset<lapis::MultiPolygon>());
        void subdivideLmus(std::string climateClassPath, int nThread);
        int getIndex(int n);
        void divideLmusThread(int& sofar, std::mutex& mut, std::unordered_map<int, int>& regionArea, lapis::Raster<int>& lmus, lapis::Raster<int>& newlmus, const int thisThread);
        void createCoreGapAndReadTaos(int nThread, double bbDbh, TaoPointGetters getters);

        Lmu createLmuThread(int& sofar, const int thisThread);
        void coreGapThread(lapis::Raster<int>& osiNum, lapis::Raster<int>& osiDen, lapis::Raster<int>& bbOsiNum, lapis::Raster<int>& bbOsiDen, const int nThread, const int thisThread, std::mutex& mut, int& sofar,
            const lapis::Raster<int>& maskr, const double canopycutoff, const double coregapdist, std::pair<lapis::coord_t, lapis::coord_t> expectedRes,
            TaoPointGetters getters, double bbDbh);
        void postGapThread(lapis::Raster<int>& osiNum, lapis::Raster<int>& osiDen, const TaoPointList& taos, const int nThread, const int thisThread, std::mutex& mut, int& sofar,
            const lapis::Raster<int>& maskr, const double canopycutoff, const double coregapdist, std::pair<lapis::coord_t, lapis::coord_t> expectedRes);

    };
}  // namespace rxtools

#endif  // !rxtools_projectarea_h