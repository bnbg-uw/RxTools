#pragma once

#ifndef licosim_projectarea_h
#define licosim_projectarea_h

#include "licosim/lmu.hpp"
#include "shapefile/shapefile.h"
#include <set>

namespace licosim {
class ProjectArea {
public:
    std::unique_ptr<lidar::ProcessedFolder> lidarDataset;
    lico::TaoList allTaos{};
    spatial::SpVectorDataset<spatial::SpMultiPolygon> projectPoly;

    spatial::Raster<double> aet;
    spatial::Raster<double> cwd;
    spatial::Raster<double> tmn;

    spatial::Raster<int> lmuRaster;
    spatial::Raster<int> lmuIds;
    spatial::Raster<int> osiDen;
    spatial::Raster<int> osiNum;
    spatial::Raster<int> bbOsiDen;
    spatial::Raster<int> bbOsiNum;

    std::unordered_map<int, int> regionType;

    ProjectArea() = default;
    //lmu param can be moderate or steep
    ProjectArea(std::string lidarDatasetPath, std::string projectPolygonPath, std::string aetPath, std::string cwdPath, std::string tmnPath, int nThread, std::string lmuPath = "", std::string lmuParam = "moderate");

    //terrain can be moderate or steep.
    //units in lidar units.
    spatial::Raster<int> createLmuRasterFromTpiAndAsp(std::unique_ptr<lidar::ProcessedFolder>& lds,
                                                        std::string terrain = "moderate",
                                                        spatial::SpVectorDataset<spatial::SpMultiPolygon> projectPoly = spatial::SpVectorDataset<spatial::SpMultiPolygon>());
    void subdivideLmus(std::string climateClassPath, int nThread);
    int getIndex(int n);
    void divideLmusThread(int& sofar, std::mutex& mut, std::unordered_map<int, int>& regionArea, spatial::Raster<int>& lmus, spatial::Raster<int>& newlmus, const int thisThread);
    void createCoreGapAndReadTaos(int nThread, double bbDbh, dbhFunction dbhFunc);
    Lmu createLmuThread(int& sofar, const int thisThread);
    void coreGapThread(spatial::Raster<int>& osiNum, spatial::Raster<int>& osiDen, spatial::Raster<int>& bbOsiNum, spatial::Raster<int>& bbOsiDen, const int nThread, const int thisThread, std::mutex& mut, int& sofar,
        const spatial::Raster<int>& maskr, const double canopycutoff, const double coregapdist, std::pair<spatial::coord_t, spatial::coord_t> expectedRes,
        dbhFunction dbhFunc, double bbDbh);
    void postGapThread(spatial::Raster<int>& osiNum, spatial::Raster<int>& osiDen, const lico::TaoList& taos, const int nThread, const int thisThread, std::mutex& mut, int& sofar,
        const spatial::Raster<int>& maskr, const double canopycutoff, const double coregapdist, std::pair<spatial::coord_t, spatial::coord_t> expectedRes);

};
}  // namespace licosim

#endif  // !licosim_projectarea_h