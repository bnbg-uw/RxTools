#pragma once

#ifndef licosim_allometry_h
#define licosim_allometry_h

#include "Shapefile/shapefile.h"
#include "LICO/LICO.hpp"
#include "licosim/utilities.hpp"
#include "OSGeo/gdal_priv.h"
#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include <math.h>

namespace licosim {

class Allometry {
public:
    std::string fiaPath = "";
    stats::FIALeastSquares model;

    Allometry() {};
    Allometry(std::string path) : fiaPath(path) {
        auto allPlots = stats::FIALeastSquares::getPlotList(fiaPath);
        model = stats::FIALeastSquares(allPlots, fiaPath, "DIA");
    }
    Allometry(std::vector<double> coeff) {
        model.intercept = coeff[0];
        model.slope = coeff[1];
        model.responseName = "DIA";
        model.init = true;
    }

    std::vector<double> getDbhFromHeightLinear(const std::vector<double>& ht) const;
    std::vector<double> getDbhFromHeightLinear(lico::adapt_type<spatial::unit_t> ht) const;
    inline spatial::unit_t getDbhFromHeightLinear(spatial::unit_t ht) const {
        return model.intercept + model.slope * ht;
    }

    std::vector<double> getDbhFromHeightAuto(const std::vector<double>& ht);
    std::vector<double> getDbhFromHeightAuto(lico::adapt_type<spatial::unit_t> ht);
    inline spatial::unit_t getDbhFromHeightAuto(spatial::unit_t ht) {
        return model.predictAsMetric(ht);
    }
};

} // namespace licosim

#endif // !licosim_allometry_h