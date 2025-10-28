
#include "licosim/allometry.hpp"

namespace licosim {
    std::vector<double> Allometry::getDbhFromHeightLinear(const std::vector<double>& ht) const {
        std::vector<double> out;
        for (auto n : ht) {
            out.push_back(model.intercept + model.slope * n);
        }
        return out;
    }

    std::vector<double> Allometry::getDbhFromHeightLinear(lico::adapt_type<spatial::unit_t> ht) const {
        std::vector<double> out;
        for (int i = 0; i < ht.size(); i++) {
            out.push_back(model.intercept + model.slope * ht[i]);
        }
        return out;
    }

    std::vector<double> Allometry::getDbhFromHeightAuto(const std::vector<double>& ht) {
        return model.predictAsMetricMulti(ht);
    }

    std::vector<double> Allometry::getDbhFromHeightAuto(lico::adapt_type<spatial::unit_t> ht) {
        return model.predictAsMetricMulti(ht);

    }
} //namespace