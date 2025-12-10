#pragma once

#ifndef rxtools_taolist_h
#define rxtools_taolist_h

#include "Vector.hpp"
#include "lico/src/GraphLico.hpp"
#include "Coordinate.hpp"

namespace rxtools {

    template<class DATASET>
    using TaoHeightGetter = std::function<lapis::coord_t(const typename DATASET::ConstFeatureType&)>;

    template<class DATASET>
    struct TaoGetters {
        lapis::lico::TaoPredicate<DATASET> predicate;
        lapis::lico::TaoCoordGetter<DATASET> xy;
        TaoHeightGetter<DATASET> height;
        lapis::lico::TaoRadiusGetter<DATASET> radius;
        lapis::lico::TaoAreaGetter<DATASET> area;
        lapis::lico::TaoDbhGetter<DATASET> dbh;
    };

    template<class DATASET>
    class TaoList {
    public:
        TaoGetters<DATASET> getters;
        DATASET taoVector;

        TaoList() = default;
        TaoList(DATASET t, TaoGetters<DATASET> g) : taoVector(t), getters(g) {}
        TaoList(std::string f, TaoGetters<DATASET>  g) : taoVector(lapis::VectorDataset<DATASET>(f)), getters(g) {}

        const typename DATASET::ConstFeatureType operator()(size_t i) {
            return taoVector.getFeature(i);
        }

        const lapis::lico::TaoNode nodeFromIndex(size_t i)  const {
            return lapis::lico::TaoNode(
                getters.xy(taoVector.getFeature(i)).x,
                getters.xy(taoVector.getFeature(i)).y,
                getters.radius(taoVector.getFeature(i)),
                getters.area(taoVector.getFeature(i)),
                1,
                getters.dbh(taoVector.getFeature(i)),
                nullptr
            );
        }

        const size_t size() const {
            return taoVector.nFeatures();
        }

        const lapis::CoordXY xy(size_t i) const {
            return getters.xy(taoVector.getFeature(i));
        }

        const lapis::coord_t x(size_t i) const {
            return xy(i).x;
        }

        const lapis::coord_t y(size_t i) const {
            return xy(i).y;
        }

        const lapis::coord_t height(size_t i) const {
            return getters.height(taoVector.getFeature(i))
        }

        const lapis::coord_t radius(size_t i) const {
            return getters.radius(taoVector.getFeature(i));
        }

        const lapis::coord_t area(size_t i) const {
            return getters.area(taoVector.getFeature(i));
        }

        const double dbh(size_t i) const  {
            return getters.dbh(taoVector.getFeature(i));
        }

        void writeCsv(std::filesystem::path path) const {
            std::ofstream out;
            out.open(path);

            out << "GridHighX,GridHighY,PolyArea,GridMaxHt,Radius,DBH\n";
            for (int i = 0; i < size(); i++) {
                out << std::setprecision(std::numeric_limits<double>::max_digits10) << x(i) << "," << y(i) << "," << area(i) << "," << height(i) << "," << radius(i) << "," << dbh(i) << "\n";
            }
            out.close();
        }
    };

    using TaoListMP = TaoList<lapis::VectorDataset<lapis::MultiPolygon>>;
    using TaoGettersMP = TaoGetters<lapis::VectorDataset<lapis::MultiPolygon>>;
} // namespace rxtools
#endif //rxtools_taolist_h