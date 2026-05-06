#pragma once

#ifndef rxtools_taolist_h
#define rxtools_taolist_h

#include "rxtools_pch.hpp"
#include "GraphLico.hpp"

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

    class TaoList {
    public:
        TaoList() = default;
        TaoList(const lapis::CoordRef& crs) : _xy(crs) {}
        
        template<class DATASET>
        TaoList(DATASET t, TaoGetters<DATASET> g);

        template<class DATASET>
        TaoList(std::string f, TaoGetters<DATASET>  g);

        const lapis::CoordRef& crs() const;
        void addTao(lapis::CoordXY xy, lapis::coord_t height, lapis::coord_t radius, lapis::coord_t area, lapis::coord_t dbh);
        const size_t size() const;

        const lapis::CoordXY xy(size_t i) const;
        const lapis::coord_t x(size_t i) const;
        const lapis::coord_t y(size_t i) const;
        const lapis::coord_t height(size_t i) const;
        const lapis::coord_t radius(size_t i) const;
        const lapis::coord_t area(size_t i) const;
        const double dbh(size_t i) const;

        const lapis::lico::TaoNode node(size_t i) const;

        void writeCsv(std::filesystem::path path) const;
        void writeShapefile(std::filesystem::path path) const;

    private:
        lapis::CoordXYVector _xy;
        std::vector<lapis::coord_t> _height;
        std::vector<lapis::coord_t> _radius;
        std::vector<lapis::coord_t> _area;
        std::vector<double> _dbh;
    };

    template<class DATASET>
    TaoList::TaoList(DATASET t, TaoGetters<DATASET> g) {
        for (auto f : t) {
            if (g.predicate(f)) {
                addTao(g.xy(f), g.height(f), g.radius(f), g.area(f), g.dbh(f));
            }
        }
    }

    template<class DATASET>
    TaoList::TaoList(std::string f, TaoGetters<DATASET>  g) {
        DATASET dataset(f);
        for (auto f : dataset) {
            if (g.predicate(f)) {
                addTao(g.xy(f), g.height(f), g.radius(f), g.area(f), g.dbh(f));
            }
        }
    }
} // namespace rxtools
#endif //rxtools_taolist_h