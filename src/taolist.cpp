#include "taolist.hpp"

namespace rxtools {
    const lapis::CoordRef& TaoList::crs() const {
        return _xy.crs;
    }

    void TaoList::addTao(lapis::CoordXY xy, lapis::coord_t height, lapis::coord_t radius, lapis::coord_t area, lapis::coord_t dbh) {
        _xy.push_back(xy);
        _height.push_back(height);
        _radius.push_back(radius);
        _area.push_back(area);
        _dbh.push_back(dbh);
    }

    const size_t TaoList::size() const {
        return _xy.size();
    }

    const lapis::CoordXY TaoList::xy(size_t i) const {
        return _xy.at(i);
    }

    const lapis::coord_t TaoList::x(size_t i) const {
        return xy(i).x;
    }

    const lapis::coord_t TaoList::y(size_t i) const {
        return xy(i).y;
    }

    const lapis::coord_t TaoList::height(size_t i) const {
        return _height.at(i);
    }

    const lapis::coord_t TaoList::radius(size_t i) const {
        return _radius.at(i);
    }

    const lapis::coord_t TaoList::area(size_t i) const {
        return _area.at(i);
    }

    const double TaoList::dbh(size_t i) const {
        return _dbh.at(i);
    }

    const lapis::lico::TaoNode TaoList::node(size_t i) const {
        return lapis::lico::TaoNode(x(i), y(i), radius(i), area(i), dbh(i));
    }

    void TaoList::writeCsv(std::filesystem::path path) const {
        std::ofstream out;
        out.open(path);

        out << "GridHighX,GridHighY,PolyArea,GridMaxHt,Radius,DBH\n" << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (int i = 0; i < size(); i++) {
            out << x(i) << "," << y(i) << "," << area(i) << "," << height(i) << "," << radius(i) << "," << dbh(i) << "\n";
        }
        out.close();
    }

    void TaoList::writeShapefile(std::filesystem::path path) const {
        lapis::VectorDataset<lapis::Point> pts(crs());
        pts.addNumericField<lapis::coord_t>("X");
        pts.addNumericField<lapis::coord_t>("Y");
        pts.addNumericField<lapis::coord_t>("Height");
        pts.addNumericField<lapis::coord_t>("Radius");
        pts.addNumericField<lapis::coord_t>("Area");
        pts.addRealField("DBH");
        for (int i = 0; i < size(); i++) {
            lapis::Point p(xy(i), crs());
            pts.addGeometry(p);

            pts.back().setNumericField("X", x(i));
            pts.back().setNumericField("Y", y(i));
            pts.back().setNumericField("Height", height(i));
            pts.back().setNumericField("Radius", radius(i));
            pts.back().setNumericField("Area", area(i));
            pts.back().setRealField("DBH", dbh(i));
        }
        pts.writeShapefile(path);
    }
}