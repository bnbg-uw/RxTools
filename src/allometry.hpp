#pragma once

#ifndef rxtools_allometry_h
#define rxtools_allometry_h

#include "rxtools_pch.hpp"
#include "utilities.hpp"
#include "Raster.hpp"

namespace rxtools::linearUnitPresets {
    const lapis::LinearUnit inch{ "inch",0.0254 };
    const lapis::LinearUnit centimeter{ "centimeter",0.01 };
}

namespace rxtools::allometry {
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;

    //string is the name of the plot, pair is x,y coords of plot.
    //TODO: Consider if PlotList can be a std::vector<std::pair<bg::point<...>, std::string>
    //      this way it could be director inserted into boost rtree rather than constructing that which is the current method.
    //      this would optimize calcKNNTree at the cost of potentially having duplicates in the plotlist.
    //      also plotTreeMap could use the pair as the key? fewer objects at the cost of complexity.
    using PlotList = std::unordered_map<std::string, lapis::CoordXY>;

    //Allows for more than one explanatory variable
    struct FIATreeList {
        std::vector<double> height;

        std::vector<std::string> names;
        std::vector<std::vector<double>> otherfields;
    };

    //abstract base class for various types allometric models we can run.
    class Model {
    public:
        lapis::LinearUnit inputUnit = lapis::linearUnitPresets::unknownLinear; //units the model expects as input
        lapis::LinearUnit outputUnit = lapis::linearUnitPresets::unknownLinear; //unitis the model expects as output

        virtual double predict(double x, const lapis::LinearUnit& thisUnit, const lapis::LinearUnit& returnUnit = lapis::linearUnitPresets::meter) const = 0;
        std::vector<double> predict(const std::vector<double>& x, const lapis::LinearUnit& thisUnit, const lapis::LinearUnit& returnUnit = lapis::linearUnitPresets::meter);

        virtual ~Model() = default;

        friend std::ostream& operator<<(std::ostream& os, const Model& m);

    protected:
        //override to make cout behavior work- pipe the output you'd like into os.
        virtual void print(std::ostream& os) const = 0;
    };

    inline std::ostream& operator<<(std::ostream& os, const Model& m) {
        m.print(os);
        return(os);
    }

    class FIAReader {
    public:
        PlotList plots;
        std::unordered_map<std::string, FIATreeList> plotTreeMap;

        using PlotPoint = bg::model::point<double, 2, bg::cs::cartesian>;
        using PlotPointName = std::pair<PlotPoint, std::string>;
        bgi::rtree<PlotPointName, bgi::quadratic<16>> tree;

        FIAReader(const std::string& fiaFolder);
        FIAReader(const std::string& fiaFolder, const lapis::Extent& e);

        const lapis::CoordRef projection();
        void projection(const lapis::CoordRef newCrs);
        void project(const lapis::CoordRef& newCrs);

        const std::size_t limitByExtent(const lapis::Extent& e);

        template<class T>
        const std::size_t limitByRasterValue(const lapis::Raster<T>& r, const T& v) {
            lapis::CoordTransform transformToExt = lapis::CoordTransform(crs, r.projection());

            PlotList newpl;

            for (const auto& plot : plots) {
                auto pt = lapis::Point(plot.second.second, plot.second.first, crs); //make Point from lon, lat.
                pt.project(transformToExt);
                auto atPlot = r.extract(pt.getX(), pt.getY());
                if (atPlot.has_value() && atPlot.value() == v) {
                    newpl.emplace(plot.first, plot.second);
                }
            }
            plots = newpl;
            return plots.size();
        }

        void makePlotTreeMap(const std::vector<std::string> colNames);
        FIATreeList collapsePlotTreeMap();
        FIATreeList collapseByPlotNames(const std::vector<std::string>& names);
        void calcKNNTree();

    private:
        std::string fiaFolder;
        lapis::CoordRef crs{ "4326" };

        const std::regex plotCsvRegex{ ".*PLOT\\.csv", std::regex_constants::icase };
        std::regex treeCsvRegex{ ".*TREE\\.csv",std::regex_constants::icase };

        const std::regex xRegex{ "\"?LAT\"?" };
        const std::regex yRegex{ "\"?LON\"?" };
        const std::regex nameRegex{ "\"?CN\"?" };

        void addPlotsFromFile(const std::string& fiaPlotFile);
    };

    
    using AllometryRaster = lapis::Raster < std::shared_ptr<Model>>;

    template<class T>
    AllometryRaster calculateAllometryOverProjectArea(lapis::Raster<T> r) {
        throw std::exception("Not implemented");
    }
} // namespace rxtools::allometry

#endif // !rxtools_allometry_h