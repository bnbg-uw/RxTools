#pragma once

#ifndef rxtools_allometry_h
#define rxtools_allometry_h

#include "rxtools_pch.hpp"
#include "utilities.hpp"
#include "LapisGis/src/Raster.hpp"

namespace rxtools::allometry {
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;

    //string is the name of the plot, pair is x,y coords of plot.
    //TODO: Consider if PlotList can be a std::vector<std::pair<bg::point<...>, std::string>
    //      this way it could be director inserted into boost rtree rather than constructing that which is the current method.
    //      this would optimize calcKNNTree at the cost of potentially having duplicates in the plotlist.
    //      also plotTreeMap could use the pair as the key? fewer objects at the cost of complexity.
    using PlotList = std::unordered_map<std::string, lapis::CoordXY>;
    using AllometryRaster = lapis::Raster < std::shared_ptr<Model>>;

    //Allows for more than one explanatory variable
    struct FIATreeList {
        std::vector<double> height;

        std::vector<std::string> names;
        std::vector<std::vector<double>> otherfields;
    };

    //abstract base class for various types allometric models we can run.
    //all models will expect input in imperial units and will return imperial units.
    //
    class Model {
    public:
        double convFactor = 0; //what to multipy the response by to make it opposite units output
        bool isModelMetric = false; //is the model expecting input and spitting out metric units or not

        virtual double predict(const double x, const bool metric = false) const = 0;
        std::vector<double> predict(const std::vector<double> x, const bool metric = false);

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
            ...
        }

        void makePlotTreeMap(const std::vector<std::string> colNames);

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

    template<class T>
    AllometryRaster calculateAllometryOverProjectArea(lapis::Raster<T> r) {
        intentional compiler error;
    }
} // namespace rxtools::allometry

#endif // !rxtools_allometry_h