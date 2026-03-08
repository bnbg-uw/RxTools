#pragma once

#ifndef rxtools_allometry_h
#define rxtools_allometry_h

#include "rxtools_pch.hpp"
#include "utilities.hpp"
#include "Raster.hpp"

namespace rxtools::allometry {
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;

    //string is the name of the plot, pair is x,y coords of plot.
    //TODO: Consider if PlotList can be a std::vector<std::pair<bg::point<...>, std::string>
    //      this way it could be directly inserted into boost rtree rather than constructing that which is the current method.
    //      this would optimize calcKNNTree at the cost of potentially having duplicates in the plotlist.
    //      also plotTreeMap could use the pair as the key? fewer objects at the cost of complexity.
    using PlotList = std::unordered_map<std::string, lapis::CoordXY>;

    //abstract base class for various types allometric models we can run.
    class Model {
    public:
        lapis::Unit inputUnit = lapis::linearUnitPresets::unknownLinear; //units the model expects as input
        lapis::Unit outputUnit = lapis::linearUnitPresets::unknownLinear; //units the model expects as output

        virtual double predict(double x, const lapis::Unit& thisUnit, const lapis::Unit& returnUnit) const = 0;
        std::vector<double> predict(const std::vector<double>& x, const lapis::Unit& thisUnit, const lapis::Unit& returnUnit);

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

    //Allows for more than one explanatory variable
    struct FIATreeList {
        std::vector<double> height;

        std::vector<std::string> names;
        std::vector<std::vector<double>> otherfields;

        void writeCsv(std::filesystem::path path) const {
            std::ofstream out;
            out.open(path);

            out << "height";
            for (int i = 0; i < names.size(); ++i) {
                out << "," << names[i];
            }
            out << "\n" << std::setprecision(std::numeric_limits<double>::max_digits10);
            for (int i = 0; i < height.size(); i++) {
                out << height.at(i);
                for (int j = 0; j < names.size(); ++j) {
                    out << "," << otherfields.at(i).at(j);
                }
                out << "\n";
            }
            out.close();
        }

        std::vector<double> get(const std::string& name) const {
            auto it = std::find(names.begin(), names.end(), name);
            if (it == names.end()) {
                throw(std::out_of_range("name is not in names of treeList"));
            }
            else {
                ptrdiff_t idx = std::distance(names.begin(), it);
                std::vector<double> out;
                for (size_t i = 0; i < otherfields.size(); ++i) {
                    out.push_back(otherfields.at(i).at(idx));
                }
                return out;
            }
        }
    };

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
            lapis::CoordTransform transformToExt = lapis::CoordTransformFactory::getTransform(crs, r.crs());

            PlotList newpl;

            for (const auto& plot : plots) {
                auto pt = lapis::Point(plot.second.x, plot.second.y, crs); //make Point from lon, lat.
                pt.projectInPlace(transformToExt);
                auto atPlot = r.extract(pt.x(), pt.y(), lapis::ExtractMethod::near);
                if (atPlot.has_value() && atPlot.value() == v) {
                    newpl.emplace(plot.first, plot.second);
                }
            }
            plots = newpl;
            std::erase_if(plotTreeMap, [&](const auto& keyval) {
                return plots.find(keyval.first) == plots.end();
                });
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

        const std::regex xRegex{ "\"?LON\"?" };
        const std::regex yRegex{ "\"?LAT\"?" };
        const std::regex nameRegex{ "\"?CN\"?" };

        void addPlotsFromFile(const std::string& fiaPlotFile);
    };

    
    using AllometryRaster = lapis::Raster < std::shared_ptr<Model>>;

    template<class T, class MODEL>
    AllometryRaster calculateAllometryOverProjectArea(lapis::Raster<T> r, FIAReader fia, std::string responseName, lapis::Unit responseUnit) {
        int id = 0;
        std::vector<std::pair<std::unordered_set<FIAReader::PlotPointName>, int>> knnToId;
        std::unordered_map<int, std::unordered_set<FIAReader::PlotPointName>> idToKnn;
        std::unordered_map<int, std::vector<lapis::cell_t>> idToCells;

        AllometryRaster out(r);

        std::cout << "Calculating KNN for each cell..." << std::endl;
        for (lapis::cell_t c = 0; c < r.ncell(); ++c) {
            if (r[c].has_value()) {
                FIAReader::PlotPoint query_point(r.xFromCell(c), r.yFromCell());

                //Need a better way to pick k.;
                std::unordered_set<FIAReader::PlotPointName> results;
                fia.tree.query(bgi::nearest(query_point, 15), std::back_inserter(results));
                
                int foundId = -1;
                for (const auto& p : knnToId) {
                    if (p.first == results) {
                        foundId = p.second;
                    }
                }
                if (foundId > -1) {
                    idToCells.at(foundId).push_back(c);
                }
                else {
                    knnToId.push_back(std::make_pair(results, id));
                    idToKnn.at(id) = results;
                    idToCells.at(id).push_back(c);
                    id++;
                }
            }
        }

        std::cout << "Calculating allometry for each unique KNN..." << std::endl;
        for (const auto& p : idToKnn) {
            FIATreeList treeList;
            for(const auto& ppv : p.second) {
                auto thisTrees = fia.plotTreeMap.at(ppv.second);
                if (!treeList.names.size()) {
                    treeList.names = thisTrees.names;
                }
                auto newSize = treeList.height.size() + thisTrees.height.size();
                treeList.height.reserve(newSize);
                treeList.otherfields.reserve(newSize);
                treeList.height.insert(treeList.height.end(), thisTrees.height.begin(), thisTrees.height.end());
                treeList.otherfields.insert(treeList.otherfields.end(), thisTrees.otherfields.begin(), thisTrees.otherfields.end());
            }
            std::shared_ptr<MODEL> mPtr = std::make_shared<MODEL>(treeList, responseName, responseUnit);
            for (const auto& cell : idToCells.at(p.first)) {
                out[cell].value() = mPtr;
                out[cell].has_value() = true;
            }
        }
        return(out);
    }
} // namespace rxtools::allometry

#endif // !rxtools_allometry_h