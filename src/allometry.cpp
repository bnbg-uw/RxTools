#include "allometry.hpp"

namespace rxtools::allometry {
    std::vector<double> Model::predict(const std::vector<double>& x, const lapis::LinearUnit& thisUnit, const lapis::LinearUnit& returnUnit) {
        std::vector<double> out;
        for (std::size_t i = 0; i < x.size(); i++) {
            out.push_back(predict(x[i], thisUnit, returnUnit));
        }
        return out;
    }

    FIAReader::FIAReader(const std::string& fiaFolder) : fiaFolder(fiaFolder) {
        for (auto fn : std::filesystem::directory_iterator(fiaFolder)) {
            if (std::regex_match(fn.path().string(), plotCsvRegex)) {
                addPlotsFromFile(fn.path().string());
            }
        }
    }

    FIAReader::FIAReader(const std::string& fiaFolder, const lapis::Extent& e) : FIAReader(fiaFolder) {
        limitByExtent(e);
    }

    const lapis::CoordRef FIAReader::projection() {
        return crs;
    }

    void FIAReader::projection(const lapis::CoordRef newCrs) {
        crs = newCrs;
    }

    void FIAReader::project(const lapis::CoordRef& newCrs) {
        const lapis::CoordTransform& t = lapis::CoordTransformFactory::getTransform(crs, newCrs);
        for (auto& plot : plots) {
            plot.second = t.transformSingleXY(plot.second.x, plot.second.y);
        }
    }

    const std::size_t FIAReader::limitByExtent(const lapis::Extent& e) {
        const lapis::CoordTransform& transformToExt = lapis::CoordTransformFactory::getTransform(crs,e.crs());

        PlotList newpl;
        for (const auto& plot : plots) {
            auto check = transformToExt.transformSingleXY(plot.second.x, plot.second.y);
            if (e.contains(check.x, check.y)) {
                newpl.emplace(plot);
            }
        }
        plots = newpl;
        std::erase_if(plotTreeMap, [&](const auto& keyval) {
            return plots.find(keyval.first) == plots.end();
            });
        return plots.size();
    }

    void FIAReader::makePlotTreeMap(const std::vector<std::string> colNames) {
        auto before = std::chrono::high_resolution_clock::now();

        std::string plotID{ "PLT_CN" };
        std::string ht{ "ACTUALHT" };
        std::string htcd{ "HTCD" };

        for (auto fn : std::filesystem::directory_iterator(fiaFolder)) {
            if (std::regex_match(fn.path().string(), treeCsvRegex)) {
                csv::CSVReader csv(fn.path().string());
                auto h = csv.get_col_names();
                if (
                    std::find(h.begin(), h.end(), plotID) == h.end() || 
                    std::find(h.begin(), h.end(), ht) == h.end() ||
                    std::find(h.begin(), h.end(), htcd) == h.end()
                    ) {
                    throw std::runtime_error(plotID + "/" + ht + "/" + htcd + " not found in " + fn.path().string() + ".");
                }
                for (auto& name : colNames) {
                    if (std::find(h.begin(), h.end(), name) == h.end()) {
                        throw std::runtime_error(name + " not found in " + fn.path().string() + ".");
                    }
                }
                
                for(auto& row : csv) {
                    try {
                        std::string plot = row[plotID].get<>();
                        if (plots.count(plot) == 0) { //not one of the plots we're using
                            continue;
                        }
                        //codes 1 and 2 correspond to measuring height; codes 3 and 4 correspond to estimating height
                        if (row[htcd].get<int>() >= 3) {
                            continue;
                        }

                        //running all type conversions first so that if any fail, the code bails before changing the vectors
                        double thisHeight = row[ht].get<double>();
                        std::vector<double> theseColumns;
                        for (int i = 0; i < colNames.size(); ++i) {
                            theseColumns.push_back(std::move(row[colNames[i]].get<double>()));
                        }
                        if (plotTreeMap.count(plot) == 0) {
                            FIATreeList add;
                            add.names = colNames;
                            add.height.push_back(thisHeight);
                            add.otherfields.push_back(theseColumns);
                            plotTreeMap.emplace(plot, add);
                        }
                        else {
                            plotTreeMap.at(plot).height.push_back(thisHeight);
                            plotTreeMap.at(plot).otherfields.push_back(theseColumns);
                        }
                    }
                    catch (...) {} //some of the lines have missing entries; there's no harm in skipping them
                }
            }
        }
        auto after = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(after - before);
        std::cout << "Made plot tree map in " << duration.count() << " seconds.\n";
    }

    FIATreeList FIAReader::collapsePlotTreeMap() {
        FIATreeList out;
        for (auto& x : plotTreeMap) {
            if (!out.names.size()) {
                out.names = x.second.names;
            }
            auto newSize = out.height.size() + x.second.height.size();
            out.height.reserve(newSize);
            out.otherfields.reserve(newSize);
            out.height.insert(out.height.end(), x.second.height.begin(), x.second.height.end());
            out.otherfields.insert(out.otherfields.end(), x.second.otherfields.begin(), x.second.otherfields.end());
        }
        return out;
    }

    FIATreeList FIAReader::collapseByPlotNames(const std::vector<std::string>& names) {
        FIATreeList out;
        for (auto& x : plotTreeMap) {
            if (std::find(names.begin(), names.end(), x.first) != names.end()) {
                if (!out.names.size()) {
                    out.names = x.second.names;
                }
                auto newSize = out.height.size() + x.second.height.size();
                out.height.reserve(newSize);
                out.otherfields.reserve(newSize);
                out.height.insert(out.height.end(), x.second.height.begin(), x.second.height.end());
                out.otherfields.insert(out.otherfields.end(), x.second.otherfields.begin(), x.second.otherfields.end());
            }
        }
        return out;
    }

    void FIAReader::calcKNNTree() {
        std::vector<PlotPointName> plotPairs;
        for (auto& plot : plots) {
            plotPairs.push_back(std::make_pair(PlotPoint(plot.second.x, plot.second.y), plot.first));
        }
        tree = bgi::rtree<PlotPointName, bgi::quadratic<16>>(plotPairs);
    }

    void FIAReader::addPlotsFromFile(const std::string& fiaPlotFile) {
        csv::CSVReader csv(fiaPlotFile);
        int nskip = 0;
        int ntot = 0;
        for(auto& row : csv) {
            ntot++;
            try {
                plots.emplace(row["CN"].get<>(), lapis::CoordXY((lapis::coord_t)row["LON"].get<double>(), (lapis::coord_t)row["LAT"].get<double>()));
            }
            catch (std::runtime_error e) {
                nskip++;
                continue;
            }
        }
        std::cout << "Loaded " << plots.size() << " plots from " << fiaPlotFile << ". Skipped " << nskip << " of " << ntot << " total.\n";
    }
} //namespace rxtools::allometry