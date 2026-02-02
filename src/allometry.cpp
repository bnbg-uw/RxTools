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

        return plots.size();
    }

    void FIAReader::makePlotTreeMap(const std::vector<std::string> colNames, int nThread) {
        auto before = std::chrono::high_resolution_clock::now();

        for (auto fn : std::filesystem::directory_iterator(fiaFolder)) {
            std::regex plotRegex{ "\"?PLT_CN\"?" }; int plotIdx = -1;
            std::regex htRegex{ "\"?ACTUALHT\"?" }; int htIdx = -1;
            std::regex htcdRegex{ "\"?HTCD\"?" }; int htcdIdx = -1;

            std::vector<std::regex> colNamesRegex; {};
            std::vector<int> colNamesIdx;
            for (int i = 0; i < colNames.size(); ++i) {
                colNamesRegex.push_back(std::regex("\"?" + colNames[i] + "\"?", std::regex_constants::icase));
                colNamesIdx.push_back(-1);
            }

            if (std::regex_match(fn.path().string(), treeCsvRegex)) {
                std::ifstream ifs{ fn.path().string() };
                auto colnames = utilities::readCSVLine(ifs);
                for (int i = 0; i < colnames.size(); ++i) {

                    if (std::regex_match(colnames[i], plotRegex)) {
                        plotIdx = i;
                    }
                    if (std::regex_match(colnames[i], htRegex)) {
                        htIdx = i;
                    }
                    if (std::regex_match(colnames[i], htcdRegex)) {
                        htcdIdx = i;
                    }
                    for (int j = 0; j < colNames.size(); ++j) {
                        if (std::regex_match(colnames[i], colNamesRegex[j])) {
                            colNamesIdx[j] = i;
                        }
                    }

                }

                if (plotIdx < 0 || htIdx < 0 || htcdIdx < 0) {
                    throw std::runtime_error(fn.path().string() + " is not formatted correctly.");
                }
                for (int i = 0; i < colNames.size(); ++i) {
                    if(colNamesIdx[i] < 0)
                        throw std::runtime_error(colNames[i] + " not found in " + fn.path().string() + ".");
                }

                while (!ifs.eof()) {
                    auto line = utilities::readCSVLine(ifs);
                    if (line.size() <= 1) {
                        continue;
                    }
                    try {
                        std::string plot = line[plotIdx];
                        if (plots.count(plot) == 0) { //not one of the plots we're using
                            continue;
                        }

                        //codes 1 and 2 correspond to measuring height; codes 3 and 4 correspond to estimating height
                        if (std::stoi(line[htcdIdx]) >= 3) {
                            continue;
                        }

                        //running all stods first so that if either fails, the code bails before changing the vectors
                        double thisHeight = std::stod(line[htIdx]);
                        std::vector<double> theseColumns;
                        for (int i = 0; i < colNames.size(); ++i) {
                            theseColumns.push_back(std::stod(line[colNamesIdx[i]]));
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
                        //f << plot << "," << thisheight << "," << thisresponse << "," << std::stod(line[spcdidx]) << "\n";
                    }
                    catch (...) {} //some of the lines have missing entries; there's no harm in skipping them
                }
            }
        }
        auto after = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(after - before);
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
        std::ifstream ifs{ fiaPlotFile };

        auto colnames = utilities::readCSVLine(ifs);
        int xidx = -1;
        int yidx = -1;
        int nameidx = -1;
        for (int i = 0; i < colnames.size(); ++i) {
            if (std::regex_match(colnames[i], xRegex)) {
                xidx = i;
            }
            else if (std::regex_match(colnames[i], yRegex)) {
                yidx = i;
            }
            else if (std::regex_match(colnames[i], nameRegex)) {
                nameidx = i;
            }
        }
        if (xidx == -1 || yidx == -1 || nameidx == -1) {
            throw std::runtime_error(fiaPlotFile + " is not formatted correctly.");
        }

        while (!ifs.eof()) {
            auto row = utilities::readCSVLine(ifs);
            if (row.size() <= 1) {
                continue;
            }
            try {
                lapis::coord_t x = std::stod(row[xidx]);
                lapis::coord_t y = std::stod(row[yidx]);
                std::string name = row[nameidx];
                plots.emplace(name, lapis::CoordXY(x, y));
            }
            catch (...) { //there will be ocassional broken lines
                continue;
            }
        }
    }

    void plotsFromFileThread(std::ifstream& ifs, std::mutex& mut) {

    }
} //namespace rxtools::allometry