#pragma once

#include "rxtools_pch.hpp"

namespace rxtools::utilities {
    inline int randomIdxFromWeights(std::vector<double>& weights, std::default_random_engine& dre) {
        double sumOfWeight = 0;
        for (int i = 0; i < weights.size(); i++) {
            sumOfWeight += weights[i];
        }
        std::uniform_real_distribution<double> urd(0., sumOfWeight);
        double v = urd(dre);

        for (int i = 0; i < weights.size(); i++) {
            if (v < weights[i])
                return i;
            v -= weights[i];    
        }
        return weights.size() - 1;
    }

    inline std::vector<std::string> readCSVLine(std::istream& file) {
        char split = ',';
        std::string line = std::string();
        std::getline(file, line);
        std::istringstream i{ line };
        std::vector<std::string> out{};
        while (!i.eof()) {
            std::string token = std::string();
            std::getline(i, token, split);
            out.push_back(token);
        }
        return out;
    }

    struct fastFuelsAllometry {
    public:
        bool init = false;

        fastFuelsAllometry() {};

        fastFuelsAllometry(std::string fiaPath) {
            init = true;
            usePlots = stats::FIALeastSquares::getPlotList(fiaPath);
            cr = stats::FIALeastSquares(usePlots, fiaPath, "CR");

            std::regex treecsvregex{ ".*TREE\\.csv", std::regex_constants::icase };
            std::vector<int> spcd;

            for (auto fn : std::filesystem::directory_iterator(fiaPath)) {
                std::regex plotregex{ "\"?PLT_CN\"?" }; int plotidx = -1;
                std::regex spcdregex{ "\"?SPCD\"?" }; int spcdidx = -1;
                
                if (std::regex_match(fn.path().string(), treecsvregex)) {
                    std::ifstream ifs{ fn.path().string() };
                    auto colnames = readCSVLine(ifs);
                    for (int i = 0; i < colnames.size(); ++i) {

                        if (std::regex_match(colnames[i], plotregex)) {
                            plotidx = i;
                        }
                        if (std::regex_match(colnames[i], spcdregex)) {
                            spcdidx = i;
                        }
                    }

                    if (plotidx < 0 || spcdidx < 0) {
                        throw std::runtime_error(fn.path().string() + " is not formatted correctly.");
                    }

                    while (!ifs.eof()) {
                        auto line = readCSVLine(ifs);
                        if (line.size() <= 1) {
                            continue;
                        }
                        try {
                            std::string plot = line[plotidx];
                            if (usePlots.count(plot) == 0) { //not one of the plots we're using
                                continue;
                            }
                            int thisspcd = std::stoi(line[spcdidx]);
                            spcd.push_back(thisspcd);
                        }
                        catch (...) {} //some of the lines have missing entries; there's no harm in skipping them
                    }
                }
            }

            std::unordered_map<int, int> spcdTable;
            for (auto v : spcd) {
                spcdTable.emplace(v, 0);
                ++spcdTable[v];
            }
            for (auto& v : spcdTable) {
                v.second = (int)std::round((double)v.second / spcd.size() * 100);
                for (size_t i = 0; i < v.second; ++i) {
                    spcdWeights.push_back(v.first);
                }
            }
            std::cout << spcdWeights.size();

        }

        double predictCbh(const double& height) {
            auto crown = cr.predict(height);
            return crown * height / 100;
        }

        int assignSpecies(const double& x, const double& y) {
            return spcdWeights[spcdHash(x * y) % spcdWeights.size()];
        }

    private:
        stats::PlotList usePlots;
        std::vector<int> spcdWeights;
        std::hash<double> spcdHash;
        stats::FIALeastSquares dia;
        stats::FIALeastSquares cr;

        std::size_t limitByExtent(const spatial::Extent& e) {
            using namespace spatial;
            auto y = CPLGetConfigOption("PROJ_DATA", nullptr);
            std::cout << y << "\n";
            y = proj_context_get_database_path(nullptr);
            std::cout << y << "\n";
            CoordRef fiacrs = CoordRef("4326");
            CoordTransform transformToExt{ fiacrs,e.projection() };

            stats::PlotList newpl;

            for (const auto& plot : usePlots) {
                auto pt = SpPoint(plot.second.first, plot.second.second, fiacrs); //make SpPoint from lon, lat.
                pt.project(transformToExt);
                if (e.contains(pt.getX(), pt.getY())) {
                    newpl.emplace(plot.first, plot.second);
                }
            }
            usePlots = newpl;
            return usePlots.size();
        }
    };
} // namepsace rxtools::utilities