#pragma once
#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "boost/filesystem/path.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/operations.hpp"
#include "Eigen/Dense"
#include "eigen/src/Eigenvalues/EigenSolver.h"
#include<boost/geometry.hpp>
#include "statistics/Statistics.hpp"
#include<regex>
#include<filesystem>


namespace licosim {
    namespace bf = boost::filesystem;
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;

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

    inline bf::path getExeLocation(int& argc, char* argv[]) {
        bf::path fullPath(bf::initial_path<bf::path>());
        if (argc > 0) {
            fullPath = bf::system_complete(bf::path(argv[0]));
            return fullPath.parent_path();
        }
        else {
            throw std::invalid_argument("You didn't call this from the cmd line? idk you fucked up somehow.");
        }
    }

    struct Pca {
    public:
        Eigen::VectorXd means;
        Eigen::VectorXd stdDevs;
        Eigen::MatrixXd eigenVectors;
        Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver;
        Eigen::MatrixXd data;
        Eigen::VectorXd importance;

        Pca() {};

        Pca(Eigen::MatrixXd inData) {
            //Calculate column means
            means = Eigen::VectorXd::Zero(inData.cols());
            for (int c = 0; c < inData.cols(); ++c) {
                for (int r = 0; r < inData.rows(); ++r)
                    means(c) += inData(r, c);
                means(c) /= inData.rows();
            }

            //Center the data
            data = Eigen::MatrixXd(inData.rows(), inData.cols());
            for (int r = 0; r < inData.rows(); ++r)
                for (int c = 0; c < inData.cols(); ++c)
                    data(r, c) = inData(r, c) - means(c);

            //Getting SD to scale the data
            stdDevs = Eigen::VectorXd::Zero(inData.cols());
            for (int r = 0; r < data.rows(); ++r)
                for (int c = 0; c < data.cols(); ++c)
                    stdDevs(c) += data(r, c) * data(r, c);
            int n = std::max(1, (int)inData.rows() - 1);
            for (int c = 0; c < stdDevs.size(); ++c)
                stdDevs(c) = std::sqrt(stdDevs(c) / n);

            //Scale the data
            for (int r = 0; r < data.rows(); ++r)
                for (int c = 0; c < data.cols(); ++c)
                    data(r, c) /= stdDevs(c);

            //New colmeans for covariance matrix
            Eigen::VectorXd scaledMeans = Eigen::VectorXd::Zero(data.cols());
            for (int c = 0; c < data.cols(); ++c) {
                for (int r = 0; r < data.rows(); ++r)
                    scaledMeans(c) += data(r, c);
                scaledMeans(c) /= data.rows();
            }

            //Making the covariance matrix.
            Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(data.cols(), data.cols());
            for (int c1 = 0; c1 < data.cols(); ++c1)
                for (int c2 = 0; c2 < data.cols(); ++c2)
                    for (int r = 0; r < data.rows(); ++r)
                        cov(c1, c2) += (data(r, c1) - scaledMeans(c1)) * (data(r, c2) - scaledMeans(c2)) / (data.rows() - 1);

            //Compute the eigenvalues and eigenvectors of the cov matrix
            eigenSolver.compute(cov);
            eigenVectors = eigenSolver.pseudoEigenvectors();
            for (int i = 0; i < 3; ++i) {
                eigenVectors.rowwise().reverse().transposeInPlace();
            }
            eigenVectors = eigenVectors.rowwise().reverse().eval();
            eigenVectors.transposeInPlace();
            data = predict(inData);

            //Calculate variable importance for weighting distances
            //this method is lifted from r's summary.prcomp:
            //https://github.com/wch/r-source/blob/5201e090c1623bb9fa04dbdeed0ab99903643f02/src/library/stats/R/prcomp.R
            //reusing colmeans vector from earlier.
            scaledMeans.setZero();
            for (int c = 0; c < data.cols(); ++c) {
                for (int r = 0; r < data.rows(); ++r)
                    scaledMeans(c) += data(r, c);
                scaledMeans(c) /= data.rows();
            }

            Eigen::VectorXd thisSd = Eigen::VectorXd::Zero(data.cols());
            double sum = 0;
            for (int c = 0; c < data.cols(); ++c) {
                for (int r = 0; r < data.rows(); ++r)
                    thisSd(c) += (data(r, c) - scaledMeans(c)) * (data(r, c) - scaledMeans(c));
                thisSd(c) /= data.rows() - 1;
                thisSd(c) = sqrt(thisSd(c));
                thisSd(c) = thisSd(c) * thisSd(c);
                sum += thisSd(c);
            }

            importance = Eigen::VectorXd(data.cols());
            for (int c = 0; c < importance.size(); ++c)
                importance(c) = thisSd(c) / sum;
        }
        

        template<class T>
        T predict(T input) {
            if (input.cols() != eigenVectors.cols())
                throw std::invalid_argument("Need data to be predicted to have same num of variables as PCA");
            for (int r = 0; r < input.rows(); ++r) {
                for (int c = 0; c < input.cols(); ++c) {
                    input(r, c) = (input(r, c) - means(c)) / stdDevs(c);
                }
            }
            return (input * eigenVectors).eval();
        }
    };

    template<size_t NDim>
    class PCA_KDTree { //A KD tree that fills itself up with the data a PCA was defined on. The distance metric used is weighted euclidean, weighted by the importance of each axis.
    public:
        typedef bg::model::point<double, NDim, bg::cs::cartesian> point;
        typedef std::pair<point, int> value;

        bgi::rtree<value, bgi::quadratic<16>> rt;
        std::vector<double> sqrt_weights;
        Pca pca;

        PCA_KDTree() {}
        PCA_KDTree(Pca p) {
            for (int i = 0; i < p.importance.size(); ++i) {
                sqrt_weights.push_back(std::sqrt(p.importance[i]));
            }
            for (int i = 0; i < p.data.rows(); ++i) {
                std::vector<double> thisValues(NDim);
                for (int j = 0; j < NDim; ++j) {
                    thisValues[j] = sqrt_weights[j] * p.data(i, j);
                }
                point thispoint;
                fillPoint<NDim-1>(thispoint, thisValues);
                rt.insert(std::make_pair(thispoint, i));
            }
            pca = p;
        }
        PCA_KDTree(Pca p, std::vector<bool> use) {
            for (int i = 0; i < p.importance.size(); ++i) {
                sqrt_weights.push_back(std::sqrt(p.importance[i]));
            }
            for (int i = 0; i < p.data.rows(); ++i) {
                if (use[i]) {
                    std::vector<double> thisValues(NDim);
                    for (int j = 0; j < NDim; ++j) {
                        thisValues[j] = sqrt_weights[j] * p.data(i, j);
                    }
                    point thispoint;
                    fillPoint<NDim - 1>(thispoint, thisValues);
                    rt.insert(std::make_pair(thispoint, i));
                }
            }
            pca = p;
        }

        std::vector<value> query(Eigen::Matrix<double, 1, NDim> p, int k) {
            auto thisValues = transform(p);
            point thispoint;
            fillPoint<NDim-1>(thispoint, thisValues);
            std::vector<value> out;
            rt.query(bgi::nearest(thispoint, k), std::back_inserter(out));
            return out;
        }

        std::vector<int> detailedQuery(Eigen::Matrix<double, 1, NDim> p, int k, bool sort = false, double maxdist = 0) {
            auto q = query(p, k);
            std::vector<int> out;
            std::vector<std::pair<int,double>> zipped;
            double sqdist = maxdist * maxdist;
            if (!sort && maxdist <= 0) {
                for (int i = 0; i < q.size(); ++i) {
                    out.push_back(q[i].second);
                }
                return out;
            }
            for (int i = 0; i < q.size(); ++i) {
                auto transformed = transform(p);
                double thisdist = calcdist<NDim - 1>(q[i].first, transformed, 0);
                if (maxdist > 0 && thisdist > sqdist) {
                    continue;
                }
                zipped.push_back(std::make_pair(q[i].second, thisdist));
            }
            if (sort) {
                std::sort(std::begin(zipped), std::end(zipped),
                    [&](const auto& a, const auto& b) {
                        return a.second < b.second;
                    });
            }
            for (int i = 0; i < zipped.size(); ++i) {
                out.push_back(zipped[i].first);
            }
            return out;
        }

    private:
       template<size_t N>
       void fillPoint(point& p, const std::vector<double>& values) {
            p.set<N>(values[N]);
            fillPoint<N - 1>(p, values);
       }
       template<>
       void fillPoint<0>(point& p, const std::vector<double>& values) {
           p.set<0>(values[0]);
       }

       std::vector<double> transform(Eigen::Matrix<double, 1, NDim> p) {
           auto transformed = pca.predict(p);
           std::vector<double> thisValues(NDim);
           for (int i = 0; i < NDim; ++i) {
               thisValues[i] = sqrt_weights[i] * transformed(0, i);
           }
           return thisValues;
       }

       template<size_t N>
       double calcdist(const point& p, const std::vector<double>& values, double currentdist) {
           double thisdist = p.get<N>() - values[N];
           thisdist *= thisdist;
           currentdist += thisdist;
           return calcdist<N - 1>(p, values, currentdist);
       }
       template<>
       double calcdist<0>(const point& p, const std::vector<double>& values, double currentdist) {
           double thisdist = p.get<0>() - values[0];
           thisdist *= thisdist;
           return currentdist + thisdist;
       }
    };

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
}