#pragma once

#include "rxtools_pch.hpp"

namespace rxtools::utilities {
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;

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

    //A KD tree that fills itself up with the data a PCA was defined on. The distance metric used is weighted euclidean, weighted by the importance of each axis.
    template<size_t NDim>
    class PCA_KDTree {
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
                fillPoint<NDim - 1>(thispoint, thisValues);
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
            fillPoint<NDim - 1>(thispoint, thisValues);
            std::vector<value> out;
            rt.query(bgi::nearest(thispoint, k), std::back_inserter(out));
            return out;
        }

        std::vector<int> detailedQuery(Eigen::Matrix<double, 1, NDim> p, int k, bool sort = false, double maxdist = 0) {
            auto q = query(p, k);
            std::vector<int> out;
            std::vector<std::pair<int, double>> zipped;
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
}