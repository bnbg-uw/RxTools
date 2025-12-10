#pragma once

#include "rxtools_pch.hpp"
#include "rasteralgos.hpp"

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

    struct Kpoint {
    public:
        double x, y;
        int cell;
        int clusterid;

        Kpoint(double x, double y, int c) : x(x), y(y), cell(c) {
            clusterid = -1;
        }
    };

    struct Kcluster {
    public:
        Kpoint centroid;
        int id;
        std::vector<Kpoint> points;

        Kcluster(int id, Kpoint centroid) : id(id), centroid(centroid) {
            points.push_back(centroid);
        };

        int size() {
            return points.size();
        }
    };


    class Kmeans {
    public:
        int k, iters;
        std::vector<Kcluster> clusters;

        Kmeans(int k, int iters) : k(k), iters(iters) {}

        void clear() {
            for (int i = 0; i < k; ++i) {
                clusters[i].points.clear();
            }
        }

        int getNearestClusterId(Kpoint p) {
            double minDist = std::numeric_limits<double>::max();
            int id = -1;

            for (int i = 0; i < k; ++i) {
                double dist = (clusters[i].centroid.x - p.x) * (clusters[i].centroid.x - p.x) + (clusters[i].centroid.y - p.y) * (clusters[i].centroid.y - p.y);

                if (dist < minDist) {
                    minDist = dist;
                    id = clusters[i].id;
                }
            }

            return id;
        }

        void run(std::vector<Kpoint>& allPoints) {
            std::random_device rd;
            std::default_random_engine rng(rd());
            std::shuffle(allPoints.begin(), allPoints.end(), rng);

            for (int i = 0; i < k; ++i) {
                allPoints[i].clusterid = i;
                clusters.push_back(Kcluster(i, allPoints[i]));
            }

            int iter = 1;
            while (true) {
                bool done = true;

                for (int i = 0; i < allPoints.size(); ++i) {
                    int currentid = allPoints[i].clusterid;
                    int nearestid = getNearestClusterId(allPoints[i]);

                    if (currentid != nearestid) {
                        allPoints[i].clusterid = nearestid;
                        done = false;
                    }
                }
                if (done)
                    break;

                clear();
                for (int i = 0; i < allPoints.size(); ++i) {
                    clusters[allPoints[i].clusterid].points.push_back(allPoints[i]);
                }

                for (int i = 0; i < k; ++i) {
                    double x = 0, y = 0;
                    for (int j = 0; j < clusters[i].size(); ++j) {
                        x += clusters[i].points[j].x;
                        y += clusters[i].points[j].y;
                    }
                    clusters[i].centroid.x = x / clusters[i].size();
                    clusters[i].centroid.y = y / clusters[i].size();
                }
                if (iter >= iters)
                    break;
                iter++;
            }
        }
    };
} // namepsace rxtools::utilities