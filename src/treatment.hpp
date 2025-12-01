#pragma once

#ifndef rxtools_treatment_h
#define rxtools_treatment_h

#include "rxtools_pch.hpp"
#include "lico/src/GraphLico.hpp"
#include "rxunit.hpp"
#include "structuresummary.hpp"

//go through and find where pass by reference can happen.
namespace rxtools {
    enum class treatmentDecision { undecided, retain, cut };
    enum class treatmentResult { success, diameterFailure, cuttingFailure };

    using td = treatmentDecision;

    struct NodePriority {
        size_t priority;
        size_t index;

        NodePriority(size_t p, size_t i) : priority(p), index(i) {}
    };

    struct NodePriorityComparator {
        bool operator()(const NodePriority& lhs, const NodePriority& rhs) {
            return lhs.priority > rhs.priority;
        }
    };

    #ifdef _DEBUG
    class test_random {
    public:
        std::queue<double> values;
        std::uniform_real_distribution<double> urd;
        test_random(double min, double max) : urd(min, max) {}
        double operator()(std::default_random_engine dre) {
            if (values.size() > 0) {
                double out = values.front();
                values.pop();
                return out;
            }
            return urd(dre);
        }
        void push(double v) {
            values.push(v);
        }
    };
    using uniform_random = test_random;
    #else
    using uniform_random = std::uniform_real_distribution<>;
    #endif

    class Treatment {
    public:
        std::default_random_engine dre;
     
        Treatment() {};
        Treatment(std::default_random_engine& re) : dre(re) {};

        //a version of doTreatment with some optimizations
        //Returns a tuple: retained trees, cut trees, treatment success
        //Treatment is "successful" if treatment is achieved without having to add any extra trees back in or being unsuccessful due to diameter limit.
        //"diameterFailure" is if the treatment was unable to be undertaken because trees > diameter limit had more than the target BA.
        //"cuttingFailure" means that too many trees were cut and some had to be added back in due to residual BA being below target BA.
        std::tuple<TaoListMP, TaoListMP, treatmentResult> doTreatment(RxUnit rx, double dbhMin, double dbhMax, lapis::coord_t maxCrown, bool intermediates = false, std::string intermediatespath = "E:/treatedcsvs/");

    private:
        // If a tree is within the limiting distance from a BB tree, then retaining that tree obligates you to retain the BB tree as well
        // This calculates 'how many trees will a given tree be "obligated" to bring along into the clump?'
        // It is a list of percolated adjacent backbone trees
        // 
        // tl is the Outer tao dataset
        // treeStatus is the total treeStatus vector for the trees in the current treatment.
        // treeIdx is a vector of idx that correspond the trees in the current clump to the their locations in the total treelist.
        // focalIdx is the idx of the focal tree in the total dataset
        // keepSet is an unordered set of the idx of the keep/obligate Trees in the total dataset.
        // What you get back is an undordered set of idx corresponding to keep trees in your clump that corresponds to idx in the total taolist.
        void obligateIdx(TaoListMP& tl, std::vector<treatmentDecision>& treeStatus, std::vector<size_t>& treeIdx, size_t& focalIdx, std::unordered_set<size_t>& keepSet);

    };
    namespace weightedRandom {
        struct weightedNode {
            double leftWeight, rightWeight, thisWeight;
        };
        struct Weights {
            std::vector<size_t> sample;
            std::vector<double> weights;
        };

        //Functions for doing the random weighted sorting of the tree list
        inline double totalWeight(std::vector<weightedNode>& nodes, const std::vector<double>& weights, size_t i = 0) {
            if (i >= nodes.size())
                return 0;
            if (weights[i] <= 0) {
                throw(std::runtime_error("Weights must be positive"));
            }
            nodes[i].thisWeight = weights[i];
            size_t left = 2 * i + 1;
            size_t right = 2 * i + 2;
            nodes[i].leftWeight = totalWeight(nodes, weights, left);
            nodes[i].rightWeight = totalWeight(nodes, weights, right);
            if (left > nodes.size()) {
                nodes[i].leftWeight = 0;
            }
            if (right > nodes.size()) {
                nodes[i].rightWeight = 0;
            }
            return nodes[i].thisWeight + nodes[i].leftWeight + nodes[i].rightWeight;
        }

        inline size_t weightedSampleWithoutReplace(std::vector<weightedNode>& nodes, const std::vector<double>& weights, uniform_random& dist, std::default_random_engine& dre, size_t i = 0) {
            size_t left = 2 * i + 1;
            size_t right = 2 * i + 2;

            if (nodes[i].leftWeight == 0 && nodes[i].rightWeight == 0) {
                nodes[i].thisWeight = 0;
                return i;
            }
            if (nodes[i].thisWeight == 0 && nodes[i].rightWeight == 0) {
                size_t chosen = weightedSampleWithoutReplace(nodes, weights, dist, dre, left);
                nodes[i].leftWeight = nodes[left].leftWeight + nodes[left].thisWeight + nodes[left].rightWeight;
                return chosen;
            }
            if (nodes[i].thisWeight == 0 && nodes[i].leftWeight == 0) {
                size_t chosen = weightedSampleWithoutReplace(nodes, weights, dist, dre, right);
                nodes[i].rightWeight = nodes[right].leftWeight + nodes[right].thisWeight + nodes[right].rightWeight;
                return chosen;
            }

            double totalWeight = nodes[i].thisWeight + nodes[i].leftWeight + nodes[i].rightWeight;
            double r = totalWeight * dist(dre);
            if (r < nodes[i].thisWeight) {
                nodes[i].thisWeight = 0;
                return i;
            }
            else if (r < (nodes[i].thisWeight + nodes[i].leftWeight)) {
                size_t chosen = weightedSampleWithoutReplace(nodes, weights, dist, dre, left);
                nodes[i].leftWeight = nodes[left].leftWeight + nodes[left].thisWeight + nodes[left].rightWeight;
                return chosen;
            }
            else {
                size_t chosen = weightedSampleWithoutReplace(nodes, weights, dist, dre, right);
                nodes[i].rightWeight = nodes[right].leftWeight + nodes[right].thisWeight + nodes[right].rightWeight;
                return chosen;
            }
        }
        Weights getWeightedOrder(lapis::lico::GraphLico& g, rxtools::StructureSummary targets, double area, std::default_random_engine& dre, std::pair<double, double> wlim = std::pair<double, double>(-20., 20.), int nw = 200);
        Weights getWeightedOrderByRatios(lapis::lico::GraphLico& g, rxtools::StructureSummary targets, double area, std::default_random_engine& dre);
    }
} //namespace rxtools

#endif  // !rxtools_treatment_h