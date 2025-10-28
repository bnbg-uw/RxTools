#pragma once

#ifndef licosim_treatment_h
#define licosim_treatment_h

#include <stdexcept>
#include <random>
#include <unordered_set>
#include "lico/LICO.hpp"
#include "licosim/rxunit.hpp"
#include "licosim/structuresummary.hpp"
#include<LICO/GraphLico.hpp>
#include<limits>
#include<random>

// Need to go through and see where I can pass by reference for memory (and some speed) optimization
namespace licosim {


enum class treatmentDecision { undecided, retain, cut };
enum class treatmentResult { success, diameterFailure, cuttingFailure };

using td = treatmentDecision;
using index_t = lico::index_t;

struct NodePriority {
    lico::index_t priority;
    lico::index_t index;

    NodePriority(lico::index_t p, lico::index_t i) : priority(p), index(i) {}
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

    lico::TaoList doTreatment(licosim::RxUnit rx, double dbhMin, double dbhMax);

    //a version of doTreatment with some optimizations
    //Returns a tuple: retained trees, cut trees, treatment success
    //Treatment is "successful" if treatment is achieved without having to add any extra trees back in or being unsuccessful due to diameter limit.
    //"diameterFailure" is if the treatment was unable to be undertaken because trees > diameter limit had more than the target BA.
    //"cuttingFailure" means that too many trees were cut and some had to be added back in due to residual BA being below target BA.
    std::tuple<lico::TaoList, lico::TaoList, treatmentResult> doTreatmentGraph(RxUnit rx, double dbhMin, double dbhMax, lico::unit_t maxCrown, bool intermediates = false, std::string intermediatespath = "E:/treatedcsvs/");

private:

    // Empirical estimation of weighting factor to align targeted BA/TPA/QMD.
    // dbh = dbh of all trees.
    // treeStatus = retain undecided cut coding for dbh
    // targets = treatment targets
    // area = lmu area
    // wLim = clamps for the weighting factors. Defaults from Sean Jeronimo. wLim.second must be greater than wLim.first
    // nw = number of weights to test in between wLim
    double calcW(std::vector<double> dbh, std::vector<treatmentDecision> treeStatus, StructureSummary targets,
                 double area, std::pair<double, double> wLim = std::pair<double, double>(-20., 20.), int nw = 200);

    // Calculate how much is needed
    std::vector<double> calcNeed(lico::TaoList tl, std::vector<double> dbh, std::vector<treatmentDecision> treeStatus,
                                 double area, double baTargs, std::vector<double> csdTargs);

    double calcBA(std::vector<double> dbh, std::vector<treatmentDecision> treeStatus, double area, bool countUndecided = false);

    double calcTPH(std::vector<treatmentDecision> treeStatus, double area, bool countUndecided = false);

    std::vector<double> calcClumpProps(lico::TaoList tl);

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
    void obligateIdx(lico::TaoList& tl, std::vector<treatmentDecision>& treeStatus, std::vector<index_t>& treeIdx, index_t& focalIdx, std::unordered_set<index_t>& keepSet);

};
namespace weightedRandom {
    struct weightedNode {
        double leftWeight, rightWeight, thisWeight;
    };
    struct Weights {
        std::vector<lico::index_t> sample;
        std::vector<double> weights;
    };

    //Functions for doing the random weighted sorting of the tree list
    inline double totalWeight(std::vector<weightedNode>& nodes, const std::vector<double>& weights, lico::index_t i = 0) {
        if (i >= nodes.size())
            return 0;
        if (weights[i] <= 0) {
            throw(std::runtime_error("Weights must be positive"));
        }
        nodes[i].thisWeight = weights[i];
        lico::index_t left = 2 * i + 1;
        lico::index_t right = 2 * i + 2;
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

    inline lico::index_t weightedSampleWithoutReplace(std::vector<weightedNode>& nodes, const std::vector<double>& weights, uniform_random& dist, std::default_random_engine& dre, lico::index_t i = 0) {
        lico::index_t left = 2 * i + 1;
        lico::index_t right = 2 * i + 2;

        if (nodes[i].leftWeight == 0 && nodes[i].rightWeight == 0) {
            nodes[i].thisWeight = 0;
            return i;
        }
        if (nodes[i].thisWeight == 0 && nodes[i].rightWeight == 0) {
            lico::index_t chosen = weightedSampleWithoutReplace(nodes, weights, dist, dre, left);
            nodes[i].leftWeight = nodes[left].leftWeight + nodes[left].thisWeight + nodes[left].rightWeight;
            return chosen;
        }
        if (nodes[i].thisWeight == 0 && nodes[i].leftWeight == 0) {
            lico::index_t chosen = weightedSampleWithoutReplace(nodes, weights, dist, dre, right);
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
            lico::index_t chosen = weightedSampleWithoutReplace(nodes, weights, dist, dre, left);
            nodes[i].leftWeight = nodes[left].leftWeight + nodes[left].thisWeight + nodes[left].rightWeight;
            return chosen;
        }
        else {
            lico::index_t chosen = weightedSampleWithoutReplace(nodes, weights, dist, dre, right);
            nodes[i].rightWeight = nodes[right].leftWeight + nodes[right].thisWeight + nodes[right].rightWeight;
            return chosen;
        }
    }
    Weights getWeightedOrder(lico::GraphLico& g, StructureSummary targets, double area, std::default_random_engine& dre, std::pair<double, double> wlim = std::pair<double, double>(-20., 20.), int nw = 200);
    Weights getWeightedOrderByRatios(lico::GraphLico& g, StructureSummary targets, double area, std::default_random_engine& dre);
}
} //namespace licosim

#endif  // !licosim_treatment_h