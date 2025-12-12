#include "treatment.hpp"

namespace rxtools {
    using ns = lapis::lico::NodeStatus;

    void Treatment::obligateIdx(TaoListMP& tl, std::vector<treatmentDecision>& treeStatus, std::vector<size_t>& treeIdx, size_t& focalIdx, std::unordered_set<size_t>& keepSet) {
        std::vector<size_t> testIdx; //These are idx in the total dataset to test as focals.
        
        //Get all the trees that are not already known as keeps, but are obligate retzins and clump with the focal tree
        for (size_t i : treeIdx) {
            if (treeStatus[i] == td::retain && keepSet.find(i) == keepSet.end()) {
                double dsq = (tl.x(focalIdx) - tl.x(i)) * (tl.x(focalIdx) - tl.x(i)) + (tl.y(focalIdx) - tl.y(i)) * (tl.y(focalIdx) - tl.y(i));
                if (dsq <= (tl.radius(focalIdx) + tl.radius(i)) * (tl.radius(focalIdx) + tl.radius(i)))
                    testIdx.push_back(i);
            }
        }

        keepSet.insert(focalIdx); //If a tree is focal we know we have to keep it.
        //Now check all these new trees for their obligates too. (We will add each of them to keepSet when doing this)
        for (size_t idx : testIdx)
            obligateIdx(tl, treeStatus, treeIdx, idx, keepSet);  //This can be optimized I believe... look into std::union_set or something maybe?
        // Since we pass keepSet by reference, it will be updated with the new taos as we go.
    }

    std::tuple<TaoListMP, TaoListMP, treatmentResult> Treatment::doTreatment(RxUnit rx, double dbhMin, double dbhMax, lapis::coord_t maxCrown, bool intermediates, std::string intermediatespath) {
        auto binMins = rx.targetStructure.binMins;
        auto binMaxs = rx.targetStructure.binMaxs;

        lapis::Alignment a{ rx.unitMask.xmin(), rx.unitMask.ymin(),(lapis::rowcol_t)std::ceil((rx.unitMask.ymax() - rx.unitMask.ymin()) / maxCrown),
            (lapis::rowcol_t)std::ceil((rx.unitMask.xmax() - rx.unitMask.xmin()) / maxCrown),maxCrown,maxCrown };
        lapis::lico::GraphLico g{ a };
        
        //absolute ba rather than per area
        double targetba = rx.targetStructure.ba * rx.areaHa;

        //current structure
        double currentba = 0;
        size_t currentntao = 0;

        //vector to convert between taolist in rxunit 
        //and the potentially subset list used in the graph because of dbhmin
        std::vector<size_t> taocrosswalk;

        //populate the graph and set up the crosswalk.
        auto tnf = lapis::lico::TaoNodeFactory<lapis::VectorDataset<lapis::MultiPolygon>>(
            rx.taos.getters.predicate,
            rx.taos.getters.xy,
            rx.taos.getters.radius,
            rx.taos.getters.area,
            rx.taos.getters.dbh
        );
        for (size_t i = 0; i < rx.taos.size(); ++i) {
            if (rx.taos.dbh(i) < dbhMin) {
                continue;
            }

            taocrosswalk.push_back(i);
            if (rx.taos.dbh(i) > dbhMax) {
                g.addTAO(tnf(rx.taos(i)), ns::on);
                currentba += g.nodes[g.nodes.size() - 1].ba;
                currentntao++;
            }
            else {
                g.addTAO(tnf(rx.taos(i)));
            }
        }

        //Vector of desired clump bins based on a lookup table from desired mcs
        std::vector<double> csdTargs = rx.targetStructure.csd;
        //converting the csd targets to be measured in basal area instead of proportion of basal area
        for (int i = 0; i < csdTargs.size(); ++i) {
            csdTargs[i] *= targetba;
        }

        //populate csd bins with existing ba.
        //should this check for node status?
        std::vector<double> csdCurrent(csdTargs.size());
        for (size_t i = 0; i < g.nodes.size(); ++i) {
            int clumpSize = g.nodes[i].clumpSize();
            for (int j = 0; j < csdCurrent.size(); ++j) {
                if (clumpSize >= binMins[j] && (clumpSize <= binMaxs[j] || j == csdCurrent.size()-1)) {
                    csdCurrent[j] += g.nodes[i].ba;
                }
            }
        }

        weightedRandom::Weights w = weightedRandom::getWeightedOrderByRatios(g, rx.targetStructure, rx.areaHa, dre);
        auto& sample = w.sample;
        auto& weights = w.weights;
        size_t sofar = 0;
        std::vector<size_t> priority(g.nodes.size());
        for (size_t i = 0; i < sample.size(); ++i) {
            priority[sample[i]] = i;
        }

        int step = 0;
        std::ofstream targetstream;
        if (intermediates) {
            rx.taos.writeCsv(intermediatespath + "-1.csv");
            targetstream.open(intermediatespath + "targets.csv");
            for (int i = 0; i < csdTargs.size(); ++i) {
                targetstream << csdTargs[i] << ",";
            }
            targetstream << "\n";
            for (int i = 0; i < csdCurrent.size(); ++i) {
                targetstream << csdCurrent[i] << ",";
            }
            targetstream << "\n";
        }

        treatmentResult result = treatmentResult::success;
        if (currentba > targetba) result = treatmentResult::diameterFailure;
        //std::cout << "start of while loop: " << sofar << " " << sample.size() << " " << currentba << " " << targetba << "\n";
        while (sofar < sample.size() && currentba < targetba) {
            size_t seedTao = sample[sofar];
            lapis::lico::TaoNode& seedNode = g.nodes[seedTao];
            //std::cout << "node status: " << (int)seedNode.status << "\n";
            if (seedNode.status != ns::off) {
                sofar++;
                continue;
            }

            sofar++;
            //int approxMaxClump = seedNode.maxClumpSizeApprox();
            int exactMaxClump = seedNode.maxClumpSizeExact();
            int maxClump = exactMaxClump; //For now, using the slow but precise calculation. Will test how far off the approximation is and possible switch which one we use
            int minClump = seedNode.peekClumpSize();

            //std::cout << "maxclump " << exactMaxClump << "minclump " << minClump << "sofar: " << sofar << "\n";
            std::vector<double> binWeights(csdTargs.size());
            double bwSum = 0;
            for (int i = 0; i < binWeights.size(); ++i) {
                //std::cout << csdCurrent[i] << " " << csdTargs[i] << " " << binMins[i] << binMaxs[i] << "\n";
                if (csdCurrent[i] > csdTargs[i]) {
                    binWeights[i] = 0;
                }
                else if (binMins[i] > maxClump) {
                    binWeights[i] = 0;
                }
                else if (binMaxs[i] < minClump) {
                    binWeights[i] = 0;
                }
                else {
                    binWeights[i] = csdTargs[i] - csdCurrent[i];
                    bwSum += binWeights[i];
                }
            }
            if (bwSum == 0) { // this TAO belongs to a clump so small or so large that it can't contribute to the bins that need more trees
                seedNode.maxClumpSizeApprox()--;
                seedNode.cutTAO();
                continue;
            }

            //seedNode.turnOn();
            //currentba += seedNode.ba;
            size_t binIdx = utilities::randomIdxFromWeights(binWeights, dre);
            std::uniform_int_distribution<int> d(binMins[binIdx], binMaxs[binIdx]);
            int targN = d(dre);
            if (targN > maxClump) {
                targN = maxClump;
            }
            int doNotExceed = binMaxs[binIdx];

            std::unordered_set<size_t> inClump{};
            //inClump.insert(seedNode.index);

            std::priority_queue<NodePriority, std::vector<NodePriority>, NodePriorityComparator> adjQueue;
            adjQueue.push(NodePriority{ priority[seedNode.index],seedNode.index });
            std::unordered_set<size_t> considered;
            considered.insert(seedNode.index);

            while (seedNode.clumpSize() < targN) {
                //pop off of adjQueue until we get a TAO which doesn't take us over the bin threshold
                size_t addTAO;
                bool giveup = false;
                bool foundone = false;
                while (!foundone) {
                    if (adjQueue.size() == 0) {
                        giveup = true;
                        break;
                    }
                    addTAO = adjQueue.top().index;
                    adjQueue.pop();
                    if (g.nodes[addTAO].peekClumpSize() <= doNotExceed) {
                        foundone = true;
                    }
                }
                if (giveup) {
                    break;
                }
                //add the TAOs adjacent to addTAO to adjQueue, and pick up any backbones adjacent to this TAO for inClump
                for (size_t adj : g.nodes[addTAO].adjList) {
                    if (g.nodes[adj].status == ns::on) {
                        inClump.insert(adj);
                    }
                    if (considered.count(adj)) {
                        continue;
                    }
                    if (g.nodes[adj].status != ns::off) {
                        continue;
                    }
                    adjQueue.push(NodePriority{ priority[adj], adj });
                    considered.insert(adj);
                }

                //updating csdCurrent to account for the clumps changing in size
                std::unordered_set<size_t> adjClumps;
                for (size_t adj : g.nodes[addTAO].adjList) {
                    if (adjClumps.count(g.nodes[adj].findAncestor().index)) {
                        continue;
                    }
                    adjClumps.insert(g.nodes[adj].findAncestor().index);
                    int size = g.nodes[adj].clumpSize();
                    double ba = g.nodes[adj].getClumpBA();
                    for (int i = 0; i < csdCurrent.size(); ++i) {
                        if (size >= binMins[i] && (size <= binMaxs[i] || i == csdCurrent.size() - 1)) {
                            csdCurrent[i] -= ba;
                        }
                    }
                }

                //turn on addTAO and update current ba and csd
                int size = g.nodes[addTAO].turnOn();
                currentba += g.nodes[addTAO].ba;
                double ba = g.nodes[addTAO].getClumpBA();
                inClump.insert(addTAO);
                for (int i = 0; i < csdCurrent.size(); ++i) {
                    if (size >= binMins[i] && (size <= binMaxs[i] || i == csdCurrent.size() - 1)) {
                        csdCurrent[i] += ba;
                    }
                }
            } //While(seedNode.clumpSize() > targN)

            //Cut the adjacent TAOs
            seedNode.maxClumpSizeApprox() -= seedNode.clumpSize();
            for (size_t idx : inClump) {
                for (size_t adj : g.nodes[idx].adjList) {
                    lapis::lico::TaoNode& n = g.nodes[adj];
                    if (n.status == ns::off) {
                        n.cutTAO();
                        seedNode.maxClumpSizeApprox()--;
                    }
                }
            }
            if (intermediates) {
                TaoListMP out;
                for (size_t i = 0; i < g.nodes.size(); ++i) {
                    if (g.nodes[i].status == ns::on) {
                        out.taoVector.addFeature(rx.taos(taocrosswalk[i]));
                    }
                }
                std::cout << intermediatespath + std::to_string(step) + ".csv" << "\n";
                out.writeCsv(intermediatespath + std::to_string(step) + ".csv");
                std::cout << step << "\n";
                step++;

                for (int i = 0; i < csdCurrent.size(); ++i) {
                    targetstream << csdCurrent[i] << ",";
                }
                targetstream << "\n";

            }
        } //while (sofar < sample.size() && currentba < targetba)
        
        if (intermediates) targetstream.close();

        //Convert it all to a TaoList
        TaoListMP keep;
        TaoListMP cut;
        for (size_t i = 0; i < g.nodes.size(); ++i) {
            if (g.nodes[i].status == ns::on) {
                keep.taoVector.addFeature(rx.taos(taocrosswalk[i]));
            }
        }

        //If BA is still too low, add in cut TAOs
        if (currentba < targetba) {
            result = treatmentResult::cuttingFailure;
            for (size_t i : sample) {
                if (g.nodes[i].status != ns::on) {
                    keep.taoVector.addFeature(rx.taos(taocrosswalk[i]));
                    currentba += g.nodes[i].ba;
                    g.nodes[i].status = ns::on;
                    if (currentba >= targetba) {
                        break;
                    }
                }
            }
        }

        // Add the off nodes to the cut list (trees below dbhmin not in the graph)
        for (size_t i = 0; i < g.nodes.size(); ++i) {
            if (g.nodes[i].status != ns::on) {
                cut.taoVector.addFeature(rx.taos(taocrosswalk[i]));
            }
        }

        //If BA is *still* too low, add in TAOs below dbhmin
        for (size_t i = 0; i < rx.taos.size(); ++i) {
            //This loop currently doesn't use weighted randomness at all, and will prioritize the tallest TAOs.
            //This might be fine because the tallest ones will be the ones which are barely below dbhmin
            if (rx.taos.dbh(i) < dbhMin) {
                if (currentba < targetba) {
                    result = treatmentResult::cuttingFailure;
                    keep.taoVector.addFeature(rx.taos(i));
                    currentba += rx.taos.dbh(i) * rx.taos.dbh(i) / 2. / 2. / 100. / 100. * M_PI;
                    if (currentba >= targetba) {
                        break;
                    }
                }
                else {
                    cut.taoVector.addFeature(rx.taos(i));
                }
            }
        }
        return std::tuple(keep, cut, result);
    }

    namespace weightedRandom {
        Weights getWeightedOrder(lapis::lico::GraphLico& g, rxtools::StructureSummary targets,
            double area, std::default_random_engine& dre, std::pair<double, double> wlim, int nw) {

            //get backbone count
            size_t nbb = 0;
            double babb = 0;
            size_t nfull = 0;
            double bafull = 0;
            for (size_t i = 0; i < g.nodes.size(); ++i) {
                double thisba = g.nodes[i].dbh / 2 / 100;
                thisba *= thisba;
                thisba *= M_PI;
                if (g.nodes[i].status == ns::on) {
                    nbb++;
                    babb += thisba;
                }
                nfull++;
                bafull += thisba;
            }
            size_t ntao = g.nodes.size();

            //Tph targets compared to area
            if (targets.tph * area > ntao)
                targets.tph = ntao / area; //if there simply aren't enough trees, reduce the tph target
            if (targets.tph * area < nbb)
                targets.tph = nbb / area; //if there's too many backbone trees, increase the target

            //Target qmd
            double qmdTarg = std::sqrt(targets.ba / targets.tph / 0.0000785);
            size_t ntaoTarg = (size_t)std::round(targets.tph * area);

            double fullqmd = std::sqrt(bafull / nfull / 0.0000785);

            double bestQMD = 0;
            std::vector<size_t> bestsample;
            std::vector<double> bestweights;

            uniform_random dist(0., 1.);

            for (double thisWeight = wlim.first; thisWeight < wlim.second; thisWeight += (wlim.second - wlim.first) / nw) {
                std::vector<double> weights(g.nodes.size(), 0.);
                for (size_t i = 0; i < g.nodes.size(); ++i) {
                    weights[i] = std::pow(g.nodes[i].dbh, thisWeight);
                }
                size_t thisntao = nbb;
                double thisba = babb;
                std::vector<weightedNode> nodes(g.nodes.size());
                std::vector<size_t> sample;
                totalWeight(nodes, weights); //populate the tree with the info it needs to make selecting a node randomly be log(n)
                
                while (thisntao < ntaoTarg && sample.size() < nodes.size()) {
                    auto x = weightedSampleWithoutReplace(nodes, weights, dist, dre);
                    sample.push_back(x);
                    thisntao++;
                    double taoba = g.nodes[sample[sample.size() - 1]].dbh / 2 / 100;
                    taoba *= taoba * M_PI;
                    thisba += taoba;
                }

                double thisqmd = std::sqrt(thisba / thisntao / 0.0000785);
                bool isbest = bestQMD == 0 || (std::abs(bestQMD / fullqmd - qmdTarg / fullqmd) > std::abs(thisqmd / fullqmd - qmdTarg / fullqmd));
                if (isbest) {
                    bestQMD = thisqmd;
                    while (sample.size() < nodes.size()) {
                        sample.push_back(weightedSampleWithoutReplace(nodes, weights, dist, dre));
                    }
                    bestsample = sample;
                    bestweights = weights;
                }
            }
            Weights w;
            w.sample = bestsample;
            w.weights = bestweights;
            return w;
        }

        Weights getWeightedOrderByRatios(lapis::lico::GraphLico& g, rxtools::StructureSummary targets, double area, std::default_random_engine& dre) {

            auto targettph = (targets.ba * area * 10000) / (targets.tph * area);
            auto targetcc =  (targets.ba * area * 10000) / (targets.cc * 10000 * area);

            std::vector<double> weights(g.nodes.size(), 0.);
            for (size_t i = 0; i < g.nodes.size(); ++i) {
                double thisba = g.nodes[i].dbh / 2 / 100;
                thisba *= thisba;
                thisba *= M_PI;
                weights[i] = std::abs(thisba - targettph) + std::abs(thisba / g.nodes[i].area - targetcc);
            }
            std::vector<weightedNode> nodes(g.nodes.size());
            totalWeight(nodes, weights);

            uniform_random dist(0., 1.);
            Weights w;
            while (w.sample.size() < nodes.size()) {
                w.sample.push_back(weightedSampleWithoutReplace(nodes, weights, dist, dre));
            }
            w.weights = weights;

            return w;
        }
    }
 
} //namespace rxtools