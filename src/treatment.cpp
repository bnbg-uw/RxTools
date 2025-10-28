#include "licosim/treatment.hpp"

namespace licosim {

    using ns = lico::NodeStatus;
    
    lico::TaoList Treatment::doTreatment(RxUnit rx, double dbhMin, double dbhMax) {
        //std::cout << "Starting Treatment:\n";
        std::vector<treatmentDecision> treeStatus{ (unsigned __int64) rx.taos.size(), td::undecided };
        std::vector<double> dbh = rx.dbhFunc(rx.taos.height());

        //Assign weights, and give initial status of retain or cut to small or big dbh's.
        std::vector<double> weights;
        //std::cout << "\tcalc w started\n";

        //mark backbone trees as retain and too-small trees as cut
        for (index_t i = 0; i < rx.taos.size(); ++i) {
            if (dbh[i] < dbhMin) {
                treeStatus[i] = td::cut;
            }
            else if (dbh[i] > dbhMax) {
                treeStatus[i] = td::retain;
            }
        }

        //assign weights to all remaining trees
        double w = calcW(dbh, treeStatus, rx.targetStructure, rx.areaHa);
        for (int i = 0; i < rx.taos.size(); ++i) {
            if (treeStatus[i] == td::cut) {
                weights.push_back(0);
            }
            else {
                weights.push_back(std::pow(dbh[i], w));
            }
        }

        std::vector<double> csdTargs = rx.targetStructure.csd;

        //std::cout << "\tGoing into treatment loop:\n";

        //Ba is most important and therefore our stopping criterion.
        //This loop happens at most ntao times
        while (calcBA(dbh, treeStatus, rx.areaHa) < rx.targetStructure.ba) {
            // Pick a seed tree, or break from the loop if there are no trees left
            index_t seed;
            bool anycandidate = false;
            for (auto& w : weights) {
                if (w != 0) {
                    anycandidate = true;
                    break;
                }
            }

            if (!anycandidate)
                break;
            else {
                //seed = seedCandidates[randomIdxFromWeights(weights, dre)]; This was incorrect to begin with...
                seed = randomIdxFromWeights(weights, dre);
            }

            // Get trees in a clump with seed tree (exclude trees that have already been cut)
            lico::TaoList seedTl;
            std::vector<int> outerIdx; // convert between our subsetted list and the total taolist.
            int seedTlIdx; //This will be the idx of our seed tree in this subsetted tao list.
            //Start with all uncut trees
            //this loop is O(n)
            for (index_t i = 0; i < rx.taos.size(); i++) {
                if (weights[i]) { //TODO: perhaps turn this into a check for the tree being marked as cut. Will first verify if that is identical logic
                    seedTl.addTAO(rx.taos[i]);
                    outerIdx.push_back(i);
                    if (i == seed)
                        seedTlIdx = outerIdx.size() - 1;
                } 
            }
            auto sd = lico::SparseDistMatrix();
            sd.nearDist(seedTl.x(), seedTl.y(), seedTl.crown(), false);
            lico::ClumpInfo clumps = sd.getClumps(seedTl.size());

            //Subset again to only trees that clump with our seed tree.
            std::vector<index_t> clOuterIdx; //This is a list of positions of trees in our focal clump in the total taolist.
            for (index_t i = 0; i < seedTl.size(); i++) {
                if (clumps.clumpID[i] == clumps.clumpID[seedTlIdx]) {
                    clOuterIdx.push_back(outerIdx[i]);
                }
            }

            int targN = std::numeric_limits<int>::max(); // Eventual target clump size.
            std::vector<double> abstractNeed = calcNeed(rx.taos, dbh, treeStatus, rx.areaHa, rx.targetStructure.ba, csdTargs); //Bins needed by currently retained trees

            // Can make this a weight random sample for the binIdx if it fails to distribute across clumps on real life data
            std::vector<double> binWeights = abstractNeed;
            double bwSum = 0;
            for (int i = 0; i < rx.targetStructure.binMins.size(); i++) {
                /*if ((binMins[i] <= clOuterIdx.size()) & (abstractNeed[i] > 0)) {
                    binIdx = i;
                    break;
                }*/
                if ((rx.targetStructure.binMins[i] > clOuterIdx.size()) || (abstractNeed[i] < 0)) { 
                    binWeights[i] = 0; //If a bin is larger than the focal clump, or if no more clumps of that size are needed, throw that bin's need out for now
                }
                bwSum += binWeights[i];

            }

            //Pick the size in the bin to pare the clump down to, randomly.
            if (bwSum) {
                int binIdx = randomIdxFromWeights(binWeights, dre);
                std::uniform_int_distribution<int> d(rx.targetStructure.binMins[binIdx], rx.targetStructure.binMaxs[binIdx]);
                targN = d(dre);
            }

            if (targN > clOuterIdx.size()) {
                targN = clOuterIdx.size();
            }
            
            //Maybe it would be better to test for obligates before getting the targ N?
            std::unordered_set<index_t> finalClIdx; //Indexes in the total tao list that we will ultimately retain.
            obligateIdx(rx.taos, treeStatus, clOuterIdx, seed, finalClIdx); //Add seed tree and obligates.
            //finalClIdx is now a set of the minimal clump it is possible for this TAO to be a part of, given existing retains

            //Keep adding trees until we reach desired size.
            while (finalClIdx.size() < targN) {
                std::unordered_set<index_t> adj;
                //IDK if this nested loop can be optimized...
                // Get all trees adjacent to current clump.
                for (index_t i : finalClIdx) {
                    for (index_t j : clOuterIdx) {
                        if (finalClIdx.find(j) == finalClIdx.end()) {
                            double dsq = (rx.taos.x()[i] - rx.taos.x()[j]) * (rx.taos.x()[i] - rx.taos.x()[j]) + (rx.taos.y()[i] - rx.taos.y()[j]) - (rx.taos.y()[i] * rx.taos.y()[j]);
                            if (dsq <= (rx.taos.crown()[i] + rx.taos.crown()[j]) * (rx.taos.crown()[i] + rx.taos.crown()[j])) { //If the 2 trees are adj
                                //And the new tree is not already cut (if it is already in the adj table, it won't be added again.
                                if (treeStatus[j] != td::cut)
                                    adj.emplace(j); //If this tree is not already in the clump, intersects with the clump, and is not already marked as cut, add it to the clump
                            }
                        }
                    }
                }

                //Get the obligations for each of the potential adjacent trees
                std::vector<std::unordered_set<index_t>> ob;
                for (auto idx : adj) {
                    std::unordered_set<index_t> keepSet = finalClIdx;
                    obligateIdx(rx.taos, treeStatus, clOuterIdx, idx, keepSet);
                    ob.push_back(keepSet);
                }
                ob.push_back(finalClIdx);

                //Get which choice(s) get us as close as possible to targN
                std::vector<int> delt;
                for (auto& set : ob)
                    delt.push_back(std::abs(targN - (index_t)set.size()));

                //Get the index(es) that get us the closest to targetN
                std::vector<int> minIdx;
                int min = std::numeric_limits<int>::max();
                for (index_t i = 0; i < delt.size(); i++) {
                    if (delt[i] < min) {
                        min = delt[i];
                        minIdx.clear();
                        minIdx.push_back(i);
                    }
                    else if (delt[i] == min)
                        minIdx.push_back(i);
                }

                index_t idx;
                if (minIdx.size() > 1) {
                    //Pick one ranomly-ish.
                    std::vector<double> w(minIdx.size(), 0);
                    for (int i = 0; i < minIdx.size();  i++) {
                        for (index_t treeIdx : ob[minIdx[i]])
                            w[i] += weights[treeIdx];
                        w[i] /= ob[minIdx[i]].size();
                    }

                    idx = randomIdxFromWeights(w, dre);
                }
                else
                    idx = minIdx[0];
                finalClIdx = ob[idx];

                //Break out if we picked the existing clump, and it was the only choice (aka future adjacent trees will be the same and we will be in an infinite loop).
                if (minIdx.size() <= 1 && (idx == ob.size() - 1))
                    break;

            }

            //Finally, update treeStatus
            for (index_t i : finalClIdx) {
                treeStatus[i] = td::retain;
                weights[i] = 0;
            }

            //Cut adjacent trees that would form a clump with a retention tree, to lock this clump at it's current size.
            std::unordered_set<int> adj;
            //IDK if this nested loop can be optimized...
            for (index_t i : finalClIdx) {
                for (index_t j : clOuterIdx) {
                    //If the new tree is not already in the clump
                    if (finalClIdx.find(j) == finalClIdx.end()) {
                        // and the 2 trees are adj
                        double dsq = (rx.taos.x()[i] - rx.taos.x()[j]) * (rx.taos.x()[i] - rx.taos.x()[j]) + (rx.taos.y()[i] - rx.taos.y()[j]) - (rx.taos.y()[i] * rx.taos.y()[j]);
                        if (dsq <= (rx.taos.crown()[i] + rx.taos.crown()[j]) * (rx.taos.crown()[i] + rx.taos.crown()[j]))
                            adj.emplace(j);
                    }
                }
            }
            for (index_t a : adj) {
                treeStatus[a] = td::cut;
                weights[a] = 0;
            }
            //std::cout << "\t\t" << calcBA(dbh, treeStatus, rx.areaHa) << "/" << rx.targetStructure.ba <<
                //" (Added " << finalClIdx.size() << " trees of " << clOuterIdx.size() << ")\n";
        }

        //std::cout << "Exiting treatment loop\n";

        //If our smart thinning algorithm cut too many trees that it backed itself into a corner and couldn't meet the BA target
        //Then add trees back randomly until we hit the target.
        std::vector<index_t> cutTreeIdx;
        std::vector<double> cutTreeWeights;
        for (index_t i = 0; i < treeStatus.size(); i++) {
            if (treeStatus[i] == td::cut) {
                cutTreeIdx.push_back(i);
                cutTreeWeights.push_back(std::pow(dbh[i], w));
            }
        }
            
        index_t available = cutTreeIdx.size();
        while (calcBA(dbh, treeStatus, rx.areaHa) < rx.targetStructure.ba && available) {
            index_t idx = randomIdxFromWeights(cutTreeWeights, dre);
            treeStatus[cutTreeIdx[idx]] = td::retain;
            cutTreeWeights[idx] = 0;
            available--;
        }

        lico::TaoList out;
        for (index_t i = 0; i < treeStatus.size(); ++i) {
            if (treeStatus[i] == td::retain)
                out.addTAO(rx.taos[i]);
        }



        return out;
    }

    double Treatment::calcW(std::vector<double> dbh, std::vector<treatmentDecision> treeStatus, StructureSummary targets,
                            double area, std::pair<double, double> wlim, int nw) {
        int n = 0; // N trees above dbhmin
        int nbb = 0; // N backbone (trees greater than dbh max cut).
        for (index_t i = 0; i < treeStatus.size(); ++i) {
            if (treeStatus[i] != td::cut) n++;
            if (treeStatus[i] == td::retain) nbb++;
        }

        //Tph targets compared to area
        if (targets.tph * area > n)
            targets.tph = n / area; //if there simply aren't enough trees, reduce the tph target
        if (targets.tph * area < nbb)
            targets.tph = nbb / area; //if there's too many backbone trees, increase the target

        //Current qmd
        double qmd0 = std::sqrt(calcBA(dbh, treeStatus, area, true) / calcTPH(treeStatus, area, true) / 0.0000785); //calculate qmd using both retain and undecided

        //Range of w's to try. nw total attempts, spread uniformly across wlim
        std::vector<double> ws;
        double step = (wlim.second - wlim.first) / (nw - 1);
        for (int i = 0; i < nw; ++i) {
            ws.push_back(wlim.first + step * i);
        }

        //Target qmd
        double qmdTarg = std::sqrt(targets.ba / targets.tph / 0.0000785);

        std::vector<double> qmd; // Vector of trial qmd's
        
        //For priming out sample vector
        std::vector<double> retainOnly;
        for (index_t i = 0; i < dbh.size(); i++) {
            if (treeStatus[i] == td::retain)
                retainOnly.push_back(dbh[i]);
        }

        //Sample options
        std::vector<double> undecidedOnly;
        for (index_t i = 0; i < dbh.size(); i++) {
            if (treeStatus[i] == td::undecided)
                undecidedOnly.push_back(dbh[i]);
        }

        //looping over plausible weights
        for (double w : ws) {
            std::vector<double> samp = retainOnly;
            
            std::vector<double> localWeights;
            for (int i = 0; i < undecidedOnly.size(); i++)
                localWeights.push_back(std::pow(undecidedOnly[i], w));

            //Unsigned int? Or do I need to check for edge cases here??
            index_t nSamp = targets.tph * area - retainOnly.size();

            //sample randomly until your tph is correct
            for (index_t i = 0; i < nSamp; i++) {
                int idx = randomIdxFromWeights(localWeights, dre);
                samp.push_back(undecidedOnly[idx]);
                localWeights[idx] = 0;
            }

            //calculate qmd for this 
            double thisQmd = 0;
            for (int i = 0; i < samp.size(); i++)
                thisQmd += samp[i] * samp[i];
            qmd.push_back(thisQmd / samp.size());
        }

        // Per Sean this should become a splined regression to avoid outliers.
        std::vector<double> qmdRatio;
        for (double q : qmd)
            qmdRatio.push_back(q / qmd0);

        //Find the weight with the qmd closest to the target
        int minIdx = 0;
        double min = std::numeric_limits<double>::max();
        for (int i = 0; i < qmdRatio.size(); i++) {
            double dist = qmdRatio[i] - qmdTarg / qmd0;
            if (std::abs(dist) < min) {
                min = dist;
                minIdx = i;
            }
        }
        return ws[minIdx];
    }

    std::vector<double> Treatment::calcNeed(lico::TaoList tl, std::vector<double> dbh, std::vector<treatmentDecision> treeStatus,
                                            double area, double baTargs, std::vector<double> csdTargs) {
        lico::TaoList subTl;
        for (index_t i = 0; i < treeStatus.size(); i++) {
            if (treeStatus[i] == td::retain)
                subTl.addTAO(tl[i]);
        }
        auto props = calcClumpProps(subTl);
        auto ba = calcBA(dbh, treeStatus, area);
        for (index_t i = 0; i < props.size(); i++) {
            props[i] = baTargs * csdTargs[i] - props[i] * ba;
        }
        return props;
    }

    double Treatment::calcBA(std::vector<double> dbh, std::vector<treatmentDecision> treeStatus, double area, bool countUndecided) {
        if (treeStatus.size() != dbh.size()) {
            throw std::invalid_argument("dbh and treestatus must have equal lengths");
        }
        double outBA = 0;
        for (int i = 0; i < dbh.size(); i++) {
            if (treeStatus[i] == td::retain || (countUndecided && treeStatus[i] == td::undecided))
                outBA += M_PI * dbh[i] * dbh[i] / 4.;
        }
        return outBA / 10000 / area;
    }

    double Treatment::calcTPH(std::vector<treatmentDecision> treeStatus, double area, bool countUndecided) {
        int n = 0;
        for (int i = 0; i < treeStatus.size(); i++) {
            if (treeStatus[i] == td::retain || (countUndecided && treeStatus[i] == td::undecided))
                n++;
        }
        return (double)n / area;
    }

    std::vector<double> Treatment::calcClumpProps(lico::TaoList tl) {
        std::vector<double> out(6, 0);
        if (tl.size() > 0) {
            auto sd = lico::SparseDistMatrix();
            sd.nearDist(tl.x(), tl.y(), tl.crown(), false);
            lico::ClumpInfo clumps = sd.getClumps(tl.size());

            for (int i = 0; i < clumps.clumpSize.size(); i++) { //TODO: don't hardcode these bins
                if (clumps.clumpSize[i] == 1)
                    out[0]++;
                else if (clumps.clumpSize[i] < 5)
                    out[1]++;
                else if (clumps.clumpSize[i] < 10)
                    out[2]++;
                else if (clumps.clumpSize[i] < 15)
                    out[3]++;
                else if (clumps.clumpSize[i] < 31)
                    out[4]++;
                else
                    out[5]++;
            }

            for (int i = 0; i < out.size(); i++)
                out[i] /= clumps.clumpSize.size();
        }
        return out;
    }

    /*std::vector<double> Treatment::getClumpProps(double mcs) {
        std::filebuf fb;
        if (!fb.open(mcsPropsDir, std::ios::in)) 
            throw std::runtime_error("Cannot open internal reference table.");
        std::istream is{ &fb };

        readCSVLine(is); //skip colnames.
        while (!is.eof()) {
            auto row = readCSVLine(is);
            if (row.size() <= 1)
                continue;
            double min = std::stod(row[0]);
            double max = std::stod(row[1]);

            //if mcs between breaks or mcs is above the last min (200 is max in table, but we treat it as inf).
            if (mcs >= min && (mcs < max || max == 200))
                return std::vector<double> { std::stod(row[2]), std::stod(row[3]), std::stod(row[4]),
                                             std::stod(row[5]), std::stod(row[6]), std::stod(row[7]) };
        }
    }*/

    void Treatment::obligateIdx(lico::TaoList& tl, std::vector<treatmentDecision>& treeStatus, std::vector<index_t>& treeIdx, index_t& focalIdx, std::unordered_set<index_t>& keepSet) {
        std::vector<int> testIdx; //These are idx in the total dataset to test as focals.
        
        //Get all the trees that are not already known as keeps, but are obligate retzins and clump with the focal tree
        for (index_t i : treeIdx) {
            if (treeStatus[i] == td::retain && keepSet.find(i) == keepSet.end()) {
                double dsq = (tl.x()[focalIdx] - tl.x()[i]) * (tl.x()[focalIdx] - tl.x()[i]) + (tl.y()[focalIdx] - tl.y()[i]) * (tl.y()[focalIdx] - tl.y()[i]);
                if (dsq <= (tl.crown()[focalIdx] + tl.crown()[i]) * (tl.crown()[focalIdx] + tl.crown()[i]))
                    testIdx.push_back(i);
            }
        }

        keepSet.insert(focalIdx); //If a tree is focal we know we have to keep it.
        //Now check all these new trees for their obligates too. (We will add each of them to keepSet when doing this)
        for (index_t idx : testIdx)
            obligateIdx(tl, treeStatus, treeIdx, idx, keepSet);  //This can be optimized I believe... look into std::union_set or something maybe?
        // Since we pass keepSet by reference, it will be updated with the new taos as we go.
    }

    std::tuple<lico::TaoList, lico::TaoList, treatmentResult> Treatment::doTreatmentGraph(RxUnit rx, double dbhMin, double dbhMax, lico::unit_t maxCrown, bool intermediates, std::string intermediatespath) {
        auto binMins = rx.targetStructure.binMins;
        auto binMaxs = rx.targetStructure.binMaxs;

        spatial::Alignment a{ rx.unitMask.xmin(), rx.unitMask.ymin(),(spatial::rowcol_t)std::ceil((rx.unitMask.ymax() - rx.unitMask.ymin()) / maxCrown),
            (spatial::rowcol_t)std::ceil((rx.unitMask.xmax() - rx.unitMask.xmin()) / maxCrown),maxCrown,maxCrown };
        lico::GraphLico g{ a };

        std::vector<double> dbh = rx.dbhFunc(rx.taos.height());
        
        //absolute ba rather than per area
        double targetba = rx.targetStructure.ba * rx.areaHa;

        //current structure
        double currentba = 0;
        lico::index_t currentntao = 0;

        //vector to convert between taolist in rxunit 
        //and the potentially subset list used in the graph because of dbhmin
        std::vector<lico::index_t> taocrosswalk;

        //populate the graph and set up the crosswalk.
        for (lico::index_t i = 0; i < rx.taos.size(); ++i) {
            if (dbh[i] < dbhMin) {
                continue;
            }
            taocrosswalk.push_back(i);
            if (dbh[i] > dbhMax) {
                g.addTAO(rx.taos[i], dbh[i], ns::on);
                currentba += g.nodes[g.nodes.size() - 1].ba;
                currentntao++;
                continue;
            }
            g.addTAO(rx.taos[i], dbh[i]);
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
        for (lico::index_t i = 0; i < g.nodes.size(); ++i) {
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
        lico::index_t sofar = 0;
        std::vector<lico::index_t> priority(g.nodes.size());
        for (lico::index_t i = 0; i < sample.size(); ++i) {
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
            lico::index_t seedTao = sample[sofar];
            lico::TaoNode& seedNode = g.nodes[seedTao];
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
            int binIdx = randomIdxFromWeights(binWeights, dre);
            std::uniform_int_distribution<int> d(binMins[binIdx], binMaxs[binIdx]);
            int targN = d(dre);
            if (targN > maxClump) {
                targN = maxClump;
            }
            int doNotExceed = binMaxs[binIdx];

            std::unordered_set<lico::index_t> inClump{};
            //inClump.insert(seedNode.index);

            std::priority_queue<NodePriority, std::vector<NodePriority>, NodePriorityComparator> adjQueue;
            adjQueue.push(NodePriority{ priority[seedNode.index],seedNode.index });
            std::unordered_set<lico::index_t> considered;
            considered.insert(seedNode.index);

            while (seedNode.clumpSize() < targN) {
                //pop off of adjQueue until we get a TAO which doesn't take us over the bin threshold
                lico::index_t addTAO;
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
                for (lico::index_t adj : g.nodes[addTAO].adjList) {
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
                std::unordered_set<lico::index_t> adjClumps;
                for (lico::index_t adj : g.nodes[addTAO].adjList) {
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
            for (lico::index_t idx : inClump) {
                for (lico::index_t adj : g.nodes[idx].adjList) {
                    lico::TaoNode& n = g.nodes[adj];
                    if (n.status == ns::off) {
                        n.cutTAO();
                        seedNode.maxClumpSizeApprox()--;
                    }
                }
            }
            if (intermediates) {
                lico::TaoList out;
                for (lico::index_t i = 0; i < g.nodes.size(); ++i) {
                    if (g.nodes[i].status == ns::on) {
                        out.addTAO(rx.taos[taocrosswalk[i]]);
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
        lico::TaoList keep;
        lico::TaoList cut;
        for (lico::index_t i = 0; i < g.nodes.size(); ++i) {
            if (g.nodes[i].status == ns::on) {
                keep.addTAO(rx.taos[taocrosswalk[i]]);
            }
        }

        //If BA is still too low, add in cut TAOs
        if (currentba < targetba) {
            result = treatmentResult::cuttingFailure;
            for (lico::index_t i : sample) {
                if (g.nodes[i].status != ns::on) {
                    keep.addTAO(rx.taos[taocrosswalk[i]]);
                    currentba += g.nodes[i].ba;
                    g.nodes[i].status = ns::on;
                    if (currentba >= targetba) {
                        break;
                    }
                }
            }
        }

        // Add the off nodes to the cut list (trees below dbhmin not in the graph)
        for (lico::index_t i = 0; i < g.nodes.size(); ++i) {
            if (g.nodes[i].status != ns::on) {
                cut.addTAO(rx.taos[taocrosswalk[i]]);
            }
        }

        //If BA is *still* too low, add in TAOs below dbhmin
        for (lico::index_t i = 0; i < rx.taos.size(); ++i) {
            //This loop currently doesn't use weighted randomness at all, and will prioritize the tallest TAOs.
            //This might be fine because the tallest ones will be the ones which are barely below dbhmin
            if (dbh[i] < dbhMin) {
                if (currentba < targetba) {
                    result = treatmentResult::cuttingFailure;
                    keep.addTAO(rx.taos[i]);
                    currentba += dbh[i] * dbh[i] / 2. / 2. / 100. / 100. * M_PI;
                    if (currentba >= targetba) {
                        break;
                    }
                }
                else {
                    cut.addTAO(rx.taos[i]);
                }
            }
        }
        return std::tuple(keep, cut, result);
    }

    namespace weightedRandom {
        Weights getWeightedOrder(lico::GraphLico& g, StructureSummary targets,
            double area, std::default_random_engine& dre, std::pair<double, double> wlim, int nw) {

            //get backbone count
            lico::index_t nbb = 0;
            double babb = 0;
            lico::index_t nfull = 0;
            double bafull = 0;
            for (lico::index_t i = 0; i < g.nodes.size(); ++i) {
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
            lico::index_t ntao = g.nodes.size();

            //Tph targets compared to area
            if (targets.tph * area > ntao)
                targets.tph = ntao / area; //if there simply aren't enough trees, reduce the tph target
            if (targets.tph * area < nbb)
                targets.tph = nbb / area; //if there's too many backbone trees, increase the target

            //Target qmd
            double qmdTarg = std::sqrt(targets.ba / targets.tph / 0.0000785);
            lico::index_t ntaoTarg = std::round(targets.tph * area);

            double fullqmd = std::sqrt(bafull / nfull / 0.0000785);

            double bestQMD = 0;
            std::vector<lico::index_t> bestsample;
            std::vector<double> bestweights;

            uniform_random dist(0., 1.);

            for (double thisWeight = wlim.first; thisWeight < wlim.second; thisWeight += (wlim.second - wlim.first) / nw) {
                std::vector<double> weights(g.nodes.size(), 0.);
                for (lico::index_t i = 0; i < g.nodes.size(); ++i) {
                    weights[i] = std::pow(g.nodes[i].dbh, thisWeight);
                }
                index_t thisntao = nbb;
                double thisba = babb;
                std::vector<weightedNode> nodes(g.nodes.size());
                std::vector<lico::index_t> sample;
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

        Weights getWeightedOrderByRatios(lico::GraphLico& g, StructureSummary targets, double area, std::default_random_engine& dre) {

            auto targettph = (targets.ba * area * 10000) / (targets.tph * area);
            auto targetcc =  (targets.ba * area * 10000) / (targets.cc * 10000 * area);

            std::vector<double> weights(g.nodes.size(), 0.);
            for (lico::index_t i = 0; i < g.nodes.size(); ++i) {
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

    
}