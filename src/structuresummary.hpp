#pragma once

#ifndef rxtools_structuresummary_h
#define rxtools_structuresummary_h

#include "rxtools_pch.hpp"
#include "Raster.hpp"
#include "allometry.hpp"
#include "GraphLico.hpp"

namespace rxtools {
    namespace bg = boost::geometry;

    struct StructureSummary : public bg::model::point<double, 4, bg::cs::cartesian> {
        double ba;
        double tph;
        double mcs;
        double osi;
        double cc;
        std::vector<double> csd { 0,0,0,0, 0, 0  };
        std::vector<int> binMins{ 1,2,5,10,15,30 };
        std::vector<int> binMaxs{ 1,4,9,14,29,35 };

        StructureSummary() : ba(0), tph(0), mcs(0), osi(0), cc(0) {
            csd = std::vector<double>{ mcsTable[0][2], mcsTable[0][3], mcsTable[0][4], mcsTable[0][5], mcsTable[0][6], mcsTable[0][7] };
        };

        StructureSummary(double ba, double tph, double mcs, double osi, double cc) : ba(ba), tph(tph), mcs(mcs), osi(osi), cc(cc) {
            if (mcs < 1) mcs = 1;
            for (int i = 0; i < 28; ++i) {
                if (mcs >= mcsTable[i][0] && (mcs < mcsTable[i][1] || mcsTable[i][1] == 200))
                    csd = std::vector<double>{ mcsTable[i][2], mcsTable[i][3], mcsTable[i][4], mcsTable[i][5], mcsTable[i][6], mcsTable[i][7] };
            }
        };

        StructureSummary(double xba, double xtph, double xmcs, double xosi, double xcc, std::vector<double> xcsd) : ba(xba), tph(xtph), mcs(xmcs), osi(xosi), cc(xcc), csd(xcsd) {

            if (csd.size() != binMins.size())
                throw std::invalid_argument("xcsd should be of size 6, when not specifying bins");
        };

        StructureSummary(double ba, double tph, double mcs, double osi, double cc, std::vector<double> csd, std::vector<int> binMins, std::vector<int> binMaxs) :
            ba(ba), tph(tph), mcs(mcs), osi(osi), csd(csd), cc(cc), binMins(binMins), binMaxs(binMaxs) {
            if (csd.size() != binMins.size() || csd.size() != binMaxs.size())
                throw std::invalid_argument("csd, binmins, and binmaxs should all have equal size.");
        };

        StructureSummary(const TaoListMP& taos,
            const lapis::Alignment& unitAlign,
            double areaHa,
            double osi = -1
            ) : osi(osi)
        {
            lapis::lico::GraphLico g{ unitAlign };
            g.addDataset(taos.taoVector, taos.nodeFactory, lapis::lico::NodeStatus::on);

            try {
                mcs = 0;
                cc = 0;
                ba = 0;
                for (size_t i = 0; i < g.nodes.size(); ++i) {
                    mcs += g.nodes.at(i).clumpSize();
                    cc += g.nodes.at(i).area;
                    ba += g.nodes.at(i).ba;
                }
                mcs /= static_cast<double>(g.nodes.size());
                cc = cc / 10000. / areaHa;
                ba = ba / areaHa;
            }
            catch (std::bad_alloc e) {
                std::cout << "ba/mcs\n";
                throw e;
            }
            catch (std::exception e) {
                std::cout << e.what();
                throw e;
            }

            tph = g.nodes.size() / areaHa;

            try {
                for (int i = 0; i < g.nodes.size(); i++) { //TODO: don't hardcode these bins
                    auto clsz = g.nodes[i].clumpSize();
                    if (clsz == 1)
                        csd[0]++;
                    else if (clsz < 5)
                        csd[1]++;
                    else if (clsz < 10)
                        csd[2]++;
                    else if (clsz < 15)
                        csd[3]++;
                    else if (clsz < 31)
                        csd[4]++;
                    else
                        csd[5]++;
                }

                for (int i = 0; i < csd.size(); i++)
                    csd[i] /= static_cast<double>(g.nodes.size());
            }
            catch (std::bad_alloc e) {
                std::cout << "csd\n";
                throw e;
            }
        }

        double operator[](int idx) const {
            if (idx == 0) return ba;
            if (idx == 1) return tph;
            if (idx == 2) return mcs;
            if (idx == 3) return osi;
            if (idx == 4) return cc;
            throw std::out_of_range("Index out of bounds");
        }

    private:
        //MCS look up table to clump size distribution from Sean Jeronimo.
        //Columns: Mcs Min, Mcs Max, 1,  2to4, 5to9, 10to14, 15to29, 30+
        double mcsTable[28][8] = {
            {1,    1.25, 0.94, 0,    0,    0,    0.06, 0},
            {1.25, 1.5,  0.69, 0,    0,    0.02, 0.29, 0},
            {1.5,  1.75, 0.6,  0,    0.02, 0.01, 0.37, 0},
            {1.75, 2,    0.53, 0,    0.03, 0.05, 0.39, 0},
            {2,    2.25, 0.49, 0,    0.08, 0.02, 0.41, 0},
            {2.25, 2.5,  0.46, 0.01, 0.1,  0.04, 0.39, 0},
            {2.5,  2.75, 0.43, 0.02, 0.14, 0.03, 0.38, 0},
            {2.75, 3,    0.41, 0.03, 0.15, 0.04, 0.36, 0.01},
            {3,    3.25, 0.37, 0.03, 0.16, 0.03, 0.4,  0.01},
            {3.25, 3.5,  0.38, 0.04, 0.17, 0.01, 0.38, 0.02},
            {3.5,  3.75, 0.35, 0.06, 0.18, 0,    0.4,  0.01},
            {3.75, 4,    0.36, 0.05, 0.18, 0.01, 0.37, 0.03},
            {4,    4.25, 0.33, 0.05, 0.21, 0.03, 0.35, 0.03},
            {4.25, 4.5,  0.33, 0.05, 0.2,  0,    0.37, 0.05},
            {4.75, 5,    0.3,  0.07, 0.2,  0.01, 0.36, 0.06},
            {5,    6,    0.29, 0.09, 0.21, 0.03, 0.31, 0.07},
            {6,    7,    0.24, 0.1,  0.21, 0.02, 0.34, 0.09},
            {7,    8,    0.25, 0.08, 0.19, 0.02, 0.31, 0.15},
            {8,    9,    0.25, 0.1,  0.18, 0.04, 0.31, 0.12},
            {9,    10,   0.25, 0.08, 0.18, 0.05, 0.27, 0.17},
            {10,   15,   0.22, 0.09, 0.17, 0.11, 0.27, 0.14},
            {15,   20,   0.24, 0.04, 0.13, 0.17, 0.27, 0.15},
            {20,   25,   0.19, 0.07, 0.14, 0.22, 0.27, 0.11},
            {25,   50,   0.15, 0.08, 0.15, 0.28, 0.24, 0.1},
            {50,   75,   0.13, 0.09, 0.13, 0.37, 0.19, 0.09},
            {75,   100,  0.15, 0.04, 0.11, 0.41, 0.2,  0.09},
            {100,  200,  0.1,  0.03, 0.08, 0.59, 0.12, 0.08}
        };
    };

    inline double calcOsi(lapis::Raster<int> chm, lapis::coord_t heightCutoff, lapis::coord_t coreDist) {
        lapis::Raster<lapis::coord_t> edt = lapis::euclideanDistanceTransform(chm, [heightCutoff](lapis::coord_t v) { return v >= heightCutoff; });
        lapis::Raster<char> isCoreGap = edt >= coreDist;
        
        int osinum = 0;
        int osiden = 0;
        for (lapis::cell_t i = 0; i < isCoreGap.ncell(); ++i) {
            if (isCoreGap[i].has_value()) {
                if (isCoreGap[i].value()) {
                    osinum++;
                }
                osiden++;
            }
        }
        return (static_cast<double>(osinum) / static_cast<double>(osiden)) * 100.;
    }
}

namespace boost {
    namespace geometry {
        namespace traits {

            template<>
            struct tag<rxtools::StructureSummary>
            {
                typedef point_tag type;
            };

            template<>
            struct coordinate_type<rxtools::StructureSummary>
            {
                typedef double type;
            };

            template<>
            struct coordinate_system<rxtools::StructureSummary>
            {
                typedef cs::cartesian type;
            };

            template<>
            struct dimension<rxtools::StructureSummary> : boost::mpl::int_<4> {};

            template<std::size_t K>
            struct access<rxtools::StructureSummary, K>
            {
                static inline double get(rxtools::StructureSummary const& p)
                {
                    return p.template get<K>();
                }

                static inline void set(rxtools::StructureSummary& p, double const& v)
                {
                    p.template set<K>(v);
                }
            };
        }
    }
} //namespace rxtools
#endif //rxtools_structuresummary_h