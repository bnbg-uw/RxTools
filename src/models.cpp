#include "models.hpp"

namespace rxtools::allometry {

    //----------------------------
    // Univaraite Linear Model
    //----------------------------
    UnivariateLinearModel::UnivariateLinearModel(const double& slope, const double& intercept, const Transform& transform, const double& rsq,
        const lapis::LinearUnit& inputUnit, const lapis::LinearUnit& outputUnit) :
        parameters(slope, intercept, transform, rsq) {
        this->inputUnit = inputUnit;
        this->outputUnit = outputUnit;
    }

    UnivariateLinearModel::UnivariateLinearModel(const FIATreeList& treeList, const std::string& responseName, 
        const lapis::LinearUnit& responseUnit, const Transform& transform) {
        inputUnit = lapis::linearUnitPresets::internationalFoot;
        outputUnit = responseUnit;

        auto it = std::find(treeList.names.begin(), treeList.names.end(), responseName);
        ptrdiff_t responseIdx;

        if (it == treeList.names.end()) {
            throw(std::out_of_range("responseName is not in names of treeList"));
        }
        else {
            responseIdx = std::distance(treeList.names.begin(), it);
        }
        
        if (transform != Transform::Suggest) {
            parameters = calcModel(treeList.height, treeList.otherfields.at(responseIdx), transform);
        }
        else {
            for (int t = static_cast<int>(Transform::None); t != static_cast<int>(Transform::Suggest); ++t) {
                auto p = calcModel(treeList.height, treeList.otherfields.at(responseIdx), static_cast<Transform>(t));
                if (p.rsq > parameters.rsq) {
                    parameters = p;
                }
            }
        }
    }

    double UnivariateLinearModel::predict(double x, const lapis::LinearUnit& thisUnit, const lapis::LinearUnit& returnUnit) const {
        auto inConverter = lapis::LinearUnitConverter(thisUnit, inputUnit);
        auto outConverter = lapis::LinearUnitConverter(outputUnit, returnUnit);
        x = inConverter(x);
        if (parameters.transform == Transform::Power)
            x = std::log(x);

        double out = parameters.intercept + parameters.slope * x;
        if (parameters.transform == Transform::Square) {
            out = std::sqrt(out);
        }
        else if (parameters.transform == Transform::Cube) {
            out = std::cbrt(out);
        }
        else if (parameters.transform == Transform::Power) {
            out = std::exp(out);
        }
        else if(parameters.transform != Transform::None) {
            throw(std::out_of_range("Unknown how to back-transform for enum value " + std::to_string(static_cast<int>(parameters.transform)) + "."));
        }
        return outConverter(out);
    }

    void UnivariateLinearModel::print(std::ostream& os) const {
        std::string tString;
        if (parameters.transform == Transform::None) {
            tString = "None";
        }
        else if (parameters.transform == Transform::Square) {
            tString = "Square";
        }
        else if (parameters.transform == Transform::Cube) {
            tString = "Cube";
        }
        else if (parameters.transform == Transform::Power) {
            tString = "Power (log-log)";
        }
        else {
            tString = "Unknown transform, int value is: " + static_cast<int>(parameters.transform);
        }

        os << "UnivariateLinearModel of the form:\n\tY = " + std::to_string(parameters.slope) + "*X + " + std::to_string(parameters.intercept) +
            "\n\t R^2: " + std::to_string(parameters.rsq) + "\tTransform: " + tString;
    }

    UnivariateLinearModel::Parameters UnivariateLinearModel::calcModel(const std::vector<double>& y, const std::vector<double>& x, const Transform& tr) const {
        Parameters p;
        p.transform = tr;

        //this is a workaround to avoid copying memory if y isn't actually being transformed
        const std::vector<double>* ypo = &y;
        std::vector<double> tfy;
        const std::vector<double>* xpo = &x;
        std::vector<double> tfx;
        if (tr != Transform::None) {
            tfy.resize(y.size());
            ypo = &tfy;
            if (tr == Transform::Power) {
                tfx.resize(x.size());
                xpo = &tfx;
            }
            for (int i = 0; i < y.size(); ++i) {
                if (tr == Transform::Square) {
                    tfy[i] = y[i] * y[i];
                }
                else if (tr == Transform::Cube) {
                    tfy[i] = y[i] * y[i] * y[i];
                }
                else if (tr == Transform::Power) {
                    tfy[i] = std::log(y[i]);
                    tfx[i] = std::log(x[i]);
                }
                else {
                    throw(std::out_of_range("Transform enum value " + std::to_string(static_cast<int>(parameters.transform)) + " is undefined here."));
                }
            }
        }

        double n = (*xpo).size();
        double sumx = 0;
        double sumx2 = 0;
        double sumxy = 0;
        double sumy = 0;
        double sumy2 = 0;

        for (int i = 0; i < n; ++i) {
            sumx += (*xpo)[i];
            sumx2 += (*xpo)[i] * (*xpo)[i];
            sumxy += (*xpo)[i] * (*ypo)[i];
            sumy += (*ypo)[i];
            sumy2 += (*ypo)[i] * (*ypo)[i];
        }
        double denom = n * sumx2 - sumx * sumx;
        p.slope = (n * sumxy - sumx * sumy) / denom;
        p.intercept = (sumy * sumx2 - sumx * sumxy) / denom;


        std::vector<double> pred; pred.resize(n);
        double predsum = 0;
        for (int i = 0; i < n; ++i) {
            pred[i] = (*xpo)[i] * p.slope + p.intercept;
            predsum += pred[i];
        }
        double ymean = sumy / n;
        double ssres = 0;
        double sstot = 0;
        for (int i = 0; i < y.size(); ++i) {
            double frommean = (*ypo)[i] - ymean;
            double frommodel = (*ypo)[i] - pred[i];
            sstot += frommean * frommean;
            ssres += frommodel * frommodel;
        }
        p.rsq = 1 - ssres / sstot;
        return p;
    }

    //--------------------
    //     FastFuels
    //--------------------
    FastFuels::FastFuels(const FIATreeList& treeList) : cr(treeList, "CR", lapis::linearUnitPresets::internationalFoot) {
        auto it = std::find(treeList.names.begin(), treeList.names.end(), "SPCD");
        ptrdiff_t spcdIdx;

        if (it == treeList.names.end()) {
            throw(std::out_of_range("responseName is not in names of treeList"));
        }
        else {
            spcdIdx = std::distance(treeList.names.begin(), it);
        }

        //double will be cast to int implicitly
        std::unordered_map<int, int> spcdTable;
        for (auto v : treeList.otherfields.at(spcdIdx)) {
            spcdTable.emplace(v, 0);
            ++spcdTable[v];
        }
        for (auto& v : spcdTable) {
            v.second = (int)std::round((double)v.second / treeList.otherfields.at(spcdIdx).size() * 100);
            for (size_t i = 0; i < v.second; ++i) {
                spcdWeights.push_back(v.first);
            }
        }
        init = true;
    }

    double FastFuels::predictCbh(const double& height) {
        auto crown = cr.predict(height, lapis::linearUnitPresets::meter, lapis::linearUnitPresets::meter);
        throw(std::exception("Bryce never finished implementing this- see the comments in the code and finish this up!"));
        return crown * height / 100;
    }

    int FastFuels::assignSpecies(const double& x, const double& y) {
        return spcdWeights[spcdHash(x * y) % spcdWeights.size()];
    }
} // namespace rxtools::allometry
