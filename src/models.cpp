#include "models.hpp"

namespace rxtools::allometry {
    namespace transforms{
        const UnivariateLinearModel::Transform none{
            "None",
            [](double x) { return x; },
            [](double y) { return y; },
            [](double y) { return y; }
        };
        const UnivariateLinearModel::Transform square{
            "Square",
            [](double x) { return x; },
            [](double y) { return y * y; },
            [](double y) { return std::sqrt(y); }
        };
        const UnivariateLinearModel::Transform cube{ 
            "Cube",
            [](double x) { return x; },
            [](double y) { return y * y * y; },
            [](double y) { return std::cbrt(y); }
        };
        const UnivariateLinearModel::Transform power{
            "Power (log-log)",
            [](double x) { return std::log(x); },
            [](double y) { return std::log(y); },
            [](double y) { return std::exp(y); }
        };
        const UnivariateLinearModel::Transform sqrt{
            "Square-root",
            [](double x) { return x; },
            [](double y) { return std::sqrt(y); },
            [](double y) { return y * y; }
        };
        const UnivariateLinearModel::Transform curt{
            "Cube-root",
            [](double x) { return x; },
            [](double y) { return std::cbrt(y); },
            [](double y) { return y * y * y; } 
        };
        const UnivariateLinearModel::Transform log{
            "Log",
            [](double x) { return x; },
            [](double y) { return std::log(y); },
            [](double y) { return std::exp(y); }
        };
    }

    /*DbhModel::UnivariateLinearModel(const FIATreeList& treeList, const std::string& responseName,
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

        std::vector<double> response;
        for (size_t i = 0; i < treeList.otherfields.size(); ++i) {
            response.push_back(treeList.otherfields.at(i).at(responseIdx));
        }

        if (transform != Transform::Suggest) {
            parameters = calcModel(response, treeList.height, transform);
        }
        else {
            for (int t = static_cast<int>(Transform::None); t != static_cast<int>(Transform::Suggest); ++t) {
                auto p = calcModel(response, treeList.height, static_cast<Transform>(t));
                if (p.rsq > parameters.rsq) {
                    parameters = p;
                }
            }
        }
    }*/

    //----------------------------
    // Univariate Linear Model
    //----------------------------
    double UnivariateLinearModel::predict(double x, const lapis::LinearUnit& thisUnit, const lapis::LinearUnit& returnUnit) const {
        auto inConverter = lapis::LinearUnitConverter(thisUnit, inputUnit);
        auto outConverter = lapis::LinearUnitConverter(outputUnit, returnUnit);
        
        x = inConverter(x);
        x = parameters.transform.applyX(x);
        double out = parameters.intercept + parameters.slope * x;
        out = parameters.transform.inverseY(out);
        return outConverter(out);
    }

    void UnivariateLinearModel::fitBestFromCandidateTransforms(const std::vector<double>& y, const std::vector<double>& x) {
        parameters.rsq = -std::numeric_limits<double>::infinity();
        for (const auto& tr : candidateTransforms()) {
            auto p = calcModel(y, x, tr);
            if (p.rsq > parameters.rsq) {
                parameters = p;
            }
        }
    }

    UnivariateLinearModel::Parameters UnivariateLinearModel::calcModel(std::vector<double> y, std::vector<double> x, const Transform& tr) const {
        Parameters p;
        p.transform = tr;

        for (int i = 0; i < y.size(); ++i) {
            x.at(i) = tr.applyX(x.at(i));
            y.at(i) = tr.applyY(y.at(i));
        }

        size_t n = x.size();
        double sumx = 0;
        double sumx2 = 0;
        double sumxy = 0;
        double sumy = 0;
        double sumy2 = 0;

        for (size_t i = 0; i < n; ++i) {
            sumx += x[i];
            sumx2 += x[i] * x[i];
            sumxy += x[i] * y[i];
            sumy += y[i];
            sumy2 += y[i] * y[i];
        }
        double denom = n * sumx2 - sumx * sumx;
        p.slope = (n * sumxy - sumx * sumy) / denom;
        p.intercept = (sumy * sumx2 - sumx * sumxy) / denom;


        std::vector<double> pred; pred.resize(n);
        double predsum = 0;
        for (int i = 0; i < n; ++i) {
            pred[i] = x[i] * p.slope + p.intercept;
            predsum += pred[i];
        }
        double ymean = sumy / n;
        double ssres = 0;
        double sstot = 0;
        for (int i = 0; i < y.size(); ++i) {
            double frommean = y[i] - ymean;
            double frommodel = y[i] - pred[i];
            sstot += frommean * frommean;
            ssres += frommodel * frommodel;
        }
        p.rsq = 1 - ssres / sstot;
        return p;
    }

    void UnivariateLinearModel::print(std::ostream& os) const {
        os << "UnivariateLinearModel of the form:\n\tY = " + std::to_string(parameters.slope) + "*X + " + std::to_string(parameters.intercept) +
            "\n\t R^2: " + std::to_string(parameters.rsq) + "\tTransform: " + parameters.transform.name;
    }

    //--------------------
    //     DbhModel
    //--------------------
    DbhModel::DbhModel(const double& slope, const double& intercept, const Transform& transform, const double& rsq,
        const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& diameterUnit) {
        parameters = Parameters(slope, intercept, transform, rsq);
        this->inputUnit = heightUnit;
        this->outputUnit = diameterUnit;
    }

    DbhModel::DbhModel(const std::vector<double>& heights, const std::vector<double>& diameters,
        const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& diameterUnit,
        std::optional<Transform> transform)
    {
        this->inputUnit = heightUnit;
        this->outputUnit = diameterUnit;
        if (transform.has_value()) {
            parameters = calcModel(diameters, heights, transform.value());
        }
        else {
            fitBestFromCandidateTransforms(diameters, heights);
        }
    }

    std::vector<UnivariateLinearModel::Transform> DbhModel::candidateTransforms() const {
        return { transforms::none, transforms::square, transforms::cube, transforms::power };
    }

    //--------------------
    //     CrownModel
    //--------------------
    CrownModel::CrownModel(const double& slope, const double& intercept, const Transform& transform, const double& rsq,
        const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& crownUnit) {
        parameters = Parameters(slope, intercept, transform, rsq);
        this->inputUnit = heightUnit;
        this->outputUnit = crownUnit;
    }

    CrownModel::CrownModel(const std::vector<double>& heights, const std::vector<double>& crowns,
        const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& crownUnit,
        std::optional<Transform> transform)
    {
        this->inputUnit = heightUnit;
        this->outputUnit = crownUnit;
        if (transform.has_value()) {
            parameters = calcModel(crowns, heights, transform.value());
        }
        else {
            fitBestFromCandidateTransforms(crowns, heights);
        }
    }

    std::vector<UnivariateLinearModel::Transform> CrownModel::candidateTransforms() const {
        return { transforms::none, transforms::square, transforms::cube, transforms::power };
    }

    //--------------------
    //     BiomassModel
    //--------------------
    BiomassModel::BiomassModel(const double& slope, const double& intercept, const Transform& transform, const double& rsq,
        const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& biomassUnit) {
        parameters = Parameters(slope, intercept, transform, rsq);
        this->inputUnit = heightUnit;
        this->outputUnit = biomassUnit;
    }

    BiomassModel::BiomassModel(const std::vector<double>& heights, const std::vector<double>& biomasses,
        const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& biomassUnit,
        std::optional<Transform> transform)
    {
        this->inputUnit = heightUnit;
        this->outputUnit = biomassUnit;
        if (transform.has_value()) {
            parameters = calcModel(biomasses, heights, transform.value());
        }
        else {
            fitBestFromCandidateTransforms(biomasses, heights);
        }
    }

    std::vector<UnivariateLinearModel::Transform> BiomassModel::candidateTransforms() const {
        return { transforms::none, transforms::sqrt, transforms::curt, transforms::log };
    }


    //--------------------
    //     FastFuels
    //--------------------
    FastFuels::FastFuels(const FIATreeList& treeList) {
        auto it = std::find(treeList.names.begin(), treeList.names.end(), "CR");
        ptrdiff_t crIdx;
        if (it == treeList.names.end()) {
            throw(std::out_of_range("\"CR\" is not in names of treeList"));
        }
        else {
           crIdx = std::distance(treeList.names.begin(), it);
        }
        auto crowns = treeList.otherfields.at(crIdx);
        cr = CrownModel(treeList.height, crowns, lapis::linearUnitPresets::internationalFoot, lapis::linearUnitPresets::internationalFoot);
        
        it = std::find(treeList.names.begin(), treeList.names.end(), "SPCD");
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
            spcdTable.emplace((int)v, 0);
            ++spcdTable[(int)v];
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
