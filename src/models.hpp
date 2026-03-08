#pragma once

#ifndef rxtools_models_h
#define rxtools_models_h

#include "allometry.hpp"

namespace rxtools::allometry {
    // This model will create a linear model (perhaps with transform) from height to the specified response variable.
    class UnivariateLinearModel : public Model {
    public:
        struct Transform {
            std::string name;
            std::function<double(double)> applyX;
            std::function<double(double)> applyY;
            std::function<double(double)> inverseY;
        };

        struct Parameters {
            double slope;
            double intercept;
            double rsq;
            Transform transform;
            Parameters(); //Declared inline after transforms because it needs to reference transforms::none
            Parameters(const double& s, const double& i, const Transform& t, const double& r) :
                slope(s), intercept(i), rsq(r), transform(t) {}
        };

        Parameters parameters;
        double predict(double x, const lapis::Unit& thisUnit, const lapis::Unit& returnUnit) const;


    protected:
        virtual std::vector<Transform> candidateTransforms() const = 0;
        void fitBestFromCandidateTransforms(const std::vector<double>& y, const std::vector<double>& x);
        Parameters calcModel(std::vector<double> y, std::vector<double> x, const Transform& tr) const;
        virtual void print(std::ostream& os) const;

    };

    namespace transforms {
        inline const UnivariateLinearModel::Transform none{
            "None",
            [](double x) { return x; },
            [](double y) { return y; },
            [](double y) { return y; }
        };
        inline const UnivariateLinearModel::Transform square{
            "Square",
            [](double x) { return x; },
            [](double y) { return y * y; },
            [](double y) { return std::sqrt(y); }
        };
        inline const UnivariateLinearModel::Transform cube{
            "Cube",
            [](double x) { return x; },
            [](double y) { return y * y * y; },
            [](double y) { return std::cbrt(y); }
        };
        inline const UnivariateLinearModel::Transform power{
            "Power (log-log)",
            [](double x) { return std::log(x); },
            [](double y) { return std::log(y); },
            [](double y) { return std::exp(y); }
        };
        inline const UnivariateLinearModel::Transform sqrt{
            "Square-root",
            [](double x) { return x; },
            [](double y) { return std::sqrt(y); },
            [](double y) { return y * y; }
        };
        inline const UnivariateLinearModel::Transform curt{
            "Cube-root",
            [](double x) { return x; },
            [](double y) { return std::cbrt(y); },
            [](double y) { return y * y * y; }
        };
        inline const UnivariateLinearModel::Transform log{
            "Log",
            [](double x) { return x; },
            [](double y) { return std::log(y); },
            [](double y) { return std::exp(y); }
        };
    }

    inline UnivariateLinearModel::Parameters::Parameters()
        : slope(0), intercept(0), rsq(0), transform(transforms::none) {
    }

    class DbhModel : public UnivariateLinearModel {
    public:
        DbhModel() {};
        DbhModel(const double& slope, const double& intercept, const Transform& transform, const double& rsq,
            const lapis::Unit& heightUnit, const lapis::Unit& diameterUnit);
        DbhModel(const std::vector<double>& heights, const std::vector<double>& diameters,
            const lapis::Unit& heightUnit, const lapis::Unit& diameterUnit,
            std::optional<Transform> transform = std::nullopt);
    protected:
        std::vector<Transform> candidateTransforms() const override;
    };

    class CrownModel : public UnivariateLinearModel {
    public:
        CrownModel() {};
        CrownModel(const double& slope, const double& intercept, const Transform& transform, const double& rsq,
            const lapis::Unit& heightUnit, const lapis::Unit& crownUnit);
        CrownModel(const std::vector<double>& heights, const std::vector<double>& crowns,
            const lapis::Unit& heightUnit, const lapis::Unit& crownUnit,
            std::optional<Transform> transform = std::nullopt);
    protected:
        std::vector<Transform> candidateTransforms() const override;
    };

    class BiomassModel : public UnivariateLinearModel {
    public:
        BiomassModel() {};
        BiomassModel(const double& slope, const double& intercept, const Transform& transform, const double& rsq,
            const lapis::Unit& heightUnit, const lapis::Unit& biomassUnit);
        BiomassModel(const std::vector<double>& heights, const std::vector<double>& biomasses,
            const lapis::Unit& heightUnit, const lapis::Unit& biomassUnit,
            std::optional<Transform> transform = std::nullopt);
    protected:
        std::vector<Transform> candidateTransforms() const override;
    };


    struct FastFuels {
    public:
        bool init = false;

        FastFuels() {};
        FastFuels(const FIATreeList& treeList);

        //assumes meter as input and output?
        //i'm going to leave this unfinished and return to it if it ever shows up again... here's a useful snippet from an email
        //"We calculate CBH as CBH = HT * CR where CR is the Compacted Crown Ratio in the FIA database.
        //There's a brief description of the quantity on page 185 (3-16) in the FIA DB documentation."
        double predictCbh(const double& height);

        int assignSpecies(const double& x, const double& y);

    private:
        std::vector<int> spcdWeights;
        std::hash<double> spcdHash;
        CrownModel cr;
    };

} // namespace rxtools::allometry


#endif // !rxtools_models_h