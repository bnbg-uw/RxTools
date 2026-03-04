#pragma once

#ifndef rxtools_models_h
#define rxtools_models_h

#include "allometry.hpp"

namespace rxtools::allometry {
    namespace transforms {
        const UnivariateLinearModel::Transform none;
        const UnivariateLinearModel::Transform square; //The response should be squared
        const UnivariateLinearModel::Transform cube; //The response should be cubed
        const UnivariateLinearModel::Transform power; //A log-log transform should be applied

        const UnivariateLinearModel::Transform sqrt; //The response should be square-rooted
        const UnivariateLinearModel::Transform curt; //The response should be cube-rooted
        const UnivariateLinearModel::Transform log; //The response should be log-transformed
    }

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
            Parameters() : slope(0), intercept(0), rsq(0), transform(transforms::none) {}
            Parameters(const double& s, const double& i, const Transform& t, const double& r) :
                slope(s), intercept(i), rsq(r), transform(t) {}
        };

        Parameters parameters;
        double predict(double x, const lapis::LinearUnit& thisUnit, const lapis::LinearUnit& returnUnit) const;


    protected:
        virtual std::vector<Transform> candidateTransforms() const = 0;
        void fitBestFromCandidateTransforms(const std::vector<double>& y, const std::vector<double>& x);
        Parameters calcModel(std::vector<double> y, std::vector<double> x, const Transform& tr) const;
        virtual void print(std::ostream& os) const;

    };

    class DbhModel : public UnivariateLinearModel {
    public:
        DbhModel() {};
        DbhModel(const double& slope, const double& intercept, const Transform& transform, const double& rsq,
            const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& diameterUnit);
        DbhModel(const std::vector<double>& heights, const std::vector<double>& diameters,
            const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& diameterUnit,
            std::optional<Transform> transform = std::nullopt);
    protected:
        std::vector<Transform> candidateTransforms() const override;
    };

    class CrownModel : public UnivariateLinearModel {
    public:
        CrownModel() {};
        CrownModel(const double& slope, const double& intercept, const Transform& transform, const double& rsq,
            const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& crownUnit);
        CrownModel(const std::vector<double>& heights, const std::vector<double>& crowns,
            const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& crownUnit,
            std::optional<Transform> transform = std::nullopt);
    protected:
        std::vector<Transform> candidateTransforms() const override;
    };

    class BiomassModel : public UnivariateLinearModel {
    public:
        BiomassModel() {};
        BiomassModel(const double& slope, const double& intercept, const Transform& transform, const double& rsq,
            const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& biomassUnit);
        BiomassModel(const std::vector<double>& heights, const std::vector<double>& biomasses,
            const lapis::LinearUnit& heightUnit, const lapis::LinearUnit& biomassUnit,
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