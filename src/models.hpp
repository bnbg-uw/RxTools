#pragma once

#ifndef rxtools_models_h
#define rxtools_models_h

#include "allometry.hpp"

namespace rxtools::allometry {

    // This model will create a linear model (perhaps with transform) from height to the specified response variable.
    class UnivariateLinearModel : public Model {
    public:
        enum class Transform {
            //Transform::None is relied upon to be ethe first element for looping purposes.
            None, //The response should not be transformed

            //Add new transforms here:
            Square, //The response should be squared
            Cube, //The response should be cubed
            Power, //A log-log transform should be applied

            //Transform::Suggest is relied upon to be the last element for looping purposes.
            Suggest //All of the above should be tried, and the model with the highest R^2 should be used
        };

        struct Parameters {
            double slope;
            double intercept;
            Transform transform;
            double rsq;
            Parameters() : slope(0), intercept(0), transform(Transform::None), rsq(0) {}
            Parameters(const double& s, const double& i, const Transform& t, const double& r) :
                slope(s), intercept(i), transform(t), rsq(r) {}
        };

        Parameters parameters;

        UnivariateLinearModel() {};
        UnivariateLinearModel(const double& slope, const double& intercept, const Transform& transform, const double& rsq,
            const lapis::LinearUnit& inputUnit, const lapis::LinearUnit& outputUnit);
        UnivariateLinearModel(const FIATreeList& treeList, const std::string& responseName, const lapis::LinearUnit& responseUnit, const Transform& transform = Transform::Suggest);

        double predict(double x, const lapis::LinearUnit& thisUnit, const lapis::LinearUnit& returnUnit = lapis::linearUnitPresets::meter) const;


    protected:
        void print(std::ostream& os) const;
    private:
        Parameters calcModel(const std::vector<double>& y, const std::vector<double>& x, const Transform& tr) const;
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
        UnivariateLinearModel cr;
    };

} // namespace rxtools::allometry


#endif // !rxtools_models_h