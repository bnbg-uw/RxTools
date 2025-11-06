#pragma once

#include "rxtools_pch.hpp"

namespace rxtools::utilities {
    inline int randomIdxFromWeights(std::vector<double>& weights, std::default_random_engine& dre) {
        double sumOfWeight = 0;
        for (int i = 0; i < weights.size(); i++) {
            sumOfWeight += weights[i];
        }
        std::uniform_real_distribution<double> urd(0., sumOfWeight);
        double v = urd(dre);

        for (int i = 0; i < weights.size(); i++) {
            if (v < weights[i])
                return i;
            v -= weights[i];    
        }
        return weights.size() - 1;
    }

    inline std::vector<std::string> readCSVLine(std::istream& file) {
        char split = ',';
        std::string line = std::string();
        std::getline(file, line);
        std::istringstream i{ line };
        std::vector<std::string> out{};
        while (!i.eof()) {
            std::string token = std::string();
            std::getline(i, token, split);
            out.push_back(token);
        }
        return out;
    }
} // namepsace rxtools::utilities