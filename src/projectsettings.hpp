//Struct for holding options parsed from cmd line. 

//TODO: handle landscape indices.

#pragma once

#ifndef licosim_projectsettings_h
#define licosim_projectsettings_h

#include "rxtools_pch.hpp"

namespace licosim {
    struct ProjectSettings {
        std::string exeParentPath;
        std::string commandLine;
        std::string outputPath;
        bool writeUnits = false;
        bool fastFuels = false;
        std::string lidarDatasetPath;
        std::string projectPolygonPath;
        std::string unitPolygonPath;
        std::string lmuRasterPath;
        std::string priorityMethod;
        std::string priorityColumn;
        std::string terrain = "moderate";
        std::string referenceDatasetPath;
        std::string fireRefPath = "/resources/reference/fire.csv";
        std::string hydroRefPath = "/resources/reference/hydro.csv";
        std::string habitatRefPath = "/resources/reference/habitat.csv";
        std::string aetPath = "/resources/biophysicalanalogs/aet.img";
        std::string cwdPath = "/resources/biophysicalanalogs/cwd.img";
        std::string tmnPath = "/resources/biophysicalanalogs/tmn.img";
        std::string climateClassPath = "/resources/biophysicalanalogs/ClimateClasses.img";
        std::string fiaPath = "/resources/fia/";
        //std::string mcsPath = "/resources/treatment/mcs_prop.csv";
        std::string defaultRefPath;
        int nThread = 1;
        int seed = -1;
        double dbhMin = 15.24;
        double dbhMax = 53.34;
        std::vector<double> allom_coefficients;
        bool allomPower = true;

        bool useFireRef = true;
        bool useHydroRef = false;
        bool useHabitatRef = false;
        bool doFireModel = false;

        bool subdivideLmus = false;
        bool overrideTargets = false;

        //Treatment behavior indices:
        //double percentLandscape = std::nan("");
        //double totalLandscape = std::nan("");
        //double clump1 = 0.5;
        //double clump2 = 0.5;  
        std::vector<double> forestType = { 1, 1, 1 };

        ProjectSettings(const std::string exeParentPath) {
            this->exeParentPath = exeParentPath;

            fiaPath = exeParentPath + fiaPath;
            fireRefPath = exeParentPath + fireRefPath;
            hydroRefPath = exeParentPath + hydroRefPath;
            habitatRefPath = exeParentPath + habitatRefPath;
            defaultRefPath = fireRefPath;
            aetPath = exeParentPath + aetPath;
            cwdPath = exeParentPath + cwdPath;
            tmnPath = exeParentPath + tmnPath;
            climateClassPath = exeParentPath + climateClassPath;
            //mcsPath = exeParentPath + mcsPath;
        }

        //TODO: Implement landscape indices.
        ProjectSettings(const std::string exeParentPath, const cxxopts::ParseResult& options) {
            this->exeParentPath = exeParentPath;

            fiaPath = exeParentPath + fiaPath;
            fireRefPath = exeParentPath + fireRefPath;
            hydroRefPath = exeParentPath + hydroRefPath;
            habitatRefPath = exeParentPath + habitatRefPath;
            //mcsPath = exeParentPath + mcsPath;
            defaultRefPath = fireRefPath;

            lidarDatasetPath = options["lidar"].as<std::string>();
            if (options.count("output"))
                outputPath = options["output"].as<std::string>();
            else
                outputPath = lidarDatasetPath;
            if (std::filesystem::exists(outputPath)) {
                auto fs = std::filesystem::path(outputPath);
                fs = fs / "licosim";
                std::filesystem::create_directory(fs);
                outputPath = fs.string();
            }
            else {
                throw std::runtime_error("Output path does not exist.");
            }

            if (options.count("writeunits")) writeUnits = true;
            if (options.count("fastfuels")) fastFuels = true;

            projectPolygonPath = options["projectpoly"].as<std::string>();
            unitPolygonPath = options["unitpoly"].as<std::string>();

            lmuRasterPath = options.count("lmu") ? options["lmu"].as<std::string>() : "";

            priorityMethod = options["priority"].as<std::string>();
            priorityColumn = options.count("column") ? options["column"].as<std::string>() : "column";
            try { if (std::stoi(priorityColumn) < 0) throw std::invalid_argument("Invalid column index"); } catch (...) {}
            if (options.count("terrain")) terrain = options["terrain"].as<std::string>();
            referenceDatasetPath = options.count("reference") ? options["reference"].as<std::string>() : "";
            if (options.count("aet"))
                aetPath = options["aet"].as<std::string>();
            else
                aetPath = exeParentPath + aetPath;

            if (options.count("cwd"))
                cwdPath = options["cwd"].as<std::string>();
            else
                cwdPath = exeParentPath + cwdPath;
            if (options.count("janmin"))
                tmnPath = options["janmin"].as<std::string>();
            else
                tmnPath = exeParentPath + tmnPath;
            climateClassPath = exeParentPath + climateClassPath;


            allom_coefficients = options.count("allom") ? options["allom"].as<std::vector<double>>() : std::vector<double>{};
            nThread = options["thread"].as<int>();

            if (options.count("seed")) seed = options["seed"].as<int>();

            useFireRef = options.count("usefire") ? options["usefire"].as<bool>() : false;
            useHydroRef = options.count("usehydro") ? options["usehydro"].as<bool>() : false;
            useHabitatRef = options.count("usehabitat") ? options["usehabitat"].as<bool>() : false;
            doFireModel = options.count("usefire") ? options["units"].as<bool>() : false;

            //if(options.count("pland"))
            //    percentLandscape = options["pland"].as<double>();
            //if (options.count("tland") && !options.count("pland"))
            //    totalLandscape = options["tland"].as<double>();
            
            dbhMin = options.count("dbhmin") ? options["dbhmin"].as<double>() : dbhMin;
            dbhMax = options.count("dbhmax") ? options["dbhmax"].as<double>() : dbhMax;

            if (options.count("cover"))
                forestType = options["cover"].as<std::vector<double>>();
            if (forestType.size() != 3)
                throw std::invalid_argument("forest type vector should be of length 3.");
            
            double forestSum = 0;
            for (auto v : forestType)
                forestSum += v;
            for (int i = 0; i < forestType.size(); i++)
                forestType[i] /= forestSum;

            overrideTargets = options.count("overridetargets") ? true : false;
            subdivideLmus = options.count("subdivideLmus") ? true : false;

        }
    };
}  // namespace licosim

#endif  // !licosim_projectsettings_h