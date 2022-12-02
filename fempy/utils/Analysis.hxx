#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <vector>
#include <string>
#include <stdio.h>
#include <map>

#include "yaml-cpp/yaml.h"

class Analysis {
    public:
    Analysis(std::string analysisName, std::string ver) {
        name = analysisName;
        version = ver;
        if (name == "DstarPi") {
            YAML::Node config = YAML::LoadFile("/home/daniel/phsw/fempy/cfg_DstarPi_HMpp13TeV2022.yml");
            selections = config[version.data()]["selections"].as<std::map<std::string, std::string>>();
            aliases = config[version.data()]["aliases"].as<std::map<std::string, std::string>>();
            
        }
    }

    void Print() {
        std::cout << "name: " << name << "version: " << version << std::endl;

        std::cout << "aliases: " << std::endl;
        for (auto alias : aliases)
            std::cout << alias.first <<"\t  --->  " << alias.second << std::endl;
        std::cout << std::endl;
            
        std::cout << "selections: " << std::endl;
        for (auto sel : selections)
            std::cout << sel.first << "  :  " << sel.second << std::endl;
        std::cout << std::endl;
        
    }

    std::string name;
    std::string version;
    std::map<std::string, std::string> aliases;
    std::map<std::string, std::string> selections;
};

#endif