/*
Copyright (C) 2018-2021 Gábor Péterffy <peterffy95@gmail.com>

This file is part of sdddstEQ.

sdddstEQ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

sdddstEQ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with sdddstEQ. If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef ARGREADER_H
#define ARGREADER_H

#include "simulation_data.h"

#include <map>
#include <memory>
#include <fstream>

class ArgReader
{
public:
    ArgReader(const int & argc, char ** argv);
    ~ArgReader();

    void printHelp();

    bool is(const std::string & option);
    int getInt(const std::string & option);
    double getDouble(const std::string & option);
    std::string getString(const std::string & option);

    std::shared_ptr<SimulationData> createSimulationData();

    std::ofstream clog;
    std::ofstream avalancheLog;
    std::ofstream speedLog;
    std::string backupPath1;
    std::string backupPath2;
    std::string resultPath;

private:
    std::map<std::string, std::string> options;
};

#endif // ARGREADER_H
