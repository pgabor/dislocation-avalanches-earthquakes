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

#include "argreader.h"
#include "fields/field.h"
#include "fields/analytic_field.h"

#include <cassert>
#include <iostream>

ArgReader::ArgReader(const int &argc, char **argv):
    clog(NULL),
    avalancheLog(NULL)
{
    for ( int i = 1; i < argc; i++)
    {
        std::string option(argv[i]);
        size_t pos = option.find("=");
        if (pos == std::string::npos)
        {
            options[option] = "";
        }
        else
        {
            options[option.substr(0, pos)] = option.substr(pos+1,option.length()-pos-1);
        }
    }

    if (argc == 1)
    {
        printHelp();
        exit(0);
    }
}

ArgReader::~ArgReader()
{
    clog.close();
    avalancheLog.close();
}

void ArgReader::printHelp()
{
    std::cerr << "DCONF=FILE_PATH\t Path to the dislocation configuration file\n";
    std::cerr << "FCONF=FILE_PATH\t Path to the fix points configuration file\n";
    std::cerr << "CLOG=FILE_PATH\t PATH to write the continous log\n";
    std::cerr << "STEPSIZE=VALUE\t Size of a step\n";
    std::cerr << "PREC=VALUE\t Simulation precisity (at least)\n";
    std::cerr << "CUTOFF_MULTIPLIER=VALUE\t The cutoff should be multiplied by this value\n";
    std::cerr << "RES_DCONF=PATH\t Path where the result dislocation conf has to be written.\n";
    std::cerr << "STEPSIZE_LIMIT=VALUE\t Maximal step size during the simulation.\n";
    std::cerr << "CONFIG_WRITE_DIRECTORY=Path\t Relative path where the configurations should be stored.\n";
    std::cerr << "WRITE_OUT_CONFIGS_WITH_LOG\t If this is set, the configuration will be written out when a log entry is created.\n";
    std::cerr << "CALCULATE_ORDER_PARAMETER\t If this is set, it will calculate the order parameter.\n";
    std::cerr << "SWITCH_TIME\t When should change to alternative cutoff.\n";
    std::cerr << "SWITCH_VALUE\t What should be the new cutoff after the switch.\n";
    std::cerr << "SAVE_SPEEDS_DURING_AVALANCHE=PATH\t Log the speed of each invidual dislocation into the specified logfile\n";
    std::cerr << "CALCULATE_STRAIN_CONTINOUSLY\t Calculates the strain during the whole simulation\n";
    std::cerr << "TRIGGER_TIME_LIMIT=VALUE\t If an avalanche is longer than this limit, the simulation will be finished.\n";
    std::cerr << "STRESS_RATE=Value\t rate for time related stress.\n";
}

bool ArgReader::is(const std::string &option)
{
    return options.count(option);
}

int ArgReader::getInt(const std::string &option)
{
    return std::stoi(options.at(option));
}

double ArgReader::getDouble(const std::string &option)
{
    return std::stod(options.at(option));
}

std::string ArgReader::getString(const std::string &option)
{
    return options.at(option);
}

std::shared_ptr<SimulationData> ArgReader::createSimulationData()
{
    std::shared_ptr<SimulationData> sD(nullptr);
    if (is("DCONF"))
    {
        sD = std::shared_ptr<SimulationData>(new SimulationData(getString("DCONF"), is("FCONF") ? getString("FCONF") : ""));
    }
    else
    {
        assert(false && "No initial data was given!");
    }

    if (is("CLOG"))
    {
        clog.close();
        clog.open(getString("CLOG"), std::fstream::app);
    }

    if (is("SAVE_SPEEDS_DURING_AVALANCHE"))
    {
        speedLog.close();
        speedLog.open(getString("SAVE_SPEEDS_DURING_AVALANCHE"), std::fstream::app);
        sD->saveSpeedsDuringAvalanche = true;
    }

    if (is("STRESS_RATE"))
    {
        sD->stressRate = getDouble("STRESS_RATE");
    }

    if (is("STEPSIZE"))
    {
        sD->stepSize = getDouble("STEPSIZE");
    }

    if (is("PREC"))
    {
        sD->prec = getDouble("PREC");
    }

    if (is("CUTOFF_MULTIPLIER"))
    {
        sD->cutOffMultiplier = getDouble("CUTOFF_MULTIPLIER");
        sD->updateCutOff();
    }

    if (is("RES_DCONF"))
    {
        resultPath = getString("RES_DCONF");
    }

    if (is("START_TIME"))
    {
      sD->simTime = getDouble("START_TIME");
    }

    if (is("START_STRAIN"))
    {
      sD->accumulatedStrain = getDouble("START_STRAIN");
    }

    if (is("AVALANCHE_LOG"))
    {
        avalancheLog.close();
        avalancheLog.open(getString("AVALANCHE_LOG"), std::fstream::app);
    }

    if (is("AVALANCHE_LIMIT"))
    {
        sD->isAvalancheLimit = true;
        sD->avalancheLimit = getDouble("AVALANCHE_LIMIT");
    }

    if (is("AVALANCHE_STRAIN_INCREASE_LIMIT"))
    {
        sD->isStrainIncreaseLimit = true;
        sD->totalAccumulatedStrainIncreaseLimit = getDouble("AVALANCHE_STRAIN_INCREASE_LIMIT");
    }

    if (is("STEPSIZE_LIMIT"))
    {
        sD->isMaxStepSizeLimit = true;
        sD->maxStepSizeLimit = getDouble("STEPSIZE_LIMIT");
    }

    sD->tau = std::unique_ptr<sdddst::Field>(new sdddst::AnalyticField);

    if (is("SWITCH_TIME") && is("SWITCH_VALUE"))
    {
        sD->switchMultiplier = getDouble("SWITCH_VALUE");
        sD->switchTime = getDouble("SWITCH_TIME");
        sD->switched = false;
    }

    if (is("FCONF"))
    {
        sD->A *= 1./sqrt(sD->dc);
        sD->KASQR *= double(sD->dc);
    }

    if (is("STOP_AT_TIME"))
    {
	sD->stopAtTime = getDouble("STOP_AT_TIME");
	sD->isStopAtTimeSet = true;
    }

    if (is("CALCULATE_STRAIN_CONTINOUSLY"))
    {
        sD->calculateStrainDuringSimulation = true;
    }

    if(is("CALCULATE_ORDER_PARAMETER"))
    {
        sD->orderParameterCalculationIsOn = true;
    }

    if (is("TRIGGER_TIME_LIMIT"))
    {
        sD->triggerLimitIsOn = true;
        sD->triggerLimitValue = getDouble("TRIGGER_TIME_LIMIT");
    }

    return sD;
}
