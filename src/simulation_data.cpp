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

#include "simulation_data.h"
#include "constants.h"

#include <fstream>
#include <cassert>
#include <iomanip>
#include <sstream>

SimulationData::SimulationData(const std::string &dislocationDataFilePath, const std::string &fixpointsDataFilePath):
    cutOffMultiplier(DEFAULT_CUTOFF_MULTIPLIER),
    cutOff(DEFAULT_CUTOFF),
    cutOffSqr(cutOff*cutOff),
    onePerCutOffSqr(1./cutOffSqr),
    prec(DEFAULT_PRECISION),
    fc(0),
    dc(0),
    ic(DEFAULT_ITERATION_COUNT),
    timeLimit(DEFAULT_TIME_LIMIT),
    stepSize(DEFAULT_STEP_SIZE),
    simTime(DEFAULT_SIM_TIME),
    computerPrecision(EPS),
    KASQR(DEFAULT_KASQR),
    A(DEFAULT_A),
    tau(nullptr),
    exV(DEFAULT_EXTERNAL_FIELD),
    Ap(nullptr),
    Ai(nullptr),
    Ax(nullptr),
    x(nullptr),
    null(NULL),
    Symbolic(nullptr),
    Numeric(nullptr),
    succesfulSteps(0),
    failedSteps(0),
    externalStressType(STRESS_FREE),
    avalancheCounter(0),
    isAvalancheLimit(false),
    avalancheLimit(0),
    totalAccumulatedStrainIncrease(0),
    isStrainIncreaseLimit(false),
    totalAccumulatedStrainIncreaseLimit(0),
    isMaxStepSizeLimit(false),
    maxStepSizeLimit(0),
    switched(true),
    switchTime(0),
    switchMultiplier(cutOffMultiplier),
    saveSpeedsDuringAvalanche(false),
    accumulatedStrain(0),
    calculateStrainDuringSimulation(false),
    orderParameterCalculationIsOn(false),
    triggerLimitIsOn(false),
    triggerLimitValue(0),
    stressRate(0),
    springConstant(1.0),
    stopAtTime(0),
    isStopAtTimeSet(false)
{
    readDislocationDataFromFile(dislocationDataFilePath);
    readFixpointDataFromFile(fixpointsDataFilePath);
    initSimulationVariables();
}

void SimulationData::readDislocationDataFromFile(const std::string &dislocationDataFilePath)
{
    if (dislocationDataFilePath.empty())
    {
        return;
    }
#ifdef RESPONSE
    std::cerr << "Dislocation data read from: " << dislocationDataFilePath << "\n";
#endif
    //Cleaning up
    deleteDislocationCountRelatedData();

    std::ifstream in(dislocationDataFilePath);
    assert(in.is_open() && "Cannot open dislocation data file!");

    // Iterating through the file
    while(!in.eof())
    {
        std::string data;
        in >> data;
        if (data == "")
        {
            break;
        }
        Dislocation tmp;
        tmp.x = std::stod(data);
        in >> tmp.y;
        in >> tmp.b;
        dislocations.push_back(std::move(tmp));
        dc++;
    }
    updateMemoryUsageAccordingToDislocationCount();
#ifdef RESPONSE
    std::cerr << "Dislocation data read ended.\n";
#endif
}

void SimulationData::writeDislocationDataToFile(const std::string &dislocationDataFilePath)
{
#ifdef RESPONSE
    std::cerr << "Dislocation data write to file: " << dislocationDataFilePath << "\n";
#endif
    std::ofstream out(dislocationDataFilePath);
    assert(out.is_open() && "Cannot open the data file to write!");
    out << std::scientific << std::setprecision(16);
    for (auto & i: dislocations)
    {
        out << i.x << " " << i.y << " " << i.b << "\n";
    }
#ifdef RESPONSE
    std::cerr << "Dislocation data write ended.\n";
#endif
}

void SimulationData::readFixpointDataFromFile(const std::string &fixpointDataFilePath)
{
    if (fixpointDataFilePath.empty())
    {
        fc = 0;
        return;
    }
#ifdef RESPONSE
    std::cerr << "Fixpoint data read from: " << fixpointDataFilePath << "\n";
#endif
    fc = 0;
    fpoints.resize(0);
    std::ifstream in(fixpointDataFilePath);
    assert(in.is_open() && "Cannot open fixpoints data file!");

    // Iterating through the file
    while (!in.eof())
    {
        std::string data;
        in >> data;
        if (data == "")
        {
            break;
        }
        FPoint tmp;
        tmp.x = std::stod(data);
        in >> tmp.y;
        fpoints.push_back(std::move(tmp));
        fc++;
    }
#ifdef RESPONSE
    std::cerr << "Fixpoint data read ended.\n";
#endif
}

void SimulationData::writeFixpointDataToFile(const std::string &fixpointDataFilePath)
{
#ifdef RESPONSE
    std::cerr << "Fixpoint data write to file: " << fixpointDataFilePath << "\n";
#endif
    std::ofstream out(fixpointDataFilePath);
    assert(out.is_open() && "Cannot open the data file to write!");
    out << std::scientific << std::setprecision(16);
    for (auto & i: fpoints)
    {
        out << i.x << " " << i.y << "\n";
    }
#ifdef RESPONSE
    std::cerr << "Fixpoint data write ended.\n";
#endif
}

void SimulationData::initSimulationVariables()
{
    updateCutOff();
}

void SimulationData::updateCutOff()
{
    double multiplier = cutOffMultiplier;
    if (cutOffMultiplier >=  1./12. * sqrt(2.0 * double(dc)))
    {
        multiplier = 1e18;
    }
    cutOff = 1./sqrt(dc) * multiplier;
    cutOffSqr = cutOff * cutOff;
    onePerCutOffSqr = 1./cutOffSqr;
}

void SimulationData::deleteDislocationCountRelatedData()
{
    dc = 0;
    dislocations.resize(0);
    g.resize(0);
    initSpeed.resize(0);
    initSpeed2.resize(0);
    speed.resize(0);
    speed2.resize(0);
    dVec.resize(0);
    bigStep.resize(0);
    firstSmall.resize(0);
    secondSmall.resize(0);
    free(Ap);
    Ap = nullptr;
    free(Ai);
    Ai = nullptr;
    free(Ax);
    Ax = nullptr;
    free(x);
    x = nullptr;
    indexes.resize(0);
}

void SimulationData::updateMemoryUsageAccordingToDislocationCount()
{
    g.resize(dc);
    initSpeed.resize(dc);
    initSpeed2.resize(dc);
    speed.resize(dc);
    speed2.resize(dc);
    dVec.resize(dc);
    bigStep.resize(dc);
    firstSmall.resize(dc);
    secondSmall.resize(dc);
    Ap = (int*) calloc(dc+1, sizeof(int));
    Ai = (int*) calloc(dc*dc, sizeof(int));
    Ax = (double*) calloc(dc*dc, sizeof(double));
    x = (double*) calloc(dc, sizeof(double));
    assert(Ap && "Memory allication for Ap failed!");
    assert(Ai && "Memory allocation for Ai failed!");
    assert(Ax && "Memory allocation for Ax failed!");
    assert(x && "Memory allocation for x failed!");
    indexes.resize(dc);
}
