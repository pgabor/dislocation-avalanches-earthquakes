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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "simulation_data.h"
#include "argreader.h"
#include "precision_handler.h"

#include <sstream>


namespace sdddst {

class Simulation
{
public:
    Simulation(int argc, char **argv);
    ~Simulation();

    void integrate(const double & stepsize, std::vector<Dislocation> &newDislocation, const std::vector<Dislocation> &old, bool useSpeed2, bool calculateInitSpeed, double externalStress);
    void calculateSpeeds(const std::vector<Dislocation> & dis, std::vector<double>  & res, double externalStress);
    void calculateG(const double & stepsize, std::vector<Dislocation> &newDislocation, const std::vector<Dislocation> &old, bool useSpeed2, bool calculateInitSpeed, bool useInitSpeedForFirstStep, double externalStress);
    void calculateJacobian(const double &stepsize, const std::vector<Dislocation> &data);
    void calculateXError();

    void calculateSparseFormForJacobian();
    void solveEQSys();

    double calculateOrderParameter(const std::vector<double> & speeds);
    double calculateStrainIncrement(const std::vector<Dislocation> & old, const std::vector<Dislocation> & newD);

    double getElement(int j, int si, int ei);

    double getSimTime();

    void run();

    bool step();

    void stepStageI();
    void stepStageII();
    void stepStageIII();

    const std::vector<Dislocation> & getStoredDislocationData();

private:
    double simStarted;
    bool succesfulStep = true;
    double implicitProbe = 0;
    double runResult = 0;
    double computationBegin = 0;
    double computationEnd = 0;
    double probeH = sD->stepSize;
    bool implicitProbeActive = false;
    bool multiplierProbeActive = false;
    double setMultiplier = sD->cutOffMultiplier;
    bool useBackupFile1 = false;
    double lastWriteTimeFinished;
    ArgReader args;
    std::shared_ptr<SimulationData> sD;
    std::string configStoragePath;
    bool saveConfigAtLog;
    bool initSpeedCalculationIsNeeded;
    std::unique_ptr<PrecisionHandler> pH;
    bool firstStepRequest;

    std::stringstream ss;

    bool avalancheInProgress;
    double energy;
    double energyAccum;
    double vsquare;
    double vsquare1;
    double vsquare2;

    double strainI;

    double calculateStress(double currentTime, double rate, double springConstant, double deformation);
};

}

#endif
