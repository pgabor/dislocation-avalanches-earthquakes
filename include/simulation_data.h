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

#ifndef SIMULATIONDATA_H
#define SIMULATIONDATA_H

#include "dislocation.h"
#include "fpoint.h"
#include "fields/field.h"
#include "fields/analytic_field.h"
#include "adiabatic_stress.h"

#include <vector>
#include <string>

enum ExternalStressType
{
    STRESS_FREE,
    CONSTANT,
    ADIABATIC
};

class SimulationData
{
public:
    /**
     * @brief SimulationData can be created with giving the file path for the dislocation and the fixpoint
     * data files
     * @param dislocationDataFilePath
     * @param fixpointsDataFilePath
     */
    SimulationData(const std::string & dislocationDataFilePath, const std::string & fixpointsDataFilePath);

    /// Data file read utilities
    void readDislocationDataFromFile(const std::string & dislocationDataFilePath);
    void writeDislocationDataToFile(const std::string & dislocationDataFilePath);
    void readFixpointDataFromFile(const std::string & fixpointDataFilePath);
    void writeFixpointDataToFile(const std::string & fixpointDataFilePath);

    /// Other utilities
    void initSimulationVariables();
    void updateCutOff();

    //////////////////
    /// DATA FIELDS
    ///

    std::vector<Dislocation> dislocations; //Valid dislocation position data -> state of the simulation at simTime
    std::vector<FPoint> fpoints; //The positions of the fix points
    std::vector<double> g;
    // Used to store initial speeds for the big step and for the first small step
    std::vector<double> initSpeed;
    // Used to store initial speeds for the second small step
    std::vector<double> initSpeed2;
    // Stores speed during the NR iteration for the big step and for the first small step
    std::vector<double> speed;
    // Stores speed during the NR iteration for the second small step
    std::vector<double> speed2;
    // Stores the d values for the integration scheme
    std::vector<double> dVec;

    // The value of the cut off multiplier
    double cutOffMultiplier;

    // The value of the cutoff
    double cutOff;

    // The value of the cutoff^2
    double cutOffSqr;

    // The value of the 1/(cutoff^2)
    double onePerCutOffSqr;

    // Precisity of the simulation
    double prec;

    // Count of the fixpoints in the system
    unsigned int fc;

    // Count of the dislocations in the system
    unsigned int dc;

    // Count of the iterations during the NR
    unsigned int ic;

    // Simulation time limit. After it is reached there should be no more calculations
    double timeLimit;

    // The current step size of the simulation
    double stepSize;

    // The current time in the simulation
    double simTime;

    // When comparing data two value will be equal if they differ less than this value
    double computerPrecision;

    // Scaling factor for fixpoint interaction strength calculation
    double KASQR;

    // Interaction strength between a fixpoint and a dislocation
    double A;

    // The dislocation data after the big step
    std::vector<Dislocation> bigStep;
    // The dislocation data after the first small step
    std::vector<Dislocation> firstSmall;
    // The dislocation data after the second small step
    std::vector<Dislocation> secondSmall;

    // The used field
    std::unique_ptr<sdddst::Field> tau;

    double exV; /// External constant stress

    // UMFPack specified sparse format stored Jacobian
    int * Ap;
    int * Ai;
    double * Ax;
    // Result data
    double * x;

    // UMFPack required variables
    double *null;
    void *Symbolic, *Numeric;

    // Diagonal indexes in the Jacobian
    std::vector<int> indexes;

    // Number of the successfuly finished steps
    size_t succesfulSteps;

    // Number of the unsuccessfuly finished steps
    size_t failedSteps;

    ////////////////////////// VERSION 0 /////////////////////////////

    ExternalStressType externalStressType;
    AdiabaticStress adiabaticStress;

    ////////////////////////// VERSION 1 /////////////////////////////

    unsigned int avalancheCounter;
    bool isAvalancheLimit;
    unsigned int avalancheLimit;
    double totalAccumulatedStrainIncrease;
    bool isStrainIncreaseLimit;
    double totalAccumulatedStrainIncreaseLimit;
    bool isMaxStepSizeLimit;
    double maxStepSizeLimit;

    bool switched;
    double switchTime;
    double switchMultiplier;

    ////////////////////////// VERSION 2 /////////////////////////////

    std::string fieldBinaryPath;

    bool saveSpeedsDuringAvalanche;

    double accumulatedStrain;
    bool calculateStrainDuringSimulation;
    bool orderParameterCalculationIsOn;
    bool triggerLimitIsOn;
    double triggerLimitValue;

    double stressRate;
    double springConstant;

    double stopAtTime;
    bool isStopAtTimeSet;

private:
    /**
     * @brief deleteDislocationCountRelatedData repsonsible to delete everything which is related to
     * the dislocations count and to reset everything to initial state regarding this
     */
    void deleteDislocationCountRelatedData();
    void updateMemoryUsageAccordingToDislocationCount();
};

#endif // SIMULATIONDATA_H
