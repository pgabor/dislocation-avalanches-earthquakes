#include "simulation.h"

#include "utility.h"
#include "constants.h"

#include <iomanip>
#include <cmath>
#include <umfpack.h>
#include <algorithm>
#include <numeric>

using namespace sdddst;

Simulation::Simulation(int argc, char ** argv) :
    simStarted(get_wall_time()),
    succesfulStep(true),
    implicitProbe(0),
    runResult(0),
    computationBegin(0),
    computationEnd(0),
    probeH(0),
    implicitProbeActive(false),
    multiplierProbeActive(false),
    setMultiplier(0),
    useBackupFile1(false),
    lastWriteTimeFinished(0),
    args(argc, argv),
    sD(args.createSimulationData()),
    configStoragePath(""),
    saveConfigAtLog(false),
    initSpeedCalculationIsNeeded(true),
    pH(new PrecisionHandler),
    firstStepRequest(true),
    avalancheInProgress(false),
    energy(0),
    energyAccum(0),
    vsquare(0)
{
    probeH = sD->stepSize;
    setMultiplier = sD->cutOffMultiplier;
    args.clog << std::scientific << std::setprecision(16);
    args.speedLog << std::scientific << std::setprecision(16);
    if (args.is("WRITE_OUT_CONFIGS_WITH_LOG"))
    {
        saveConfigAtLog = true;
    }
    if (args.is("CONFIG_WRITE_DIRECTORY"))
    {
        configStoragePath = args.getString("CONFIG_WRITE_DIRECTORY");
    }

    // Precisity settings
    if (args.is("PREC"))
    {
        pH->setMinPrecisity(args.getDouble("PREC"));
    }
    pH->setSize(sD->dc);
    ss << std::setprecision(16) << std::scientific;
}

Simulation::~Simulation()
{
}


void Simulation::integrate(const double &stepsize, std::vector<Dislocation> &newDislocation, const std::vector<Dislocation> & old,
                           bool useSpeed2, bool calculateInitSpeed, double externalStress)
{
    calculateJacobian(stepsize, newDislocation);
    calculateSparseFormForJacobian();
    for (size_t i = 0; i < sD->ic; i++)
    {
        if (i > 0)
        {
            calculateG(stepsize, newDislocation, old, useSpeed2, false, false, externalStress);
        }
        else
        {
            calculateG(stepsize, newDislocation, old, useSpeed2, calculateInitSpeed, sD->externalStressType == STRESS_FREE ? true : false, externalStress);
        }
        solveEQSys();
        for (size_t j = 0; j < sD->dc; j++)
        {
            newDislocation[j].x -= sD->x[j];
        }
    }
    umfpack_di_free_numeric (&sD->Numeric) ;
}

void Simulation::calculateSpeeds(const std::vector<Dislocation> &dis, std::vector<double> &res, double externalStress)
{
    std::fill(res.begin(), res.end(), 0);

    for (size_t i = 0; i < sD->dc; i++)
    {
        for (size_t j = i+1; j < sD->dc; j++)
        {
            double dx = dis[i].x - dis[j].x;
            normalize(dx);

            double dy = dis[i].y - dis[j].y;
            normalize(dy);

            double tmp = dis[i].b * dis[j].b * sD->tau->xy(dx, dy);

            double r2 = dx*dx+dy*dy;
            pH->updateTolerance(r2, i);
            pH->updateTolerance(r2, j);

            res[i] +=  tmp;
            res[j] -=  tmp;
        }

        for (size_t j = 0; j < sD->fc; j++)
        {
            double dx = dis[i].x - sD->fpoints[j].x;
            normalize(dx);

            double dy = dis[i].y - sD->fpoints[j].y;
            normalize(dy);

            double xSqr = X2(dx);
            double ySqr = X2(dy);
            double rSqr = xSqr + ySqr;
            double expXY = exp(-sD->KASQR * rSqr);
            res[i] -= 2.0 * sD->A * X(dx) * X(dy) * ((1.-expXY)/rSqr- sD->KASQR * expXY) / rSqr * dis[i].b;

            pH->updateTolerance(rSqr, i);
        }

        res[i] += dis[i].b * externalStress;
    }
}

void Simulation::calculateG(const double &stepsize, std::vector<Dislocation> &newDislocation, const std::vector<Dislocation> &old,
                            bool useSpeed2, bool calculateInitSpeed, bool useInitSpeedForFirstStep, double externalStress)
{
    std::vector<double> * isp = &(sD->initSpeed);
    std::vector<double> * csp = &(sD->speed);
    if (useSpeed2)
    {
        isp = &(sD->initSpeed2);
        csp = &(sD->speed2);
    }

    if (calculateInitSpeed)
    {
        calculateSpeeds(old, *isp, externalStress);
    }

    if (useInitSpeedForFirstStep)
    {
        csp = isp;
    }
    else
    {
        calculateSpeeds(newDislocation, *csp, externalStress);
    }

    for (size_t i = 0; i < sD->dc; i++)
    {
        sD->g[i] = newDislocation[i].x - (1.0+sD->dVec[i]) * 0.5 * stepsize * (*csp)[i] - old[i].x - (1.0 - sD->dVec[i]) * 0.5 * stepsize * (*isp)[i];
    }
}

double Simulation::getElement(int j, int si, int ei)
{
    int len = ei - si;
    if (len > 1)
    {
        int tmp = len /2;
        double a;

        if (sD->Ai[si+tmp] > j)
        {
            a = getElement(j, si, si + tmp);
            if (a != 0.0)
            {
                return a;
            }
        }
        else
        {
            a = getElement(j, si + tmp, ei);
            if (a != 0.0)
            {
                return a;
            }
        }
    }
    else
    {
        if (sD->Ai[si] == j)
        {
            return sD->Ax[si];
        }
        return 0;
    }
    return 0;
}

double Simulation::getSimTime()
{
    return sD->simTime;
}

void Simulation::calculateJacobian(const double & stepsize, const std::vector<Dislocation> & data)
{
    int totalElementCounter = 0;

    for (unsigned int j = 0; j < sD->dc; j++)
    {
        // Previously calculated part
        for (unsigned int i = 0; i < j; i++)
        {
            double v = getElement(j, sD->Ap[i], sD->Ap[i+1]);
            if (v != 0.0)
            {
                sD->Ai[totalElementCounter] = i;
                sD->Ax[totalElementCounter++] = v;
            }
        }
        // Add the diagonal element (it will be calculated later and the fpoints now)
        sD->Ai[totalElementCounter] = j;
        double tmp = 0;
        double dx;
        double dy;
        for (size_t l = 0; l < sD->fc; l++)
        {
            dx = data[j].x - sD->fpoints[l].x;
            normalize(dx);
            dy = data[j].y - sD->fpoints[l].y;
            normalize(dy);

            if (pow(sqrt(dx * dx + dy * dy) - sD->cutOff, 2) < 36.8 * sD->cutOffSqr)
            {
                double multiplier = 1;
                if (dx * dx + dy * dy > sD->cutOffSqr)
                {
                    multiplier = exp(-pow(sqrt(dx*dx+dy*dy)-sD->cutOff, 2) * sD->onePerCutOffSqr);
                }
                tmp -= data[j].b * (- sD->A * cos(0.2e1 * M_PI * dx) / M_PI * sin(0.2e1 * M_PI * dy) * ((0.1e1 - pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 +
                                                                                                                                        (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1))) /
                                                                                                        ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                                                                                                         pow(M_PI, -0.2e1) / 0.2e1) - sD->KASQR * pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) *
                                                                                                                                                                         pow(M_PI, -0.2e1) / 0.2e1 +
                                                                                                                                                                         (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                                                                                                                                                                         pow(M_PI, -0.2e1) / 0.2e1))) /
                                    ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1)
                                    - sD->A * sin(0.2e1 * M_PI * dx) * pow(M_PI, -0.2e1) * sin(0.2e1 * M_PI * dy) * (pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 +
                                                                                                                                            (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1)) *
                                                                                                                     sD->KASQR * sin(0.2e1 * M_PI * dx) / M_PI * log(M_E) / ((0.1e1 - cos(0.2e1 * M_PI * dx)) *
                                                                                                                                                                             pow(M_PI, -0.2e1) / 0.2e1 +
                                                                                                                                                                             (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                                                                                                                                                                             pow(M_PI, -0.2e1) / 0.2e1) -
                                                                                                                     (0.1e1 - pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) /
                                                                                                                                                     0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                                                                                                                                                     pow(M_PI, -0.2e1) / 0.2e1))) *
                                                                                                                     pow((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 +
                                                                                                                         (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1, -0.2e1) *
                                                                                                                     sin(0.2e1 * M_PI * dx) / M_PI + sD->KASQR * sD->KASQR * pow(M_E, -sD->KASQR *
                                                                                                                                                                                 ((0.1e1 - cos(0.2e1 * M_PI * dx)) *
                                                                                                                                                                                  pow(M_PI, -0.2e1) /
                                                                                                                                                                                  0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                                                                                                                                                                                  pow(M_PI, -0.2e1) / 0.2e1)) *
                                                                                                                     sin(0.2e1 * M_PI * dx) / M_PI * log(M_E)) / ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) /
                                                                                                                                                                  0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                                                                                                                                                                  pow(M_PI, -0.2e1) / 0.2e1) / 0.2e1 + sD->A *
                                    pow(sin(0.2e1 * M_PI * dx), 0.2e1) * pow(M_PI, -0.3e1) * sin(0.2e1 * M_PI * dy) * ((0.1e1 - pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 +
                                                                                                                                                       (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1))) /
                                                                                                                       ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                                                                                                                        pow(M_PI, -0.2e1) / 0.2e1) - sD->KASQR * pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) *
                                                                                                                                                                                        pow(M_PI, -0.2e1) / 0.2e1 +
                                                                                                                                                                                        (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                                                                                                                                                                                        pow(M_PI, -0.2e1) / 0.2e1))) *
                                    pow((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1, -0.2e1) / 0.2e1) *  multiplier;
            }
        }
        sD->Ax[totalElementCounter++] = - tmp * stepsize;
        // Totally new part
        for (unsigned int i = j+1; i < sD->dc; i++)
        {
            dx = data[i].x - data[j].x;
            normalize(dx);

            dy = data[i].y - data[j].y;
            normalize(dy);

            if (pow(sqrt(dx * dx + dy * dy) - sD->cutOff, 2) < 36.8 * sD->cutOffSqr)
            {
                double multiplier = 1;
                if (dx * dx + dy * dy > sD->cutOffSqr)
                {
                    multiplier = exp(-pow(sqrt(dx*dx+dy*dy)-sD->cutOff, 2) * sD->onePerCutOffSqr);
                }
                sD->Ai[totalElementCounter] = i;
                sD->Ax[totalElementCounter++] = stepsize * data[i].b * data[j].b * sD->tau->xy_diff_x(dx, dy) * multiplier;
            }
        }
        sD->Ap[j+1] = totalElementCounter;
    }

    for (unsigned int j = 0; j < sD->dc; j++)
    {
        double subSum = 0;
        for (int i = sD->Ap[j]; i < sD->Ap[j+1]; i++)
        {
            if (sD->Ai[i] == int(j))
            {
                sD->indexes[j] = i;
            }
            subSum += sD->Ax[i];
        }

        subSum *= -1.;
        sD->Ax[sD->indexes[j]] = subSum;
        if (subSum > 0)
        {
            subSum = 1./subSum;
            subSum += 1.;
            subSum *= subSum;
            sD->dVec[j] = 1./subSum;
        }
        else
        {
            sD->dVec[j] = 0.;
        }
    }

    for (unsigned int j = 0; j < sD->dc; j++)
    {
        for (int i = sD->Ap[j]; i < sD->Ap[j+1]; i++)
        {
            sD->Ax[i] *= (1.0+sD->dVec[sD->Ai[i]]) * 0.5;
        }
        sD->Ax[sD->indexes[j]] += 1.0;
    }
}

void Simulation::calculateSparseFormForJacobian()
{
    (void) umfpack_di_symbolic (sD->dc, sD->dc, sD->Ap, sD->Ai, sD->Ax, &(sD->Symbolic), sD->null, sD->null);
    (void) umfpack_di_numeric (sD->Ap, sD->Ai, sD->Ax, sD->Symbolic, &(sD->Numeric), sD->null, sD->null);
    umfpack_di_free_symbolic (&(sD->Symbolic));
}

void Simulation::solveEQSys()
{
    (void) umfpack_di_solve (UMFPACK_A, sD->Ap, sD->Ai, sD->Ax, sD->x, sD->g.data(), sD->Numeric, sD->null, sD->null) ;
}

void Simulation::calculateXError()
{
    for (size_t i = 0; i < sD->dc; i++)
    {
        double tmp = fabs(sD->bigStep[i].x - sD->secondSmall[i].x);
        pH->updateError(tmp, i);
    }
}

double Simulation::calculateOrderParameter(const std::vector<double> &speeds)
{
    double orderParameter = 0;
    for (size_t i = 0; i < sD->dc; i++)
    {
        orderParameter += sD->dislocations[i].b * speeds[i];
    }
    return orderParameter;
}

double Simulation::calculateStrainIncrement(const std::vector<Dislocation> &old, const std::vector<Dislocation> &newD)
{
    double result = 0;
    for (size_t i = 0; i < old.size(); i++)
    {
        result += newD[i].b * (newD[i].x - old[i].x);
    }
    return result;
}

void Simulation::run()
{
    while (sD->isStopAtTimeSet == false || (sD->simTime < sD->stopAtTime && sD->isStopAtTimeSet))
    {
        step();
    }

    sD->writeDislocationDataToFile(args.resultPath);
}

bool Simulation::step()
{
    stepStageI();
    stepStageII();
    stepStageIII();
    return succesfulStep;
}

void Simulation::stepStageI()
{
    energyAccum = 0;
    if (!sD->switched && sD->switchTime < sD->simTime)
    {
        sD->switched = true;
        sD->cutOffMultiplier = sD->switchMultiplier;
        sD->updateCutOff();
    }

    double sumAvgSp = 0;
    vsquare = 0;
    if (firstStepRequest)
    {
        lastWriteTimeFinished = get_wall_time();
        calculateSpeeds(sD->dislocations, sD->initSpeed, calculateStress(sD->simTime, sD->stressRate, sD->springConstant, sD->accumulatedStrain));
        initSpeedCalculationIsNeeded = false;
        sumAvgSp = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0.0, [](double a, double b){return a + fabs(b);}) / double(sD->dc);
        vsquare = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0.0, [](double a, double b){return a + b*b;});

        // First log line
        args.clog << sD->simTime << " " << sD->succesfulSteps << " " << sD->failedSteps << " " << 0 << " " << sumAvgSp << " " << sD->cutOff << " " << "-";
        args.clog << " " << calculateStress(sD->simTime, sD->stressRate, sD->springConstant, sD->accumulatedStrain);

        args.clog << " " << "-";

        if (sD->calculateStrainDuringSimulation)
        {
            args.clog << " " << sD->accumulatedStrain;
        }
        else
        {
            args.clog << " -";
        }
        args.clog << " " << vsquare << " " << energy;

        args.clog << "\n";

        if (saveConfigAtLog)
        {
            sD->writeDislocationDataToFile(configStoragePath + "/" + std::to_string(sD->simTime) + ".dconf");
        }


        firstStepRequest = false;
    }

    // Reset the variables for the integration
    sD->bigStep = sD->dislocations;

    // Handle the adiabatic stress protocol
    if (sD->externalStressType == ADIABATIC)
    {
        sD->adiabaticStress.createNewBasedOn(0);
    }

    /////////////////////////////////
    /// Integrating procedure begins

    integrate(sD->stepSize, sD->bigStep, sD->dislocations, false, initSpeedCalculationIsNeeded, calculateStress(sD->simTime, sD->stressRate, sD->springConstant, sD->accumulatedStrain));

    // This can not get before the first integration step
    succesfulStep = false;
}

void Simulation::stepStageII()
{

    sD->firstSmall = sD->dislocations;

    integrate(0.5*sD->stepSize, sD->firstSmall, sD->dislocations, false, false, calculateStress(sD->simTime, sD->stressRate, sD->springConstant, sD->accumulatedStrain));

    strainI = calculateStrainIncrement(sD->dislocations, sD->firstSmall);
}

void Simulation::stepStageIII()
{
    double sumAvgSp = 0;
    sD->secondSmall = sD->firstSmall;

    integrate(0.5 * sD->stepSize, sD->secondSmall, sD->firstSmall, true, true, calculateStress(sD->simTime+0.5*sD->stepSize, sD->stressRate, sD->springConstant, sD->accumulatedStrain + strainI));
    if (sD->externalStressType == ADIABATIC)
    {
            sD->adiabaticStress.addStrainIncrement(2, calculateStrainIncrement(sD->firstSmall, sD->secondSmall));
    }

    vsquare1 = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0.0, [](double a, double b){return a + b*b;});
    vsquare2 = std::accumulate(sD->initSpeed2.begin(), sD->initSpeed2.end(), 0.0, [](double a, double b){return a + b*b;});

    energyAccum = (vsquare1+vsquare2) * 0.5 * sD->stepSize * 0.5;

    calculateXError();

    /// Precision related error handling
    if (pH->getMaxErrorRatioSqr() < 1.0)
    {
        succesfulStep = true;
        initSpeedCalculationIsNeeded = true;

        sD->accumulatedStrain += strainI + calculateStrainIncrement(sD->firstSmall, sD->secondSmall);

        sD->dislocations.swap(sD->secondSmall);
        for (size_t i = 0; i < sD->dc; i++)
        {
            normalize(sD->dislocations[i].x);
        }
        sD->simTime += sD->stepSize;
        sD->succesfulSteps++;

        bool writeConfigToFileIfSuccess = false;

        if (sD->externalStressType == ADIABATIC)
        {
            writeConfigToFileIfSuccess = sD->adiabaticStress.triggered(2);
            if (sD->adiabaticStress.avalanche(2))
            {
                sD->avalancheCounter++;
                sD->totalAccumulatedStrainIncrease += sD->adiabaticStress.getStrainIncrement(2);
                args.avalancheLog << sD->avalancheCounter << " " << sD->adiabaticStress.getTriggerTime(2) << " " <<
                                     sD->adiabaticStress.getEndTriggerTime(2) << " " << sD->adiabaticStress.getStrainIncrement(2) << " " <<
                                     sD->adiabaticStress.getBaseStress(2) << " " << sD->totalAccumulatedStrainIncrease << std::endl;
                avalancheInProgress = false;
            }
            sD->adiabaticStress.deleteCase(0);
            sD->adiabaticStress.deleteCase(0);
        }

        double orderParameter = 0;
        if (args.is("CALCULATE_ORDER_PARAMETER"))
        {
            orderParameter = calculateOrderParameter(sD->speed);
        }


        double current_wall_time = get_wall_time();

        calculateSpeeds(sD->dislocations, sD->initSpeed, calculateStress(sD->simTime, sD->stressRate, sD->springConstant, sD->accumulatedStrain));
        initSpeedCalculationIsNeeded = false;
        sumAvgSp = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0.0, [](double a, double b){return a + fabs(b);}) / double(sD->dc);
        vsquare = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0.0, [](double a, double b){return a + b*b;});

        energyAccum += (vsquare2+vsquare)* 0.5 * sD->stepSize * 0.5;
        args.clog << sD->simTime << " " << sD->succesfulSteps << " " << sD->failedSteps << " " << pH->getMaxErrorRatioSqr() << " " << sumAvgSp << " " << sD->cutOff << " ";

        if (!sD->switched && sumAvgSp < 0.5)
        {
            sD->switched = true;
            sD->cutOffMultiplier = sD->switchMultiplier;
            sD->updateCutOff();
        }

        if (sD->orderParameterCalculationIsOn)
        {
            args.clog << orderParameter;
        }
        else
        {
            args.clog << "-";
        }

        args.clog << " " << calculateStress(sD->simTime, sD->stressRate, sD->springConstant, sD->accumulatedStrain);

        args.clog << " " << current_wall_time - lastWriteTimeFinished;

        if (sD->calculateStrainDuringSimulation)
        {
            args.clog << " " << sD->accumulatedStrain;
        }
        else
        {
            args.clog << " -";
        }

        energy += energyAccum;

        args.clog << " " << vsquare << " " << energy;

        args.clog << "\n";


        if (saveConfigAtLog && sD->succesfulSteps % 5000 == 0)
        {
            ss.str("");
            ss << sD->succesfulSteps;
            sD->writeDislocationDataToFile(configStoragePath + "/" + ss.str() + ".dconf");
        }

        if (sD->saveSpeedsDuringAvalanche && writeConfigToFileIfSuccess)
        {
            args.speedLog << sD->simTime;
            for (auto & a: sD->initSpeed)
            {
                args.speedLog << " " << a;
            }
            args.speedLog << "\n";
        }

        lastWriteTimeFinished = get_wall_time();

    }
    else
    {
        sD->failedSteps++;
    }

    sD->stepSize = pH->getNewStepSize(sD->stepSize);
    pH->reset();

    if (sD->isMaxStepSizeLimit && sD->maxStepSizeLimit < sD->stepSize)
    {
        sD->stepSize = sD->maxStepSizeLimit;
    }

    if (!succesfulStep && sD->externalStressType == ADIABATIC)
    {
        sD->adiabaticStress.deleteCase(2);
        sD->adiabaticStress.deleteCase(1);
    }
}

const std::vector<Dislocation> &Simulation::getStoredDislocationData()
{
    return sD->dislocations;
}

double Simulation::calculateStress(double currentTime, double rate, double springConstant, double deformation)
{
    return currentTime * rate - deformation * springConstant;
}
